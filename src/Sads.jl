module Sads

using Distributions
using KernelDensity: kde

abstract type CrossSection end

include("crosssections.jl")
include("kalman.jl")
include("gvf.jl")
include("ensemble.jl")

"""
Estimate discharge by assimilating observed water surface elevations along a river channel.

"""
function estimate(qwbm, H, W, x, rₚ, ri; nₚ=Uniform(0.02, 0.05), nsamples=1000, nens=100)
    zm, zs, dm, ds = priors(qwbm, H, W, x, nₚ, rₚ, nsamples, nens)
    Qa = assimilate(H, W, x, maximum(W, dims=2), maximum(H, dims=2),
                    LogNormal(log(qwbm/sqrt(1+dm^2)), log(1+dm^2)), nₚ, rₚ, Normal(zm, zs), nens, ri)
    Qa
end

"""
Assimilate SWOT observations for river reach.

"""
function assimilate(H, W, x, wbf, hbf, Qₚ, nₚ, rₚ, zₚ, nens, ri; ϵₒ=0.01, cv_thresh=1.5)
    Qa = zeros(length(ri)-1, size(H, 2))
    S = diff(H, dims=1) ./ diff(x)
    S = [S[1, :]'; S]
    ze = zeros(length(x), nens)
    Qe, ne, re, ze[1, :] = lhs_ensemble(nens, Qₚ, nₚ, rₚ, zₚ)
    for t in 1:size(H,2)
        he = gvf_ensemble!(H[:, t], W[:, t], S[:, t], x, hbf, wbf, Qe, ne, re, ze)
        i = findall(he[1, :] .> 0)
        h = ze .+ he .* ((re .+ 1) ./ re)'
        X = repeat(Qe[i]', outer=length(ri)-1)
        XA = h[:, i]
        d = H[:, t]
        E = rand(Normal(ϵₒ, 1e-6), length(d), length(i)) .* rand([-1, 1], length(d), length(i))
        A = letkf(X, d, XA, E, [[j] for j in 1:length(ri)-1], [collect(ri[j]:ri[j+1]-1) for j in 1:length(ri)-1], diagR=true)
        Qa[:, t] = mean(abs.(A), dims=2)  # could also use the absolute of the mean
    end
    Qa
end

"""
Estimate the prior probability distributions of bed elevation and discharge using rejection sampling.

"""
function priors(qwbm, H, W, x, nₚ, rₚ, nsamples, nens)
    zbnds, qbnds = prior_bounds(qwbm, nsamples, H, W, x, nₚ, rₚ)
    dm, ds = discharge_prior(qwbm, nens, nsamples, H, W, x, nₚ, rₚ, Normal(minimum(H[1, :]), 0.1))
    zm, zs = bed_elevation_prior(qwbm, nens, nsamples, H, W, x, nₚ, rₚ, zbnds, [dm, ds])
    zm, zs, dm, ds
end

"""
Estimate the prior probability distribution of discharge.

"""
function discharge_prior(qwbm, nens, nsamples, H, W, x, nₚ, rₚ, zₚ)
    @info "Estimating discharge prior PDF"
    S = diff(H, dims=1) ./ diff(x)
    S = [S[1, :]'; S]
    obs = [std(sample(H[end, :], nens)) for _ in 1:nsamples]
    Fobs = kde(obs)
    h = zeros(nsamples)
    δ = lhs_ensemble(nsamples, Uniform(0.5, 2.0))[1]
    t = sample(collect(1:size(H, 2)), nsamples)
    for s in 1:nsamples
        μ = log(qwbm / sqrt(1 + δ[s]^2))
        σ = log(1 + δ[s]^2)
        Qₚ = LogNormal(μ, σ)
        Qe, ne, re, ze = prior_ensemble(x, Qₚ, nₚ, rₚ, zₚ, nens)
        he = gvf_ensemble!(H[:, t[s]], W[:, t[s]], S[:, t[s]], x, maximum(H, dims=2), maximum(W, dims=2),
                          Qe, [mean(ne) for i in 1:nens], [mean(re) for i in 1:nens], ze)
        i = findall(he[1, :] .> 0)
        h[s] = std(he[end, i] .* (mean(re) .+ 1) ./ mean(re) .+ ze[end, i])
    end
    δ = δ[.!isnan.(h)]
    h = h[.!isnan.(h)]
    Fmod = kde(h)
    L = 1
    accepted = [rand(Uniform(0, L)) * pdf(Fmod, s) <= pdf(Fobs, s) for s in  h]
    if length(findall(accepted)) < 2
        δₘ, δₛ = δ[findmin(abs.(h .- mean(obs)))[2]], 0.01
    else
        δₘ, δₛ = mean(δ[accepted]), std(δ[accepted])
    end
    δₘ, δₛ
end

"""
Estimate the prior probability distribution of bed elevation.

"""
function bed_elevation_prior(qwbm, nens, nsamples, H, W, x, nₚ, rₚ, zbnds, dpars)
    @info "Estimating bed elevation prior PDF"
    S = diff(H, dims=1) ./ diff(x)
    S = [S[1, :]'; S]
    obs = [h for h in H[end, :]]
    Fobs = kde(obs)
    h = zeros(nsamples)
    zd, dd = lhs_ensemble(nsamples, Uniform(zbnds...), Normal(dpars...))
    for s in 1:nsamples
        zₚ = Normal(zd[s], 0.1)
        Qₚ = LogNormal(log(qwbm/sqrt(1+dd[s]^2)), log(1+dd[s]^2))
        Qe, ne, re, ze = prior_ensemble(x, Qₚ, nₚ, rₚ, zₚ, nens)
        he = gvf_ensemble!(mean(H, dims=2), mean(W, dims=2), mean(S, dims=2), x, maximum(H, dims=2),
                           maximum(W, dims=2), Qe, [mean(ne) for i in 1:nens], [mean(re) for i in 1:nens], ze)
        i = findall(he[1, :] .> 0)
        h[s] = mean(he[end, i] .* (mean(re) .+ 1) ./ mean(re) .+ ze[end, i])
    end
    zd = zd[.!isnan.(h)]
    dd = dd[.!isnan.(h)]
    h = h[.!isnan.(h)]
    Fmod = kde(h)
    L = 1
    accepted = [rand(Uniform(0, L)) * pdf(Fmod, s) <= pdf(Fobs, s) for s in h]
    if length(findall(accepted)) < 2
        zₘ, zₛ = zd[findmin(abs.(h .- mean(obs)))[2]], 0.1
    else
        zₘ, zₛ = mean(zd[accepted]), std(zd[accepted])
    end
    zₘ, zₛ
end

"""
Estimate the bounds of the prior distributions.

"""
function prior_bounds(qwbm, nsamples, H, W, x, nₚ, rₚ)
    @info "Estimating prior distribution bounds"
    S = diff(H, dims=1) ./ diff(x)
    S = [S[1, :]'; S]
    # We will use arbitrary values that are large enough to represent the uninformative priors
    zₚ = Uniform(minimum(H[1, :]) - 30, minimum(H[1, :]))
    Qₚ = Uniform(qwbm / 10, qwbm * 10)
    ze = zeros(length(x), nsamples)
    Qe, ze[1, :] = lhs_ensemble(nsamples, Qₚ, zₚ)
    re = rand(rₚ, nsamples)
    ne = rand(nₚ, nsamples)
    he = gvf_ensemble!(mean(H, dims=2), mean(W, dims=2), mean(S, dims=2), x, maximum(H, dims=2),
                      maximum(W, dims=2), Qe, [mean(ne) for i in 1:nsamples], [mean(re) for i in 1:nsamples], ze)
    i = findall(he[1, :] .> 0)
    h = he[end, i] .* (mean(re) .+ 1) ./ mean(re) .+ ze[end, i]
    j = i[(h .> minimum(H[end, :])) .& (h .< maximum(H[end, :]))]
    zbnds = [minimum(ze[1, j]), maximum(ze[1, j])]
    qbnds = [minimum(Qe[j]), maximum(Qe[j])]
    zbnds, qbnds
end

end
