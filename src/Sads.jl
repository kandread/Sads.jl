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
    zm, zs, dm, ds, qb = priors(qwbm, H, W, x, nₚ, rₚ, nsamples, nens)
    Qₚ = Truncated(LogNormal(log(qwbm/sqrt(1+dm^2)), log(1+dm^2)), qb[1], qb[2])
    Qa, A0, n = assimilate(H, W, x, maximum(W, dims=2), maximum(H, dims=2),
                           Qₚ, nₚ, rₚ, Normal(zm, 0.1), nens, ri)
    Qa, A0, n
end

"""
Estimate bed slope by assimilating observed water surface elevations.

"""
function bed_slope(S, H, W, x, hbf, wbf, Qₚ, nₚ, rₚ, zₚ, nens; ϵₒ=0.01)
    S0 = [mean(S[j, :]) > 0 ? mean(S[j, :]) : mean(S[j, :][S[j,:] .> 0]) for j in 1:size(S, 1)]
    S0[isnan.(S0)] .= minimum(S[S .> 0])
    Se = rand.(Normal(1., 0.2), length(x), nens)
    Se = [Se[i, j] .* S0[i] for i in 1:length(x), j in 1:nens]
    ze = zeros(length(x), nens)
    Qe, ne, re, ze[1, :] = lhs_ensemble(nens, Qₚ, nₚ, rₚ, zₚ)
    he = gvf_ensemble!(mean(H, dims=2), mean(W, dims=2), Se, x, hbf, wbf, Qe, ne, re, ze)
    i = findall(he[1, :] .> 0)
    h = ze .+ he .* ((re .+ 1) ./ re)'
    X = Se[:, i]
    XA = h[:, i]
    d = mean(H, dims=2)[:, 1]
    E = rand(Normal(ϵₒ, 1e-6), length(d), length(i)) .* rand([-1, 1], length(d), length(i))
    A = letkf(X, d, XA, E, [collect(1:length(x))], [collect(1:length(d))], diagR=true)
    mean(A, dims=2)[:, 1]
end

"""
Estimate parameters for discharge estimation, i.e. A₀ and n.

"""
function estimate_Q_params(H, W, S, ri, ze, re, ne, Qa)
    nr, nt = size(Qa)
    nens = length(ne)
    A0 = zeros(nr)
    n = zeros(nr)
    Wbf = maximum(W, dims=2)
    Hbf = maximum(H, dims=2)
    Hmin = H[:, findmin(H[1, :])[2]]
    Wmin = W[:, findmin(H[1, :])[2]]
    ybf = Hbf .- ze
    ymin = Hmin .- ze
    ybf[ybf .< 0] .= 0.0
    ymin[ymin .< 0] .= 0.0
    A0e = Wbf .* (ymin ./ ybf).^(1 ./ re') .* (re ./ (re .+ 1))' .* ymin
    # A0e = Wmin .* (re ./ (re .+ 1))' .* ymin
    A0e = [mean(A0e[ri[j]:ri[j+1], e]) for j in 1:nr, e in 1:nens]
    dA = (H .- Hmin) .* (W .+ Wmin) ./ 2
    dA = reshape([mean(dA[ri[j]:ri[j+1], t]) for j in 1:nr, t in 1:nt]', nr*nt, 1)
    Wr = reshape([mean(W[ri[j]:ri[j+1], t]) for j in 1:nr, t in 1:nt]', nr*nt, 1)
    Sr = reshape([mean(S[ri[j]:ri[j+1], t]) for j in 1:nr, t in 1:nt]', nr*nt, 1)
    A = A0e .+ dA
    A[A .< 0] .= 0.0
    Qe = (1 ./ ne') .* A.^(5/3) .* Wr.^(-2/3) .* Sr.^(1/2)
    X = zeros(nr*2, nens)
    X[1:2:end, :] = A0e
    X[2:2:end, :] = repeat(ne', outer=nr)
    XA = Qe
    d = reshape(Qa, nr*nt)
    E = rand(Normal(0.1*mean(Qa), 1e-6), length(d), nens) .* rand([-1, 1], length(d), nens)
    A = letkf(X, d, XA, E, [[2*j-1;2*j] for j in 1:nr], [collect(nt*j-1:nt*j) for j in 1:nr], diagR=true)
    A0 = mean(A[1:2:end, :], dims=2)
    n = mean(A[2:2:end, :], dims=2)
    A0, n
end

"""
Assimilate SWOT observations for river reach.

- `H`: water surface elevation
- `W`: water surface width
- `x`: downstream distance for each cross section
- `wbf`: bankfull width
- `hbf`: bankfull depth
- `Qₚ`: prior probability distribution for discharge
- `nₚ`: prior probability distribution for roughness coefficient
- `rₚ`: prior probability distribution for channel shape parameter
- `zₚ`: prior distribution for downstream bed elevation
- `nens`: ensemble size
- `ri`: reach definition indices
- `ϵₒ`: observation error standard deviation
- `cv_thresh`: threshold for using dynamic or time-constant bed slope

"""
function assimilate(H, W, x, wbf, hbf, Qₚ, nₚ, rₚ, zₚ, nens, ri; ϵₒ=0.01, logQ=false)
    Qa = zeros(length(ri)-1, size(H, 2))
    S = diff(H, dims=1) ./ diff(x)
    S = [S[1, :]'; S]
    So = bed_slope(S, H, W, x, hbf, wbf, Qₚ, nₚ, rₚ, zₚ, nens)
    ze = zeros(length(x), nens)
    Qe, ne, re, ze[1, :] = lhs_ensemble(nens, Qₚ, nₚ, rₚ, zₚ)
    for t in 1:size(H, 2)
        rri = [ri[1:end-1]; length(x) + 1]
        he = gvf_ensemble!(H[:, t], W[:, t], So, x, hbf, wbf, Qe, ne, re, ze)
        i = findall(he[1, :] .> 0)
        h = ze .+ he .* ((re .+ 1) ./ re)'
        X = repeat(Qe[i]', outer=length(ri)-1)
        XA = h[:, i]
        d = H[:, t]
        E = rand(Normal(ϵₒ, 1e-6), length(d), length(i)) .* rand([-1, 1], length(d), length(i))
        if logQ
            A = letkf(log.(X), d, XA, E, [[j] for j in 1:length(ri)-1],
                      [collect(rri[j]:rri[j+1]-1) for j in 1:length(ri)-1], diagR=true)
            Qa[:, t] = mean(exp.(A), dims=2)
        else
            A = letkf(X, d, XA, E, [[j] for j in 1:length(ri)-1],
                      [collect(rri[j]:rri[j+1]-1) for j in 1:length(ri)-1], diagR=true)
            A[A .< 0] .= 0.0
            Qa[:, t] = mean(A, dims=2)  # could also use the absolute of the mean
        end
    end
    # remove outliers
    outlier = findall(mean(Qa,dims=2)[:, 1] ./ mean(Qa) .> 1.5)
    Qa[outlier, :] .= mean(Qa[[i for i in setdiff(Set(1:size(Qa, 1)), Set(outlier))], :], dims=1)
    # estimate SWOT parameters
    A0, n = estimate_Q_params(H, W, S, ri, ze, re, ne, Qa)
    Qa, A0, n
end

"""
Assimilate SWOT observations and apply hydraulic geometry constraints to estimate discharge.

"""
function assimilateHG(H, W, x, wbf, hbf, Qₚ, nₚ, rₚ, zₚ, nens, ri, ϵₒ=0.01, loc=0)
    Qa = zeros(length(ri)-1, size(H, 2))
    S = diff(H, dims=1) ./ diff(x)
    S = [S[1, :]'; S]
    So = bed_slope(S, H, W, x, hbf, wbf, Qₚ, nₚ, rₚ, zₚ, nens)
    ze = zeros(length(x), nens)
    Qe, ne, re, ze[1, :] = lhs_ensemble(nens, Qₚ, nₚ, rₚ, zₚ)
    for t in 1:size(H, 2)
        rri = [ri[1:end-1]; length(x) + 1]
        he = gvf_ensemble!(H[:, t], W[:, t], So, x, hbf, wbf, Qe, ne, re, ze)
        h = he .* ((re .+ 1) ./ re)' .+ ze
        i = findall([!((any(diff(h[:, e]) .< 0)) | (he[1, e] < 0))  for e in 1:size(h, 2)]) # i = findall(he[1, :] .> 0)
        X = repeat(Qe[i]', outer=length(ri)-1)
        if loc == 0
            XA = zeros(length(x), length(i))
            for j in 1:length(x)
                py = polyfit(log.(Qe[i]), log.(he[j, i]), 1)
                c = exp(py[0])
                f = py[1]
                XA[j, :] = c .* Qe[i]'.^f .* ((re[i] .+ 1) ./ re[i])'
            end
            XA = XA .+ ze[:, i]
            d = H[:, t]
        else
            py = polyfit(log.(Qe[i]), log.(he[loc, i]), 1)
            c=exp(py[0])
            f=py[1]
            XA = c .* X.^f .* ((re[i] .+ 1) ./ re[i])' .+ ze[loc, i]
            d = [H[loc, t]]
        end
        E = rand(Normal(0.1, 1e-6), length(d), length(i)) .* rand([-1, 1], length(d), length(i))
        A = letkf(X, d, XA, E, [[j] for j in 1:length(ri)-1], [collect(rri[j]:rri[j+1]-1) for j in 1:length(ri)-1], diagR=true)
        A[A .<= 0] .= 0.0
        A = reshape(mean(A, dims=1), 1, length(i))
        Qc = zeros(nens)
        Qc .= mean(A)
        Qc[i] = A
        he = gvf_ensemble!(H[:, t], W[:, t], So, x, hbf, wbf, Qc, ne, re, ze)
        i = findall(he[1, :] .> 0)
        h = ze .+ he .* ((re .+ 1) ./ re)'
        X = repeat(Qc[i]', outer=length(ri)-1)
        XA = h[:, i]
        d = H[:, t]
        E = rand(Normal(ϵₒ, 1e-6), length(d), length(i)) .* rand([-1, 1], length(d), length(i))
        A = Sads.letkf(X, d, XA, E, [[j] for j in 1:length(ri)-1],
                       [collect(rri[j]:rri[j+1]-1) for j in 1:length(ri)-1], diagR=true)
        A[A .< 0] .= 0.0;
        Qa[:, t] = mean(A, dims=2)  # could also use the absolute of the mean
    end
    # remove outliers
    outlier = findall(mean(Qa, dims=2)[:, 1] ./ mean(Qa) .> 1.5)
    #Qa[outlier, :] .= mean(Qa[[i for i in setdiff(Set(1:size(Qa, 1)), Set(outlier))], :], dims=1)
    Qa
end

"""
Estimate the prior probability distributions of bed elevation and discharge using rejection sampling.

"""
function priors(qwbm, H, W, x, nₚ, rₚ, nsamples, nens)
    zbnds, qbnds = prior_bounds(qwbm, nsamples, H, W, x, nₚ, rₚ)
    dm, ds = discharge_prior(qwbm, nens, nsamples, H, W, x, nₚ, rₚ, Normal(minimum(H[1, :]), 0.1))
    zm, zs = bed_elevation_prior(qwbm, nens, nsamples, H, W, x, nₚ, rₚ, zbnds, [dm, ds], qbnds)
    zm, zs, dm, ds, zbnds, qbnds
end

"""
Estimate the prior probability distribution of discharge.

"""
function discharge_prior(qwbm, nens, nsamples, H, W, x, nₚ, rₚ, zₚ)
    S = diff(H, dims=1) ./ diff(x)
    S = [S[1, :]'; S]
    S0 = [mean(S[j, :]) > 0 ? mean(S[j, :]) : maximum(S[j, :]) for j in 1:size(S, 1)]
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
        he = gvf_ensemble!(H[:, t[s]], W[:, t[s]], S0, x, maximum(H, dims=2), maximum(W, dims=2), Qe, [mean(ne) for i in 1:nens], [mean(re) for i in 1:nens], ze)
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
function bed_elevation_prior(qwbm, nens, nsamples, H, W, x, nₚ, rₚ, zbnds, dpars, qbnds)
    S = diff(H, dims=1) ./ diff(x)
    S = [S[1, :]'; S]
    S0 = [mean(S[j, :]) > 0 ? mean(S[j, :]) : maximum(S[j, :]) for j in 1:size(S, 1)]
    # S0 = [mean(S[j, :]) > 0 ? mean(S[j, :]) : minimum(S[j, :][S[j, :] .> 0]) for j in 1:size(S, 1)]
    dm, ds = dpars
    Qₚ = Truncated(LogNormal(log(qwbm/sqrt(1+dm^2)), log(1+dm^2)), qbnds[1], qbnds[2])
    zₚ = Uniform(zbnds...)
    Qe, ne, re, ze = prior_ensemble(x, Qₚ, nₚ, rₚ, zₚ, nsamples)
    So = bed_slope(S, H, W, x, maximum(H, dims=2), maximum(W, dims=2), Qₚ, nₚ, rₚ, zₚ, nens)
    he = gvf_ensemble!(mean(H, dims=2), mean(W, dims=2), So, x, maximum(H, dims=2), maximum(W, dims=2), Qe, ne, re, ze)
    h  = he .* ((re .+ 1) ./ re)' .+ ze
    i = findall([!((any(diff(h[:, e]) .< 0)) | (he[1, e] < 0))  for e in 1:size(h, 2)])
    X = reshape(ze[1, i], 1, length(i))
    XA  = he[:,i] .* ((re[i] .+ 1) ./ re[i])' .+ ze[:,i]
    d = mean(H, dims=2)[:, 1]
    E = rand(Normal(0.1, 1e-6), length(d), length(i)) .* rand([-1, 1], length(d), length(i))
    A = letkf(X, d, XA, E, [[1]], [collect(1:length(d))], diagR=true)
    zm = mean(A)
    zs = std(A)
    zm, zs
end

"""
Estimate the bounds of the prior distributions.

"""
function prior_bounds(qwbm, nsamples, H, W, x, nₚ, rₚ)
    S = diff(H, dims=1) ./ diff(x)
    S = [S[1, :]'; S]
    S0 = [mean(S[j, :]) > 0 ? mean(S[j, :]) : maximum(S[j, :]) for j in 1:size(S, 1)]
    # We will use arbitrary values that are large enough to represent the uninformative priors
    zₚ = Uniform(minimum(H[1, :]) - 20, minimum(H[1, :]))
    Qₚ = Uniform(qwbm / 10, qwbm * 10)
    ze = zeros(length(x), nsamples)
    Qe, ze[1, :] = lhs_ensemble(nsamples, Qₚ, zₚ)
    re = rand(rₚ, nsamples)
    ne = rand(nₚ, nsamples)
    ze = zeros(length(x), nsamples)
    Qe, ze[1, :] = lhs_ensemble(nsamples, Qₚ, zₚ)
    he = gvf_ensemble!(mean(H, dims=2), mean(W, dims=2), S0, x, maximum(H, dims=2), maximum(W, dims=2), Qe, [mean(ne) for i in 1:nsamples], [mean(re) for i in 1:nsamples], ze)
    i = findall(he[1, :] .> 0)
    h = he[end, i] .* (mean(re) .+ 1) ./ mean(re) .+ ze[end, i]
    j = i[(h .> minimum(H[end, :])) .& (h .< maximum(H[end, :]))]
    zbnds = [minimum(ze[1, j]), maximum(ze[1, j])]
    qbnds = [minimum(Qe[j]), maximum(Qe[j])]
    zbnds, qbnds
end

"""
Estimate At-a-station Hydraulic Geometry coefficients.

"""
function estimate_ahg(r, n, wbf, ybf, S)
    p = 2 / 3
    q = 1 / 2
    κ = S^q / n
    δ = 1 + r + r * p
    R = (1 + r) / r
    α = wbf.^((r + r * p) / δ) .* (ybf ./ R).^(-(1 + p) / δ) .* κ.^(-1 / δ)
    c = wbf.^(-r / δ) .* (ybf / R).^(1 / δ) .* κ.^(-r / δ)
    k = wbf.^(-r * p / δ) .* (ybf / R).^(p / δ) .* κ.^((1 + r) / δ)
    b = 1 / δ
    f = r / δ
    m = p * r / δ
    return α, c, k, b, f, m
end

end
