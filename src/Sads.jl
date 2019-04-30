module Sads

abstract type CrossSection end

include("crosssections.jl")
include("kalman.jl")
include("gvf.jl")
include("ensemble.jl")

"""
Estimate discharge by assimilating observed water surface elevations along a river channel.

"""
function estimate()

end

"""
Estimate the prior probability distributions of bed elevation and discharge using rejection sampling.

"""
function priors(qwbm, H, W, x, nₚ, rₚ, nsamples, nens)
    zbnds, qbnds = prior_bounds(qwbm, nsamples, H, W, x, nₚ, rₚ)
end

"""
Estimate the prior probability distribution of discharge.

"""
function discharge_prior()

end

"""
Estimate the prior probability distribution of bed elevation.

"""
function bed_elevation_prior()

end

"""
Estimate the bounds of the prior distributions.

"""
function prior_bounds(qwbm, nsamples, H, W, x, nₚ, rₚ)
    S = diff(H, dims=1) ./ diff(x)
    S = [S[1, :]'; S]
    # We will use arbitrary values that are large enough to represent the uninformative priors
    zₚ = Uniform(minimum(H[1, :]) - 30, minimum(H[1, :]))
    Qₚ = Uniform(qwbm / 10, qwbm * 10)
    ze = zeros(length(x), nsamples)
    Qe, ze[1, :] = lhs_ensemble(nsamples, Qₚ, zₚ)
    re = rand(rₚ, nsamples)
    ne = rand(nₚ, nsamples)
    he = gvf_ensemble(mean(H, dims=2), mean(W, dims=2), mean(S, dims=2), x, maximum(H, dims=2),
                      maximum(W, dims=2), Qe, [mean(ne) for i in 1:nsamples], [mean(re) for i in 1:nsamples], ze)
    i = findall(he[1, :] .> 0)
    h = he[end, i] .* (mean(re) .+ 1) ./ mean(re) .+ ze[end, i]
    j = i[(h .> minimum(H[end, :])) .& (h .< maximum(H[end, :]))]
    zbnds = [minimum(ze[1, j]), maximum(ze[1, j])]
    qbnds = [minimum(Qe[j]), maximum(Qe[j])]
end

end
