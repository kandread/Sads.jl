# module Swot

using Sads
using Distributions
using NCDatasets


"Load hydraulic variables from Pepsi-1 NetCDF file."
function load_data(ncfile)
    ds = Dataset(ncfile)
    xs = NCDatasets.group(ds, "XS_Timeseries")
    ri = NCDatasets.group(ds, "River_Info")
    rbnd = ri["rch_bnd"][:]
    qwbm = ri["QWBM"][1]
    x = xs["X"][:, 1]
    ri = [findmin(sqrt.((r .- x).^2))[2] for r in rbnd]
    x = x[end] .- x; x = x[end:-1:1]
    ri = length(x) .- ri .+ 1; ri = ri[end:-1:1]
    W = xs["W"][:]; W = W[end:-1:1, :]
    H = xs["H"][:]; H = H[end:-1:1, :]
    Q = xs["Q"][:]; Q = Q[end:-1:1, :]
    close(ds)
    Q, H, W, x, qwbm, ri
end

"Read data from SWOT NetCDF file."
function read_data(ncfile)
    ds = Dataset(ncfile)
    x = ds["XS_90m"][:] .* 90.0 .- 90
    xs = NCDatasets.group(ds, "XS_Timseries")
    W = xs["W"][:]
    H = xs["H_1km"][:]
    Q = xs["Q"][:]
    S = abs.(NCDatasets.group(ds, "Reach_Timeseries")["S_90m"][:]) ./ length(x)
    close(ds)
    j = findall((sum(isnan.(H), dims=1)[1, :] .== 0) .& (all(W .> 0, dims=1)[1, :]))
    # H = sort(H[:, j], dims=1)
    W = W[:, j]; W = W[end:-1:1, :]
    Hmin = minimum(H[:, j], dims=1)
    H = zeros(size(W))
    H[1, :] = Hmin
    for c in 2:size(H, 1)
        for t in 1:size(H, 2)
            H[c, t] = H[c-1, t] + S[t] .* (x[c] - x[c-1])
        end
    end
    Q = Q[:, j]; Q = Q[end:-1:1, :]
    qwbm = mean(abs.(Q))
    ri = [1; length(x)]
    Q, H, W, x, qwbm, ri
end

"Write output data to NetCDF file."
function write_data(ncfile, H, W, S, Qa, A0, n)
    nr, nt = size(Qa)
    ds = Dataset(ncfile, "c")
    defDim(ds, "nx", nr)
    defDim(ds, "nt", nt)
    Hvar = defVar(ds, "H", Float32, ("nx", "nt"))
    Hvar.attrib["units"] = "m"
    Hvar.attrib["comments"] = "Water Surface Elevation"
    Wvar = defVar(ds, "W", Float32, ("nx", "nt"))
    Wvar.attrib["units"] = "m"
    Wvar.attrib["comments"] = "Width"
    Svar = defVar(ds, "S", Float32, ("nx", "nt"))
    Svar.attrib["units"] = "m/m"
    Svar.attrib["comments"] = "Water Surface Slope"
    Qvar = defVar(ds, "Q", Float32, ("nx", "nt"))
    Qvar.attrib["units"] = "m3/s"
    Qvar.attrib["comments"] = "Discharge"
    nvar = defVar(ds, "n", Float32, ("nx",))
    nvar.attrib["units"] = "unitless"
    nvar.attrib["comments"] = "Manning's roughness"
    A0var = defVar(ds, "Abase", Float32, ("nx",))
    Hvar[:, :] = H
    Wvar[:, :] = W
    Svar[:, :] = S
    Qvar[:, :] = Qa
    A0var[:] = A0
    nvar[:] = n
    close(ds)
end

"Calculate Nash-Sutcliffe efficiency."
function nse(o, m)
    1.0 .- sum((m .- o).^2) ./ sum((o .- mean(o)).^2)
end

"Calculate relative root mean squared error."
function rrmse(o, m)
    sqrt(mean(((o .- m) ./ o).^2))
end

"Calculate normalized root mean squared error."
function nrmse(o, m)
    sqrt(mean(((o .- m) ./ mean(o)).^2))
end

"Compare theoretical with estimated AHG coefficients."
function compare_ahg(Qa, H, W, x, qwbm, rₚ, nₚ=Uniform(0.02, 0.05), nens=100, nsamples=1000)
    est = zeros(length(x), 6)
    the = zeros(length(x), 6)
    wbf, hbf, hmin = maximum(W, dims=2), maximum(H, dims=2), minimum(H[1, :])
    S = diff(H, dims=1) ./ diff(x)
    S = [S[1, :]'; S]
    zm, zs, dm, ds, zb, qb = Sads.priors(qwbm, H, W, x, nₚ, rₚ, nsamples, nens)
    Qₚ = Truncated(LogNormal(log(qwbm/sqrt(1+dm^2)), log(1+dm^2)), qb[1], qb[2])
    zₚ = Normal(zm, 0.1)
    So = Sads.bed_slope(S, H, W, x, hbf, wbf, Qₚ, nₚ, rₚ, zₚ, nens)
    ze = zeros(length(x), size(H,2))
    ze[1,:] .= zm
    for j in 2:length(x) ze[j, :] = ze[j-1, :] .+ So[j] * (x[j] - x[j-1]); end
    ybf = hbf .- ze[:, 1]
    n = mean(rand(nₚ, nsamples))
    r = mean(rand(rₚ, nsamples))
    he = Sads.gvf_ensemble!(mean(H,dims=2), mean(W,dims=2), So, x, maximum(H, dims=2), maximum(W, dims=2), Qa, [n for _ in 1:size(H,2)], [r for _ in 1:size(H,2)], ze)
    i = findall(he[1, :] .> 0)
    h = he .* (r+1)/r .+ ze
    S = diff(h, dims=1) ./ diff(x)
    S = [S[1, :]'; S]
    for xi=1:length(x)
        py = polyfit(log.(Qa[i]), log.(he[xi,i]), 1)
        w = wbf[end] * (he[xi, i] ./ ybf[xi]).^(1/r)
        pw = polyfit(log.(Qa[i]), log.(w), 1)
        sf = [S[xi,t] > 0 ? S[xi,t] : minimum(S[xi,S[xi,:] .> 0]) for t in 1:size(S,2)]
        U = (1/n) .* he[end,i].^(2/3) .* sf[i].^(1/2)
        pu = polyfit(log.(Qa[i]), log.(U), 1)
        a, c, k, b, f, m = Sads.estimate_ahg(r, n, wbf[xi], maximum(he[xi,i]), mean(sf))
        est[xi, :] = [exp(pw[0]), exp(py[0]), exp(pu[0]), pw[1], py[1], pu[1]]
        the[xi, :] = [a, c, k, b, f, m]
    end
    est, the
end

"Main driver routine."
function main()
    ncpath = ARGS[1]
    rₚ = eval(Symbol(ARGS[2]))(parse(Float64, ARGS[3]), parse(Float64, ARGS[4]))
    ncfiles = joinpath.(abspath(ncpath), readdir(ncpath))
    nt = maximum([length(Dataset(ncfile)["Time steps"]) for ncfile in ncfiles])
    Qa = zeros(length(ncfiles), nt)
    A0 = zeros(length(ncfiles))
    n = zeros(length(ncfiles))
    H = zeros(length(ncfiles), nt)
    W = zeros(length(ncfiles), nt)
    S = zeros(length(ncfiles), nt)
    for (i, ncfile) in enumerate(ncfiles)
        Q, h, w, x, qwbm, ri = read_data(ncfile)
        qa, a0, mann = Sads.estimate(qwbm, h, w, x, rₚ, ri)
        H[i, 1:size(h, 2)] = mean(h, dims=1)
        W[i, 1:size(w, 2)] = mean(w, dims=1)
        S[i, 1:size(w, 2)] = mean(diff(h, dims=1) ./ diff(x), dims=1)
        Qa[i, 1:size(h, 2)], A0[i], n[i] = mean(qa, dims=1), mean(a0), mean(mann)
    end
    write_data("$(ncpath)/sads_output.nc", H, W, S, Qa, A0, n)
end

main()

# Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
#     main(ARGS)
#     return 0
# end

# end
