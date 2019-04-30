#!/usr/bin/env julia

using Sads
using Distributions
using NCDatasets

"Load hydraulic variables from NetCDF file."
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

"Write output data to NetCDF file."
function write_data(ncfile, Qa)
    ds = Dataset(ncfile, "a")
    xs = NCDatasets.group(ds, "XS_Timeseries")
    ri = NCDatasets.group(ds, "River_Info")
    rc = NCDatasets.group(ds, "Reach_Timeseries")
    rbnd = ri["rch_bnd"][:]
    x = xs["X"][:, 1]
    ri = [findmin(sqrt.((r .- x).^2))[2] for r in rbnd]
    x = x[end] .- x; x = x[end:-1:1]
    ri = length(x) .- ri .+ 1; ri = ri[end:-1:1]
    Qo = zeros(length(x), size(Qa, 2))
    for r in 1:length(ri)-1
        for rr in ri[r]:ri[r+1]-1
            Qo[rr, :] .= Qa[r, :]
        end
    end
    Qxs = xs["Q"]
    Qxs[:] = Qo[end:-1:1, :]
    Qr = rc["Q"]
    Qr[:] = Qa[end:-1:1, :]
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

"Main driver routine."
function main()
    ncfile = ARGS[1]
    rₚ = eval(Symbol(ARGS[2]))(parse(Float64, ARGS[3]), parse(Float64, ARGS[4]))
    Q, H, W, x, qwbm, ri = load_data(ncfile)
    Qa = Sads.estimate(qwbm, H, W, x, rₚ, ri)
    write_data(ncfile, Qa)
end

main()
