# module Swot

using Sads
using Distributions
using NCDatasets

function load_severn(ncfile)
    res = 90.0
    ds = Dataset(ncfile)
    xs = ds.group["XS_Timseries"]
    W = xs["W"][:]; W = W[end:-1:1, :]
    H = xs["H_1km"][:]; H = H[end:-1:1, :]
    x = collect(0.0:res:(size(H, 1)-1)*res)
    ri = [1;length(x)]
    Q = xs["Q"][:]; Q = Q[end:-1:1, :]
    qwbm = mean(Q[Q .> 0])
    close(ds)
    Qr = zeros(size(Q))
    Wr = zeros(size(W))
    Hr = zeros(size(H))
    for t in 1:size(H,2)
        j = sortperm(H[:,t])
        Qr[:,t] .= Q[j,t]
        Wr[:,t] .= W[j,t]
        Hr[:,t] .= H[j,t]
    end
    Qr = Qr[:, findall(all(.!isnan.(H), dims=1)[1,: ])]
    Wr = Wr[:, findall(all(.!isnan.(H), dims=1)[1,: ])]
    Hr = Hr[:, findall(all(.!isnan.(H), dims=1)[1,: ])]
    Wr[Wr .== 0] .= mean(Wr) 
    Hr[Hr .== 0] .= mean(Hr)
    Qr, Hr, Wr, x, qwbm, ri
end

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
function write_data(ncfile, H, W, x, ri, Qa, A0, n)
    nr, nt = size(Qa)
    ds = Dataset(ncfile, "c")
    defDim(ds, "reach", nr)
    defDim(ds, "time", nt)
    Hr = [mean(H[ri[j]:ri[j+1], t]) for j in 1:nr, t in 1:nt]
    Wr = [mean(W[ri[j]:ri[j+1], t]) for j in 1:nr, t in 1:nt]
    Sr = [(H[ri[j+1], t] - H[ri[j], t]) / (x[ri[j+1]] - x[ri[j]]) for j in 1:nr, t in 1:nt]
    Hvar = defVar(ds, "H", Float32, ("reach", "time"))
    Hvar.attrib["units"] = "m"
    Hvar.attrib["comments"] = "Water Surface Elevation"
    Wvar = defVar(ds, "W", Float32, ("reach", "time"))
    Wvar.attrib["units"] = "m"
    Wvar.attrib["comments"] = "Width"
    Svar = defVar(ds, "S", Float32, ("reach", "time"))
    Svar.attrib["units"] = "m/m"
    Svar.attrib["comments"] = "Water Surface Slope"
    Qvar = defVar(ds, "Q", Float32, ("reach", "time"))
    Qvar.attrib["units"] = "m3/s"
    Qvar.attrib["comments"] = "Discharge"
    nvar = defVar(ds, "n", Float32, ("reach", "time"))
    nvar.attrib["units"] = "unitless"
    nvar.attrib["comments"] = "Manning's roughness"
    A0var = defVar(ds, "A0", Float32, ("reach", "time"))
    Hvar[:, :] = Hr
    Wvar[:, :] = Wr
    Svar[:, :] = Sr
    Qvar[:, :] = Qa
    A0var[:, :] = A0
    nvar[:, :] = n
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
    Q, H, W, x, qwbm, ri = load_severn(ncfile)
    Qa, A0, n = Sads.estimate(qwbm, H, W, x, rₚ, ri)
    write_data("sads_$ncfile", H, W, x, ri, Qa, A0, n)
end

main()

# Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
#     main(ARGS)
#     return 0
# end

# end
