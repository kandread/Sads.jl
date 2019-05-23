using Pkg

Pkg.update()
Pkg.add(Pkg.PackageSpec(url="https://github.com/kandread/Sads.jl"))
Pkg.add("PackageCompiler")
Pkg.add("Distributions")
Pkg.add("NCDatasets")


