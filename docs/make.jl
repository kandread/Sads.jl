using Documenter, Sads

makedocs(;
    modules=[Sads],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/kandread/Sads.jl/blob/{commit}{path}#L{line}",
    sitename="Sads.jl",
    authors="Kostas Andreadis",
    assets=[],
)

deploydocs(;
    repo="github.com/kandread/Sads.jl",
)
