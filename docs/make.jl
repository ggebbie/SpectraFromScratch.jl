using SpectraFS
using Documenter

DocMeta.setdocmeta!(SpectraFS, :DocTestSetup, :(using SpectraFS); recursive=true)

makedocs(;
    modules=[SpectraFS],
    authors="G Jake Gebbie",
    repo="https://github.com/ggebbie/SpectraFS.jl/blob/{commit}{path}#{line}",
    sitename="SpectraFS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ggebbie.github.io/SpectraFS.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ggebbie/SpectraFS.jl",
    devbranch="main",
)
