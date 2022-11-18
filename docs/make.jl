using SpectraFromScratch
using Documenter

DocMeta.setdocmeta!(SpectraFromScratch, :DocTestSetup, :(using SpectraFS); recursive=true)

makedocs(;
    modules=[SpectraFromScratch],
    authors="G Jake Gebbie",
    repo="https://github.com/ggebbie/SpectraFS.jl/blob/{commit}{path}#{line}",
    sitename="SpectraFromScratch.jl",
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
    repo="github.com/ggebbie/SpectraFromScratch.jl",
    devbranch="main",
)
