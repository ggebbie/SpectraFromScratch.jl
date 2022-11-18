using Revise, PkgTemplates

t = Template(; user="ggebbie", dir="~/projects", authors="G Jake Gebbie", julia=v"1.6", plugins=[ License(; name="MIT"), Git(; manifest=true, ssh=true), GitHubActions(; x86=false, extra_versions=["1.6","1.7","nightly"]), Codecov(), Documenter{GitHubActions}(), Develop(), ], )

t("SpectraFromScratch.jl")
