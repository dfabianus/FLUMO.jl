using FLUMO-modeling-paper
using Documenter

DocMeta.setdocmeta!(FLUMO-modeling-paper, :DocTestSetup, :(using FLUMO-modeling-paper); recursive=true)

makedocs(;
    modules=[FLUMO-modeling-paper],
    authors="Fabian Müller",
    repo="https://github.com/dfabianus/FLUMO-modeling-paper.jl/blob/{commit}{path}#{line}",
    sitename="FLUMO-modeling-paper.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://dfabianus.github.io/FLUMO-modeling-paper.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/dfabianus/FLUMO-modeling-paper.jl",
    devbranch="master",
)
