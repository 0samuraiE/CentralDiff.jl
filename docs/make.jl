using CentralDiff
using Documenter

DocMeta.setdocmeta!(CentralDiff, :DocTestSetup, :(using CentralDiff); recursive=true)

makedocs(;
    modules=[CentralDiff],
    authors="Taketo Tominaga",
    sitename="CentralDiff.jl",
    format=Documenter.HTML(;
        canonical="https://0samuraiE.github.io/CentralDiff.jl", edit_link="master", assets=String[]
    ),
    pages=["Home" => "index.md"],
)

deploydocs(; repo="github.com/0samuraiE/CentralDiff.jl", devbranch="master")
