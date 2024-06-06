using IsometricPEPS
using Documenter

DocMeta.setdocmeta!(IsometricPEPS, :DocTestSetup, :(using IsometricPEPS); recursive=true)

makedocs(;
    modules=[IsometricPEPS],
    authors="Huanhai Zhou",
    sitename="IsometricPEPS.jl",
    format=Documenter.HTML(;
        canonical="https://fliingelephant.github.io/IsometricPEPS.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/fliingelephant/IsometricPEPS.jl",
    devbranch="main",
)
