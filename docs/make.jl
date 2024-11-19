using MatrixFuns
using Documenter

DocMeta.setdocmeta!(MatrixFuns, :DocTestSetup, :(using MatrixFuns); recursive=true)

makedocs(;
    modules=[MatrixFuns],
    authors="xquan818 <xuequan818@gmail.com> and contributors",
    sitename="MatrixFuns.jl",
    format=Documenter.HTML(;
        canonical="https://xuequan818.github.io/MatrixFuns.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Background" => [
            "Computing Matrix Functions" => "background/matfun.md",
            "Divided Differences" => "background/divdiff.md",
            "Fréchet Derivatives" => "background/frechet.md",
            "Limilations" => "background/limit.md",
        ],
        "Examples" => "examples.md",
        "API" => "api.md"],
)

deploydocs(;
    repo="github.com/xuequan818/MatrixFuns.jl",
    devbranch="master",
)
