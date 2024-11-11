using MatrixFuns
using Documenter

DocMeta.setdocmeta!(MatrixFuns, :DocTestSetup, :(using MatrixFuns); recursive=true)

makedocs(;
    modules=[MatrixFuns],
    authors="xquan818 <840169780@qq.com> and contributors",
    sitename="MatrixFuns.jl",
    format=Documenter.HTML(;
        canonical="https://xuequan818.github.io/MatrixFuns.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Background" => [
            "Computing Matrix Functions" => "background/matfun.md",
            "Divided difference" => "background/divdiff.md",
            "Fréchet derivative" => "background/frechet.md"
        ],
        "Examples" => "examples.md",
        "API" => "api.md"],
)

deploydocs(;
    repo="github.com/xuequan818/MatrixFuns.jl",
    devbranch="master",
)
