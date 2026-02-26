using Documenter
using PenguinSolverCore

makedocs(
    modules = [PenguinSolverCore],
    authors = "PenguinxCutCell contributors",
    sitename = "PenguinSolverCore.jl",
    format = Documenter.HTML(
        canonical = "https://PenguinxCutCell.github.io/PenguinSolverCore.jl",
        repolink = "https://github.com/PenguinxCutCell/PenguinSolverCore.jl",
        collapselevel = 2,
    ),
    pages = [
        "Home" => "index.md",
        "Callbacks and Scheduling" => "callbacks.md",
        "API Reference" => "reference.md",
    ],
    pagesonly = true,
    warnonly = true,
    remotes = nothing,
)

if get(ENV, "CI", "") == "true"
    deploydocs(
        repo = "github.com/PenguinxCutCell/PenguinSolverCore.jl",
        push_preview = true,
    )
end
