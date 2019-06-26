using Documenter, IteratedIntegrals

makedocs(
    sitename = "IteratedIntegrals.jl",
    modules = [IteratedIntegrals],
    pages = [
        "Home" => "index.md",
        "Index" => "functions.md",
        "References" => "references.md"
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true")
)

deploydocs(
    repo = "github.com/fkastner/IteratedIntegrals.jl.git"
)
