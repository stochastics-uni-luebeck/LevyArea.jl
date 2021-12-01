using Documenter, LevyArea

# Make sure that doctests in docstrings have the package available
DocMeta.setdocmeta!(LevyArea, :DocTestSetup, :(using LevyArea); recursive=true)

makedocs(
    sitename = "LevyArea.jl",
    modules = [LevyArea],
    authors = "Felix Kastner <kastner.felix@gmail.com>",
    repo = "https://github.com/stochastics-uni-luebeck/LevyArea.jl/blob/{commit}{path}#L{line}",
    pages = [
        "Home" => "index.md",
        "Index" => "functions.md",
        "References" => "references.md"
    ],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://stochastics-uni-luebeck.github.io/LevyArea.jl",
        assets = String[],
    ),
    strict = true
)

deploydocs(;
    repo = "github.com/stochastics-uni-luebeck/LevyArea.jl",
    push_preview = true
)