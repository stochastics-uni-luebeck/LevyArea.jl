using Documenter, LevyArea

# Plots (GR) workaround
# see https://discourse.julialang.org/t/generation-of-documentation-fails-qt-qpa-xcb-could-not-connect-to-display/60988
ENV["GKSwstype"] = "100"

# Make sure that doctests in docstrings have the package available
DocMeta.setdocmeta!(LevyArea, :DocTestSetup, :(using LevyArea); recursive=true)

makedocs(
    sitename = "LevyArea.jl",
    modules = [LevyArea],
    authors = "Felix Kastner <kastner.felix@gmail.com>",
    repo = "https://github.com/stochastics-uni-luebeck/LevyArea.jl/blob/{commit}{path}#L{line}",
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Function reference" => "functions.md",
        "Algorithms" => "algorithms.md",
        "References" => "references.md"
    ],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://stochastics-uni-luebeck.github.io/LevyArea.jl",
        assets = String[],
        ansicolor=true,
    ),
    strict = true
)

deploydocs(;
    repo = "github.com/stochastics-uni-luebeck/LevyArea.jl",
    devbranch = "main",
    push_preview = true
)
