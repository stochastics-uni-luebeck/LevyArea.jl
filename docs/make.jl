using Documenter, IteratedIntegrals

# Make sure that doctests in docstrings have the package available
DocMeta.setdocmeta!(IteratedIntegrals, :DocTestSetup, :(using IteratedIntegrals); recursive=true)

makedocs(
    sitename = "IteratedIntegrals.jl",
    modules = [IteratedIntegrals],
    authors = "Felix Kastner <kastner.felix@gmail.com>",
    repo = "https://github.com/fkastner/IteratedIntegrals.jl/blob/{commit}{path}#L{line}",
    pages = [
        "Home" => "index.md",
        "Index" => "functions.md",
        "References" => "references.md"
    ],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://fkastner.github.io/IteratedIntegrals.jl",
        assets = String[],
    ),
    strict = true
)

deploydocs(;
    repo = "github.com/fkastner/IteratedIntegrals.jl",
    push_preview = true
)