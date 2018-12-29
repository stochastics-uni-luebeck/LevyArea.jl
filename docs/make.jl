using Documenter, SRK

makedocs(
    sitename = "SRK.jl",
    modules = [SRK],
    pages = [
        "Home" => "index.md"
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true")
)

deploydocs(
    repo = "github.com/mrbauff/SRK.jl.git"
)
