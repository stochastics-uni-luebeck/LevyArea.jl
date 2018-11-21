using Documenter, SRK

makedocs(
    sitename = "SRK.jl",
    modules = [SRK],
    pages = [
        "Home" => "index.md"
    ]
)

deploydocs(
    repo = "github.com/mrbauff/SRK.jl.git"
)
