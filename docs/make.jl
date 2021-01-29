using Documenter
using Krippendorff

push!(LOAD_PATH,"../src/")
makedocs(
    sitename = "Krippendorff.jl",
    format = Documenter.HTML(),
    pages = [
            "Krippendorffs Alpha" => "index.md",
            "Important Functions" => "importantFunctions.md",
            "About" => "about.md"
         ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/FPGro/Krippendorff.jl.git",
    devbranch = "main"
)
