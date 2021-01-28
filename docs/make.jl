using Documenter
using Krippendorff

push!(LOAD_PATH,"../src/")
makedocs(
    sitename = "Krippendorff",
    format = Documenter.HTML(),
    pages = [
            "Index" => "index.md",
            "An other page" => "anotherPage.md",
         ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/FPGro/Krippendorff.jl.git",
    devbranch = "main"
)
