using Documenter
using SymmetricSparseMatrices

makedocs(
    sitename = "SymmetricSparseMatrices",
    format = Documenter.HTML(),
    modules = [SymmetricSparseMatrices]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
