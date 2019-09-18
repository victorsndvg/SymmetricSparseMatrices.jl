using Documenter
using SymmetricSparseMatrices

makedocs(
    sitename = "SymmetricSparseMatrices",
    format = Documenter.HTML(),
    modules = [SymmetricSparseMatrices]
)

deploydocs(repo = "github.com/victorsndvg/SymmetricSparseMatrices.jl.git")
