module SymmetricSparseMatricesTests

using Test

@time @testset "SymmetricSparseMatrix" begin include("SymmetricSparseMatrix/runtests.jl") end

end # module
