using SymmetricSparseMatrices
using Test

@testset "SymmetricSparseMatrix" begin
    @test size(SymmetricSparseMatrices.SymmetricSparseMatrix()) == (0,0)
    @test nnz(SymmetricSparseMatrices.SymmetricSparseMatrix()) == 0
    @test getindex(SymmetricSparseMatrices.SymmetricSparseMatrix(1,1,[1,2],[1],[1]),1,1) == 1
    @test nnz(SymmetricSparseMatrices.SymmetricSparseMatrix(1,1,[1,2],[1],[1])) == 1
    @test SymmetricSparseMatrices.SymmetricSparseMatrix(5,5,[1,4,6,9,12,14], [1,2,4,1,2,3,4,5,1,3,4,2,5], [1,-1,-3,-2,5,4,6,4,-4,2,7,8,-5])*ones(5) == [-3, 4, 14, 10, -1]
end
