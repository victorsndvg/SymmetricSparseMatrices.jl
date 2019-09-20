using SymmetricSparseMatrices
using SparseArrays
using Test

@testset "SymmetricSparseMatrix" begin
    # SymmetricSparseMatrix
    @test size(SymmetricSparseMatrix()) == (0,0)
    @test nnz(SymmetricSparseMatrix()) == 0
    @test getindex(SymmetricSparseMatrix(1,[1,2],[1],[1]),1,1) == 1
    @test nnz(SymmetricSparseMatrix(1,[1,2],[1],[1])) == 1
    @test SymmetricSparseMatrix(5,[1,4,6,9,12,14], [1,2,4,1,2,3,4,5,1,3,4,2,5], [1,-1,-3,-2,5,4,6,4,-4,2,7,8,-5])*ones(5) == [-3, 4, 14, 10, -1]

    # Conversions
    a = sparse([1,2,3,4],[1,2,3,4],[1,2,3,4])
    b = Matrix(a)
    c = symsparse(b)
    @test a == b == c

    a = sparse([1,2,3,4,1,2,3,4],[1,2,3,4,4,3,2,1],[1,2,3,4,4,3,2,1])
    b = Matrix(a)
    c = symsparse(b)
    @test a != c != b
end
