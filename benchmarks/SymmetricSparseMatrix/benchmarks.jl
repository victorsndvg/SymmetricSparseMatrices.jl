using SymmetricSparseMatrices
using SparseArrays
using BenchmarkTools
using Statistics
using Test

# Define some input data for benchmarks
ZeroVector           = Vector{Int}()
OneVector            = ones(1) 
FiveVector           = ones(5)
ZeroSymCSR           = SymmetricSparseMatrix()
OneSymCSR            = SymmetricSparseMatrix(1,[1,2],[1],[1])
NonTriangInputSymCSR = SymmetricSparseMatrix(5,[1,4,6,9,12,14], [1,2,4,1,2,3,4,5,1,3,4,2,5], [1,-1,-3,-2,5,4,6,4,-4,2,7,8,-5])
DiagCSC              = sparse([1,2,3,4],[1,2,3,4],[1,2,3,4])
DiagDense            = Matrix(DiagCSC)
UpperTriangCSC       = sparse([1,2,3,4,1,2,3,4],[1,2,3,4,4,3,2,1],[1,2,3,4,4,3,2,1])
UpperTriangDense     = Matrix(UpperTriangCSC)

# Define a parent BenchmarkGroup to contain our Suite
const Suite = BenchmarkGroup(["SymmetricSparseMatrix"])

# Add some benchmarks to the "Product" group
Suite["Product"]                   = BenchmarkGroup(["Product"])
Suite["Product"]["Zero"]           = @benchmarkable $ZeroSymCSR*$ZeroVector
Suite["Product"]["One"]            = @benchmarkable $OneSymCSR*$OneVector
Suite["Product"]["NonTriangInput"] = @benchmarkable $NonTriangInputSymCSR*$FiveVector

# Add some benchmarks to the "Conversion" group
Suite["Conversion"]                     = BenchmarkGroup(["Conversion"])
Suite["Conversion"]["DiagCSC"]          = @benchmarkable symsparse($DiagCSC)
Suite["Conversion"]["DiagDense"]        = @benchmarkable symsparse($DiagDense)
Suite["Conversion"]["UpperTriangCSC"]   = @benchmarkable symsparse($UpperTriangCSC)
Suite["Conversion"]["UpperTriangDense"] = @benchmarkable symsparse($UpperTriangDense)

# If a cache of tuned parameters already exists, use it, otherwise, tune and cache
# the benchmark parameters. Reusing cached parameters is faster and more reliable
# than re-tuning `Suite` every time the file is included.
ParamsPath = joinpath(dirname(@__FILE__), "params.json")

if isfile(ParamsPath)
    loadparams!(Suite, BenchmarkTools.load(ParamsPath)[1], :evals)
else
    tune!(Suite)
    BenchmarkTools.save(ParamsPath, params(Suite))
end

# Run benchmark Suite
Results = run(Suite)

# Load reference benchmark Suite results
ReferenceResultsPath = joinpath(dirname(@__FILE__), "results.json")

if isfile(ReferenceResultsPath)
    ReferenceResults =  BenchmarkTools.load(ReferenceResultsPath)[1]
else
    ReferenceResults = Results
    BenchmarkTools.save(ReferenceResultsPath, Results)
end

# Compare current results against reference results
Comparisson = judge(mean(Results), mean(ReferenceResults))
Regressions = regressions(Comparisson)

@test length(Regressions) == 0

