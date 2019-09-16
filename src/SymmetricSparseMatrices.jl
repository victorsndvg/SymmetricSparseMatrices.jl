module SymmetricSparseMatrices

import Base: size, getindex, show, *

"""
    SymmetricSparseMatrix{T}
Symmetric Sparse matrix implementation.
Extension from two-dimensional arrays (or array-like types) with
elements of type `T`. Alias for [`AbstractArray{T,2}`](@ref).
"""
struct SymmetricSparseMatrix{T<:Real} <: AbstractMatrix{T}
    m      ::Int                # Number of rows
    n      ::Int                # Number of columns
    nnz    ::Int                # Number of nonzero entries
    rowptr ::Vector{<:Integer}  # Row i is in rowptr[i]:(rowptr[i+1]-1) 
    colval ::Vector{<:Integer}  # Col indices of stored values
    nzval  ::Vector{T}          # Stored values, typically nonzeros
end # struct SymmetricSparseMatrix

"""
    show(io::IO, A::SymmetricSparseMatrix)
SymetriSparseMatrix show method implementation.
"""
show(io::IO, A::SymmetricSparseMatrix) = dump(A::SymmetricSparseMatrix)
#show(io::IO, m::MIME"text/plain", A::SymmetricSparseMatrix) = show(io::IO, A::SymmetricSparseMatrix) 

"""
    size(A::SymmetricSparseMatrix) = (A.m, A.n)
SymetriSparseMatrix size method implementation.
"""
size(A::SymmetricSparseMatrix) = (A.m, A.n)

"""
    getindex(A::SymmetricSparseMatrix, x, y) 
SymetriSparseMatrix getindex method implementation.
"""
function getindex(A::SymmetricSparseMatrix, x, y) 
    x_=min(x,y)
    y_=max(x,y)
    lbound=A.rowptr[x_]
    ubound=A.rowptr[x_+1]-1
    if ubound - lbound  == 0
        return 0
    else 
        index = findfirst(isequal(y_), A.colval[lbound:ubound])
        return index == nothing ? 0 : A.nzval[index+lbound-1]
    end

end

"""
    getindex(A::SymmetricSparseMatrix, x, y) 
SymetriSparseMatrix product method implementation.
"""
function *(A::SymmetricSparseMatrix, v::Vector{T}) where {T<:Real}
    b = zeros(T, size(v))
    for x in 1:A.m
        lbound=A.rowptr[x]
        ubound=A.rowptr[x+1]-1
        for i in lbound:ubound
            y = A.colval[i]
            if x==y
                b[x] += A.nzval[i] * v[x]
            elseif y>x
                b[x] += A.nzval[i] * v[x]
                b[y] += A.nzval[i] * v[y] 
            end
        end
    end
    
    return b
end


end # module
