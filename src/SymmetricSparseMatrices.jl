module SymmetricSparseMatrices

import Base: size, getindex, show, *

export SymmetricSparseMatrix
export nnz, getindex, show, size, *

"""
    SymmetricSparseMatrix{T<:Real, I<:Integer}
Symmetric Sparse matrix implementation.
Extension from two-dimensional arrays (or array-like types) with
elements of type `T` and `I`. Alias for [`AbstractArray{T,2}`](@ref).
"""
    struct SymmetricSparseMatrix{T<:Real, I<:Integer} <: AbstractMatrix{T}
        m      ::Int                # Number of rows
        n      ::Int                # Number of columns
        rowptr ::Vector{I}          # Row i is in rowptr[i]:(rowptr[i+1]-1) 
        colval ::Vector{I}          # Col indices of stored values
        nzval  ::Vector{T}          # Stored values, typically nonzeros

        function SymmetricSparseMatrix{T,I}(m::Integer, n::Integer, rowptr::Vector{I}, colval::Vector{I}, nzval::Vector{T}) where {T<:Real,I<:Integer}
            rowptr[end]-1 == length(nzval) == length(colval) || throw(ArgumentError(string("Length of colval and nzval arrays does not match with rowptr[end]-1 ", "(",rowptr[end]-1,")")))
            new(Int(m), Int(n), rowptr, colval, nzval)
        end
    end # SymmetricSparseMatrix


"""
    SymmetricSparseMatrix(m::Integer, n::Integer, rowptr::Vector, colval::Vector, nzval::Vector)
SymetricSparseMatrix outter constructor method with full signature.
"""
function SymmetricSparseMatrix(m::Integer, n::Integer, rowptr::Vector, colval::Vector, nzval::Vector)
    Tv = eltype(nzval)
    Ti = promote_type(eltype(rowptr), eltype(colval))
    SymmetricSparseMatrix{Tv,Ti}(m, n, rowptr, colval, nzval)
end

"""
    SymmetricSparseMatrix()
SymetricSparseMatrix outter constructor method to create empty .
"""
SymmetricSparseMatrix() = SymmetricSparseMatrices.SymmetricSparseMatrix(0,0,[1],Vector{Int}(),Vector{Int}())

"""
    show(io::IO, A::SymmetricSparseMatrix)
SymetriSparseMatrix show method implementation.
"""
    show(io::IO, A::SymmetricSparseMatrix) = dump(A::SymmetricSparseMatrix)
#show(io::IO, m::MIME"text/plain", A::SymmetricSparseMatrix) = show(io::IO, A::SymmetricSparseMatrix) 

"""
    size(A::SymmetricSparseMatrix)
Returns the size (#rows, #columns) of the SymmetricSparseMatrix.
"""
    size(A::SymmetricSparseMatrix) = (A.m, A.n)

"""
    nnz(A::SymmetricSparseMatrix)
Returns the number of stored (filled) elements in a SymmetricSparseMatrix.
"""
    nnz(A::SymmetricSparseMatrix) = A.rowptr[end]-1

"""
    getindex(A::SymmetricSparseMatrix{T}, x, y) where {T<:Real}
Return the (x, y) value of the SymmetricSparseMatrix.
"""
    function getindex(A::SymmetricSparseMatrix{T,I}, x::Integer, y::Integer) where {T, I}
        if !(x in 1:A.m) || !(y in 1:A.n) throw(BoundsError(string(size(A), " SymmetricSparseMatrix at index ", "[",x,",", y,"]"))) end
        x_= min(x,y)
        y_= max(x,y)
        lbound=A.rowptr[x_]
        ubound=A.rowptr[x_+1]-1
        if ubound - lbound  < 0
            return zero(T)
        else 
            index = findfirst(isequal(y_), A.colval[lbound:ubound])
            return index == nothing ? zero(T) : A.nzval[index+lbound-1]
        end
    end # getindex

"""
    *(A::SymmetricSparseMatrix, v::Vector{T}) where {T<:Real}
Return the A*v product where `A` is SymmetricSparseMatrix{T,I} and `v` a Vector{T}.
"""
    function *(A::SymmetricSparseMatrix, v::Vector{T}) where {T<:Real}
        b = size(v) == (A.m ,) ? zero(v) : throw(ArgumentError(string("v::Vector{T} size must match with SymmetricSparseMatrix number of rows", "(",A.m,")")))
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
    end # *

end # module
