"""
    LMPS(left::AbstractTensor, tensors::Vector{AbstractTensor})

LMPS represents a Matrix Product State (MPS) with a left boundary tensor of order 2 (`left`) and a vector of tensors of order 3 (`tensors`).

# Arguments
- `left::AbstractTensor`: The left boundary tensor of order 2.
- `tensors::Vector{AbstractTensor}`: The vector of tensors of order 3.

"""
mutable struct LMPS{T}
    left::Array{T, 2}
    tensors::Vector{Array{T, 3}}
end

"""
    left_canonicalize!(mps::LMPS, atol::Real=0)

Perform left canonicalization on an LMPS.

# Arguments
- `mps::LMPS`: The LMPS to be left canonicalized.
- `atol::Real=0`: Truncation tolerance for performing SVD.

# Output
An LMPS which consists of an isometric left boundary tensor of order 2 and (n-1) left canonical tensors of order 3.
"""
function left_canonicalize!(mps::LMPS, atol::Real=0)
    # Perform left canonicalization
    l, r = mps.left, mps.tensors[1] 
    ml, mr = reshape(l, :, size(l)[end]), reshape(r, size(r, 1), :)
    u, s, v, trunc = truncated_svd(ml * mr, size(ml, 2), atol)
    mps.left = u

    mps.tensors[1] = reshape(Diagonal(s) * v, size(v, 1), size(r, 2), size(r, 3))
    for i = 1:length(mps.tensors)-1
        left, right = mps.tensors[i], mps.tensors[i+1]
        mleft, mright = reshape(left, :, size(left)[end]), reshape(right, size(right, 1), :)

        u, s, v, trunc = truncated_svd(mleft * mright, size(mleft, 2), atol)
        mps.tensors[i] = reshape(u, size(left, 1), size(left, 2), size(u, 2))
        mps.tensors[i+1] = reshape(Diagonal(s) * v, size(v, 1), size(right, 2), size(right, 3))
    end

    pop!(mps.tensors)
    
    return mps
end