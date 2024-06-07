"""
    truncated_svd(tmat::AbstractMatrix, dmax::Int, atol::Real=0)

Perform singular value decomposition (SVD) for `tmat` and keep only the largest `dmax` singular values (not less than `atol`) and associated columns and rows.

# Arguments
- `tmat::AbstractMatrix`: The matrix to perform SVD on.
- `dmax::Int`: The maximum number of singular values to keep.
- `atol::Real=0`: The minimum threshold for singular values.

# Returns
- `U`: The left singular vectors.
- `S`: The singular values.
- `Vt`: The right singular vectors.
- `sum_discarded`: The sum of discarded singular values.
"""
function truncated_svd(tmat::AbstractMatrix, dmax::Int, atol::Real=0)
    u, s, v = LinearAlgebra.svd(tmat)
    dmax = min(searchsortedfirst(s, atol, rev=true), dmax, length(s))
    return u[:, 1:dmax], s[1:dmax], v'[1:dmax, :], sum(s[dmax+1:end])
end