function truncated_svd(tmat::AbstractMatrix, dmax::Int, atol::Real)
    u, s, v = LinearAlgebra.svd(tmat)
    dmax = min(searchsortedfirst(s, atol, rev=true), dmax, length(s))
    return u[:, 1:dmax], s[1:dmax], v'[1:dmax, :], sum(s[dmax+1:end])
end