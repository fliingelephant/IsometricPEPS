σᶻ = [1 0; 0 -1];
σˣ = [0 1; 1 0];
I = [1 0; 0 1];

"""
    transverse_Ising_1d_seperate_effH(mps::LMPS, h::Float64, J::Float64)

This function generates the effective Hamiltonian for a 1d transverse Ising model by contracting ZZ and X terms seperately.

# Arguments
- `mps::LMPS`: The LMPS representing the state of the system.
- `J::Float64`: The coupling strength between adjacent spins.
- `h::Float64`: The strength of the transverse magnetic field.

# Returns
- `eff_H::Matrix{ComplexF64}`: The effective Hamiltonian matrix.
"""
function transverse_Ising_1d_seperate_effH(mps::LMPS, h::Float64, J::Float64)
    mpo_zz = reshape(- J * kron(σᶻ, σᶻ), 2, 2, 2, 2)
    mpo_x = - h * σˣ

    # ZZ terms
    # - J Z_{left} Z_{1}
    #=
    j                   l
    4                   3
        4-leg tensor
    2                   1
    i                   k
    =#
    eff_H = ein"ia, jb, akc, bld, kilj->cd"(mps.left, conj(mps.left), mps.tensors[1], conj(mps.tensors[1]), mpo_zz)
    for i in 2:length(mps.tensors)
        eff_H = ein"(ab, aic), bid->cd"(eff_H, mps.tensors[i], conj(mps.tensors[i]))
    end
    
    # - J Z_{i} Z_{i+1}
    for i in 1:length(mps.tensors)-1
        #=
        k                   l
        4                   3
            4-leg tensor
        2                   1
        i                   j
        =#
        tmp_eff_H = ein"(((aib, akc), jilk), cle), bjd->de"(mps.tensors[i], conj(mps.tensors[i]), mpo_zz, conj(mps.tensors[i+1]), mps.tensors[i+1])
        for j in i+2:length(mps.tensors)
            tmp_eff_H = ein"(ab, aic), bid->cd"(tmp_eff_H, mps.tensors[j], conj(mps.tensors[j]))
        end
        eff_H += tmp_eff_H
    end

    # X terms
    # - h X_{left}
    tmp_eff_H = ein"ia, jb, ij->ab"(mps.left, conj(mps.left), mpo_x)
    for i in 1:length(mps.tensors)
        tmp_eff_H = ein"(ab, aic), bid->cd"(tmp_eff_H, mps.tensors[i], conj(mps.tensors[i]))
    end
    eff_H += tmp_eff_H

    # - h X_{i}
    for i in 1:length(mps.tensors)
        tmp_eff_H = ein"(aib, ajc), ij->bc"(mps.tensors[i], conj(mps.tensors[i]), mpo_x)
        for j in i+1:length(mps.tensors)
            tmp_eff_H = ein"(ab, aic), bid->cd"(tmp_eff_H, mps.tensors[j], conj(mps.tensors[j]))
        end
        eff_H += tmp_eff_H
    end

    return eff_H
end

function transverse_Ising_1d_seperate_energy(mps::LMPS, h::Float64, J::Float64)
    return real(ein"ii->"(transverse_Ising_1d_seperate_effH(mps, J, h))[])
end

"""
    transverse_Ising_1d_effH(mps::LMPS, h::Float64, J::Float64)

This function generates the effective Hamiltonian for a 1d transverse Ising model by contracting ZZ and X terms together for every two sites.

# Arguments
- `mps::LMPS`: The LMPS representing the state of the system.
- `J::Float64`: The coupling strength between adjacent spins.
- `h::Float64`: The strength of the transverse magnetic field.

# Returns
- `eff_H::Matrix{ComplexF64}`: The effective Hamiltonian matrix.
"""
function transverse_Ising_1d_effH(mps::LMPS, h::Float64, J::Float64)
    mpo_left = reshape(- J * kron(σᶻ, σᶻ) - h * kron(σˣ, I) - h * kron(I, σˣ), 2, 2, 2, 2)
    mpo = reshape(- J * kron(σᶻ, σᶻ) - h * kron(I, σˣ), 2, 2, 2, 2)

    # - J Z_{left} Z_{1} - h X_{left} - h X_{1}
    #=
    j                   l
    4                   3
        4-leg tensor
    2                   1
    i                   k
    =#
    eff_H = ein"ia, jb, akc, bld, kilj->cd"(mps.left, conj(mps.left), mps.tensors[1], conj(mps.tensors[1]), mpo_left)
    for i in 2:length(mps.tensors)
        eff_H = ein"(ab, aic), bid->cd"(eff_H, mps.tensors[i], conj(mps.tensors[i]))
    end

    # - J Z_{i} Z_{i+1} - h X_{i+1}
    for i in 1:length(mps.tensors)-1
        #=
        k                   l
        4                   3
            4-leg tensor
        2                   1
        i                   j
        =#
        tmp_eff_H = ein"(((aib, akc), jilk), cle), bjd->de"(mps.tensors[i], conj(mps.tensors[i]), mpo, conj(mps.tensors[i+1]), mps.tensors[i+1])
        for j in i+2:length(mps.tensors)
            tmp_eff_H = ein"(ab, aic), bid->cd"(tmp_eff_H, mps.tensors[j], conj(mps.tensors[j]))
        end

        eff_H += tmp_eff_H
    end

    return eff_H
end

function transverse_Ising_1d_energy(mps::LMPS, h::Float64, J::Float64)
    return real(ein"ii->"(transverse_Ising_1d_effH(mps, J, h))[])
end

"""
    steifel_update!(mps::LMPS, G::@NamedTuple{left::Matrix{T}, tensors::Vector{Array{T, 3}}}, τ::T) where T

Update the LMPS `mps` using the gradient `G` and the updating rate `τ`. The gradient `G` should have the same shape as `mps`.

# Arguments
- `mps::LMPS`: The LMPS to update.
- `G::@NamedTuple{left::Matrix{T}, tensors::Vector{Array{T, 3}}}`: The gradient of `mps`, which shares the same shape as `mps`.
- `τ::T`: The updating rate.

# Description
This function calculates the projected gradient on the tangent plane and updates `mps` with respect to the updating rate `τ`.

"""
function steifel_update!(mps::LMPS, G::@NamedTuple{left::Matrix{T}, tensors::Vector{Array{T, 3}}}, τ::Float64) where T
    left_mat = mps.left
    grad_left_mat = G.left
    projected_grad_left_mat = grad_left_mat * left_mat' - left_mat * grad_left_mat'

    mps.left = inv((Diagonal(ones(T, size(projected_grad_left_mat)[end])) + τ / 2 * projected_grad_left_mat)) *
        (Diagonal(ones(T, size(projected_grad_left_mat)[end])) - τ / 2 * projected_grad_left_mat) *
        left_mat

    for i in 1:length(mps.tensors)
        mps_mat = reshape(mps.tensors[i], size(mps.tensors[i], 1) * size(mps.tensors[i], 2), size(mps.tensors[i], 3))
        grad_mps_mat = reshape(G.tensors[i], size(G.tensors[i], 1) * size(G.tensors[i], 2), size(G.tensors[i], 3))
        projected_grad_mat = grad_mps_mat * mps_mat' - mps_mat * grad_mps_mat'

        Y_mat = inv((Diagonal(ones(T, size(projected_grad_mat)[end])) + τ / 2 * projected_grad_mat)) *
            (Diagonal(ones(T, size(projected_grad_mat)[end])) - τ / 2 * projected_grad_mat) *
            mps_mat

        mps.tensors[i] = reshape(Y_mat, size(mps.tensors[i]))
    end
end

function transverse_Ising_1d_exact_GSE(L::Int64)
    return 1 - csc(pi/(2*(2*L+1)))
end