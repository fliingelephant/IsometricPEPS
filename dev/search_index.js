var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = IsometricPEPS","category":"page"},{"location":"#IsometricPEPS","page":"Home","title":"IsometricPEPS","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for IsometricPEPS.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [IsometricPEPS]","category":"page"},{"location":"#IsometricPEPS.LMPS","page":"Home","title":"IsometricPEPS.LMPS","text":"LMPS(left::AbstractTensor, tensors::Vector{AbstractTensor})\n\nLMPS represents a Matrix Product State (MPS) with a left boundary tensor of order 2 (left) and a vector of tensors of order 3 (tensors).\n\nArguments\n\nleft::AbstractTensor: The left boundary tensor of order 2.\ntensors::Vector{AbstractTensor}: The vector of tensors of order 3.\n\n\n\n\n\n","category":"type"},{"location":"#IsometricPEPS.left_canonicalize!","page":"Home","title":"IsometricPEPS.left_canonicalize!","text":"left_canonicalize!(mps::LMPS, atol::Real=0)\n\nPerform left canonicalization on an LMPS.\n\nArguments\n\nmps::LMPS: The LMPS to be left canonicalized.\natol::Real=0: Truncation tolerance for performing SVD.\n\nOutput\n\nAn LMPS which consists of an isometric left boundary tensor of order 2 and (n-1) left canonical tensors of order 3.\n\n\n\n\n\n","category":"function"},{"location":"#IsometricPEPS.steifel_update!-Union{Tuple{T}, Tuple{LMPS, @NamedTuple{left::Matrix{T}, tensors::Array{Array{T, 3}, 1}}, Float64}} where T","page":"Home","title":"IsometricPEPS.steifel_update!","text":"steifel_update!(mps::LMPS, G::@NamedTuple{left::Matrix{T}, tensors::Vector{Array{T, 3}}}, τ::T) where T\n\nUpdate the LMPS mps using the gradient G and the updating rate τ. The gradient G should have the same shape as mps.\n\nArguments\n\nmps::LMPS: The LMPS to update.\nG::@NamedTuple{left::Matrix{T}, tensors::Vector{Array{T, 3}}}: The gradient of mps, which shares the same shape as mps.\nτ::T: The updating rate.\n\nDescription\n\nThis function calculates the projected gradient on the tangent plane and updates mps with respect to the updating rate τ.\n\n\n\n\n\n","category":"method"},{"location":"#IsometricPEPS.transverse_Ising_1d_effH-Tuple{LMPS, Float64, Float64}","page":"Home","title":"IsometricPEPS.transverse_Ising_1d_effH","text":"transverse_Ising_1d_effH(mps::LMPS, h::Float64, J::Float64)\n\nThis function generates the effective Hamiltonian for a 1d transverse Ising model by contracting ZZ and X terms together for every two sites.\n\nArguments\n\nmps::LMPS: The LMPS representing the state of the system.\nJ::Float64: The coupling strength between adjacent spins.\nh::Float64: The strength of the transverse magnetic field.\n\nReturns\n\neff_H::Matrix{ComplexF64}: The effective Hamiltonian matrix.\n\n\n\n\n\n","category":"method"},{"location":"#IsometricPEPS.transverse_Ising_1d_seperate_effH-Tuple{LMPS, Float64, Float64}","page":"Home","title":"IsometricPEPS.transverse_Ising_1d_seperate_effH","text":"transverse_Ising_1d_seperate_effH(mps::LMPS, h::Float64, J::Float64)\n\nThis function generates the effective Hamiltonian for a 1d transverse Ising model by contracting ZZ and X terms seperately.\n\nArguments\n\nmps::LMPS: The LMPS representing the state of the system.\nJ::Float64: The coupling strength between adjacent spins.\nh::Float64: The strength of the transverse magnetic field.\n\nReturns\n\neff_H::Matrix{ComplexF64}: The effective Hamiltonian matrix.\n\n\n\n\n\n","category":"method"},{"location":"#IsometricPEPS.truncated_svd","page":"Home","title":"IsometricPEPS.truncated_svd","text":"truncated_svd(tmat::AbstractMatrix, dmax::Int, atol::Real=0)\n\nPerform singular value decomposition (SVD) for tmat and keep only the largest dmax singular values (not less than atol) and associated columns and rows.\n\nArguments\n\ntmat::AbstractMatrix: The matrix to perform SVD on.\ndmax::Int: The maximum number of singular values to keep.\natol::Real=0: The minimum threshold for singular values.\n\nReturns\n\nU: The left singular vectors.\nS: The singular values.\nVt: The right singular vectors.\nsum_discarded: The sum of discarded singular values.\n\n\n\n\n\n","category":"function"}]
}