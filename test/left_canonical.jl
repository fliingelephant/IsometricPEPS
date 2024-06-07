using IsometricPEPS

using LinearAlgebra
using OMEinsum
using Test
using Zygote

@testset "left canonical" begin
    D = 128;
    n = 10;
    χ = 10;

    mps = LMPS(rand(ComplexF64, 2, 2),
        push!([rand(ComplexF64, log2(D) < i ? D : 2^i, 2, log2(D) < i+1 ? D : 2^(i+1)) for i in 1:n],
        rand(ComplexF64, D, 2, χ))
    );
    @test length(mps.tensors) == n + 1

    isometric_mps = left_canonicalize!(mps);
    @test length(isometric_mps.tensors) == n
    @test size(isometric_mps.tensors[end], 3) == 2 * χ

    contracted = ein"ia, ib->ab"(isometric_mps.left, conj(isometric_mps.left));
    @test contracted ≈ Diagonal(ones(ComplexF64, size(contracted, 1)))

    for i = 1:length(isometric_mps.tensors)
        contracted = ein"(ab, aic), bid->cd"(contracted, isometric_mps.tensors[i], conj(isometric_mps.tensors[i]));
        @test contracted ≈ Diagonal(ones(ComplexF64, size(contracted, 1)))
    end

    for _ in 1:5
        G = Zygote.gradient(x -> begin
            energy = transverse_Ising_1d_energy(x, 1.0, 1.0);
    
            return energy
        end, isometric_mps)[1];
    
        steifel_update!(isometric_mps, G, tau);
    end

    contracted = ein"ia, ib->ab"(isometric_mps.left, conj(isometric_mps.left));
    @test contracted ≈ Diagonal(ones(ComplexF64, size(contracted, 1)))

    for i = 1:length(isometric_mps.tensors)
        contracted = ein"(ab, aic), bid->cd"(contracted, isometric_mps.tensors[i], conj(isometric_mps.tensors[i]));
        @test contracted ≈ Diagonal(ones(ComplexF64, size(contracted, 1)))
    end
end