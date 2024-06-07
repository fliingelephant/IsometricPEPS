using IsometricPEPS

using Test

D = 128;
n = 10;
χ = 10;

@testset "evaluation" begin
    mps = LMPS(rand(ComplexF64, 2, 2),
        push!([rand(ComplexF64, log2(D) < i ? D : 2^i, 2, log2(D) < i+1 ? D : 2^(i+1)) for i in 1:n],
        rand(ComplexF64, D, 2, χ))
    );

    left_canonicalize!(mps);
    @test transverse_Ising_1d_seperate_effH(mps, 1.0, 1.0) ≈ transverse_Ising_1d_effH(mps, 1.0, 1.0)
    @test transverse_Ising_1d_seperate_energy(mps, 1.0, 1.0) ≈ transverse_Ising_1d_energy(mps, 1.0, 1.0)

    @test transverse_Ising_1d_seperate_effH(mps, 2.0, 3.5) ≈ transverse_Ising_1d_effH(mps, 2.0, 3.5)
    @test transverse_Ising_1d_seperate_energy(mps, 2.0, 3.5) ≈ transverse_Ising_1d_energy(mps, 2.0, 3.5)
end