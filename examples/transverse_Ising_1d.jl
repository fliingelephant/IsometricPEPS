using IsometricPEPS

using LinearAlgebra
using Test
using Zygote

D = 128;    # max bond dimension
n = 10;    # number of bulk tensors
χ = 10;   # dimension of the rightmost leg when generating the MPS
tau = 5 * 1e-3;   # updating rate

# generate a random MPS with a left boundary tensor of order 2, n tensors of order 3, and a right boundary tensor of order 3
mps = LMPS(rand(ComplexF64, 2, 2),
    push!([rand(ComplexF64, log2(D) < i ? D : 2^i, 2, log2(D) < i+1 ? D : 2^(i+1)) for i in 1:n],
    rand(ComplexF64, D, 2, χ))
);

# left canonicalize the MPS and discard the rightmost tensor
# a left boundary tensor of order 2, n tensors of order 3, corresponding to an (n+1)-site OBC system
# 2 * χ: dimension of the rightmost leg, corresponding to the projected subspace
left_canonicalize!(mps);

# solution for 1d transverse Ising with 11 sites at the critical point, obtained by ED
exact_sepctrum_11_sites = [-13.65364354,
    -13.38067389,
    -12.83981949,
    -12.56684984,
    -12.31412509,
    -12.04115544, 
    -11.81338339, 
    -11.54041374, 
    -11.50030104, 
    -11.34692226, 
    -11.22733139, 
    -11.0739526, 
    -10.99955934, 
    -10.92343097, 
    -10.72658969, 
    -10.65046132, 
    -10.55079838, 
    -10.5330982, 
    -10.47386494, 
    -10.27782873
];

# exact solution for 1d transverse Ising with 11 sites at the critical point
@test transverse_Ising_1d_exact_GSE(11) ≈ exact_sepctrum_11_sites[1]

min_trace = sum(exact_sepctrum_11_sites);

for _ in 1:1000
    G = Zygote.gradient(x -> begin
        energy = transverse_Ising_1d_energy(x, 1.0, 1.0);
        @show energy

        # non-physical
        if (energy < min_trace)
            println("WRONG")
        end

        return energy
    end, mps)[1];

    steifel_update!(mps, G, tau);
end
