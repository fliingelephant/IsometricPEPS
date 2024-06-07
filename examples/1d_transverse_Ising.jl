using IsometricPEPS

D = 128;
n = 10;
χ = 10;
tau = 1e-2;

mps = LMPS(rand(ComplexF64, 2, 2),
        push!([rand(ComplexF64, log2(D) < i ? D : 2^i, 2, log2(D) < i+1 ? D : 2^(i+1)) for i in 1:n],
        rand(ComplexF64, D, 2, χ))
    );

@time left_canonicalize!(mps);

