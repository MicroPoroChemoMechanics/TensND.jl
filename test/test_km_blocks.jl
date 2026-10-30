using TensND.StaticArrays: SMatrix, SVector
# Kelvin-Mandel block algebra about a normal (`src/km_blocks.jl`): the partition,
# the acoustic tensor on a triclinic stiffness (the only case that sees the
# (3, 2, 1) slot reversal), the plane pseudo-inverse against the SVD, the Schur
# inverse against `inv`, then duals and symbols.
@testsection "Kelvin-Mandel blocks about a normal" begin
    rng = Xoshiro(20261026)
    rand_km6(rng) = (A = randn(rng, 6, 6); SMatrix{6, 6}(A * transpose(A) + 8I))
    iso_km(k, μ) = SMatrix{6, 6}(mandel(TensISO{3}(3k, 2μ)))

    @testset "partition, blocks and embeddings" begin
        @test MANDEL_IP == SVector(1, 2, 6) && MANDEL_OP == SVector(3, 4, 5)
        M = rand_km6(rng)
        @test mandel_op_block(mandel_op_embed(mandel_op_block(M))) == mandel_op_block(M)
        @test mandel_ip_block(mandel_ip_embed(mandel_ip_block(M))) == mandel_ip_block(M)
        @test iszero(mandel_ip_block(mandel_op_embed(mandel_op_block(M))))
        @test iszero(mandel_op_block(mandel_ip_embed(mandel_ip_block(M))))
        @test mandel_from_blocks(mandel_ip_block(M), M[MANDEL_IP, MANDEL_OP], M[MANDEL_OP, MANDEL_IP], mandel_op_block(M)) == M
    end

    @testset "acoustic tensor: the (3,2,1) reversal on a triclinic stiffness" begin
        M = rand_km6(rng)
        C = from_mandel(Matrix(M))
        arr = Array(get_array(C))
        K_naive = SMatrix{3, 3}(arr[3, :, 3, :])                   # Kₐᵦ = C₃ₐ₃ᵦ
        @test acoustic_tensor(M) ≈ K_naive atol = 1.0e-12
        @test !isapprox(K_naive, Diagonal(K_naive); atol = 1.0e-6)  # the reversal is exercised
        @test Array(get_array(acoustic_tensor(C, (0, 0, 1)))) ≈ K_naive atol = 1.0e-12
        # Any direction: the tensor form against the frame form.
        n = collect(random_axis(rng))
        ℬ = frame_from_axis(n)
        Kn = acoustic_tensor(SMatrix{6, 6}(mandel(C, ℬ)))
        @test Array(get_array(change_tens(acoustic_tensor(C, n), ℬ))) ≈ Kn atol = 1.0e-10
        k, μ = 2.0, 0.8
        λ = k - 2μ / 3
        @test acoustic_tensor(iso_km(k, μ)) ≈ Diagonal([μ, μ, λ + 2μ]) atol = 1.0e-12
    end

    @testset "plane pseudo-inverses and the flat Hill tensor" begin
        k, μ = 2.0, 0.8
        λ = k - 2μ / 3
        P = plane_pinv(iso_km(k, μ))
        # The flat-inclusion Hill tensor of an isotropic medium, normal e₃.
        @test P[3, 3] ≈ 1 / (λ + 2μ) && P[4, 4] ≈ 1 / (2μ) && P[5, 5] ≈ 1 / (2μ)
        @test iszero(mandel_ip_block(P))
        M = rand_km6(rng)
        Mop = mandel_op_embed(mandel_op_block(M))
        @test plane_pinv(M) ≈ pinv(Matrix(Mop)) atol = 1.0e-10   # against the SVD
        @test plane_pinv(plane_pinv(M)) ≈ Mop atol = 1.0e-10
        Pm = plane_pinv(M)
        @test Pm * M * Pm ≈ Pm atol = 1.0e-10                     # ℙ:ℂ:ℙ = ℙ
        K3 = [4.0 0.3 0.2; 0.3 3.0 0.1; 0.2 0.1 2.0]
        @test plane_pinv2(K3) ≈ [0 0 0; 0 0 0; 0 0 1 / 2.0]
        # compliance_op_block undoes the acoustic map on the compliance side.
        X = [1.0 0.2 0.1; 0.2 2.0 0.3; 0.1 0.3 3.0]
        Sx = mandel_op_embed(compliance_op_block(X))
        Ka = acoustic_tensor(plane_pinv(Sx))
        @test Ka ≈ inv(X) atol = 1.0e-12
        # And it is n ⊗ˢ X ⊗ˢ n, X symmetric or not, in a random frame whose
        # third axis is n: against the symmetrized products of the components.
        for Y in (X, [1.0 0.7 -0.2; 0.1 2.0 0.4; 0.5 -0.3 3.0])
            ℬ = random_basis(rng)
            n̂ = TensND.tens_basis(ℬ, 3)
            Yt = Tens(Y, ℬ)
            H = from_mandel(mandel_op_embed(compliance_op_block(Y)), ℬ)
            @test Array(components_canon(H)) ≈ Array(components_canon(n̂ ⊗ˢ Yt ⊗ˢ n̂)) atol = 1.0e-14
        end
    end

    @testset "Schur inverse" begin
        for _ in 1:5
            M = rand_km6(rng)
            @test mandel_inv_schur(M) * M ≈ I atol = 1.0e-10
            @test mandel_inv_schur(M) ≈ inv(Matrix(M)) atol = 1.0e-10
        end
    end

    @testset "duals and symbols" begin
        M0 = rand_km6(rng)
        dM = rand_km6(rng)
        # d(M⁻¹) = −M⁻¹ dM M⁻¹, the analytic derivative of the inverse.
        Mi = inv(Matrix(M0))
        @test ForwardDiff.derivative(t -> mandel_inv_schur(M0 + t * dM), 0.0) ≈ -Mi * dM * Mi atol = 1.0e-10
        Bop = mandel_op_block(M0)
        @test ForwardDiff.derivative(t -> plane_pinv(M0 + t * dM), 0.0) ≈
            mandel_op_embed(-inv(Bop) * mandel_op_block(dM) * inv(Bop)) atol = 1.0e-10
        @syms k::positive μ::positive
        Ks = acoustic_tensor(SMatrix{6, 6}(mandel(TensISO{3}(3k, 2μ))))
        @test simplify(Ks[3, 3] - (k + 4μ / 3)) == 0 && simplify(Ks[1, 1] - μ) == 0
        Ps = plane_pinv(SMatrix{6, 6}(mandel(TensISO{3}(3k, 2μ))))
        @test simplify(Ps[4, 4] - 1 / (2μ)) == 0
        TensND.Symbolics.@variables a b
        Mn = mandel_from_blocks(SMatrix{3, 3}(a * I), zeros(3, 3), zeros(3, 3), SMatrix{3, 3}(b * I))
        sub(x) = Float64(TensND.Symbolics.value(TensND.Symbolics.substitute(x, Dict(a => 3.0, b => 2.0); fold = Val(true))))
        @test sub(mandel_inv_schur(Mn)[3, 3]) ≈ 1 / 2.0 && sub(mandel_inv_schur(Mn)[1, 1]) ≈ 1 / 3.0
    end

    @testset "deviatoric pseudo-inverse" begin
        J, K = SMatrix{6, 6}(mandel(TensISO{3}(1, 0))), SMatrix{6, 6}(mandel(TensISO{3}(0, 1)))
        # The deviatoric stiffness of an incompressible isotropic medium.
        @test deviatoric_pinv(2.0 * 0.8 * K) ≈ K / (2 * 0.8) atol = 1.0e-14
        # A random anisotropic operator on deviators: the four Moore-Penrose
        # conditions, the projector onto the deviators as the identity there.
        for _ in 1:5
            M = K * rand_km6(rng) * K
            P = deviatoric_pinv(M)
            @test M * P ≈ K atol = 1.0e-10
            @test P * M ≈ K atol = 1.0e-10
            @test M * P * M ≈ M atol = 1.0e-10
            @test P * M * P ≈ P atol = 1.0e-10
            @test maximum(abs, P * J) < 1.0e-12
            @test P ≈ pinv(Matrix(M)) atol = 1.0e-10
        end
        # Derivative: d(M⁺) = −M⁺ dM M⁺ on deviators, as for an inverse.
        M0, dM = K * rand_km6(rng) * K, K * rand_km6(rng) * K
        P0 = deviatoric_pinv(M0)
        @test ForwardDiff.derivative(t -> deviatoric_pinv(M0 + t * dM), 0.0) ≈ -P0 * dM * P0 atol = 1.0e-10
        # Symbols: the incompressible isotropic case in closed form.
        @syms μs::positive
        Ps = deviatoric_pinv(SMatrix{6, 6}(mandel(TensISO{3}(0, 2μs))))
        @test simplify(Ps[1, 1] - 2 / (3 * 2μs)) == 0 && simplify(Ps[4, 4] - 1 / (2μs)) == 0 && simplify(Ps[1, 2] + 1 / (3 * 2μs)) == 0
    end
end

# The 3×3 inverses by cofactors and the acoustic kernel, against the generic
# inverse and the tensor contraction, with derivatives and symbols.
@testset "inv3, sym3_inv, acoustic_components" begin
    rng = Xoshiro(20261030)
    for _ in 1:50
        K = randn(rng, 3, 3) + 3I
        @test inv3(K) ≈ inv(K)
        @test inv3(K) isa SMatrix{3, 3, Float64}
        S = K + transpose(K)
        @test collect(sym3_inv(S[1, 1], S[2, 2], S[3, 3], S[1, 2], S[1, 3], S[2, 3])) ≈ vec(inv(S))
        ℂ = rand(rng, Tens{4, 3}; definite = true)
        ζ = randn(rng, 3)                          # a direction of any length
        A = acoustic_components(Array(components_canon(ℂ)), ζ)
        Kζ = components_canon(acoustic_tensor(ℂ, ζ))
        @test collect(A) ≈ [Kζ[1, 1], Kζ[2, 2], Kζ[3, 3], Kζ[1, 2], Kζ[1, 3], Kζ[2, 3]]
    end
    # Derivatives: d(K⁻¹)/dx = -K⁻¹ (dK/dx) K⁻¹ for K(x) = K₀ + x E.
    K0 = SMatrix{3, 3}(4.0, 1.0, 0.5, 1.0, 3.0, 0.2, 0.5, 0.2, 2.0)
    E = SMatrix{3, 3}(1.0, 0.3, 0.0, 0.3, 0.0, 0.1, 0.0, 0.1, 0.5)
    @test ForwardDiff.derivative(x -> inv3(K0 + x * E), 0.0) ≈ -inv(K0) * E * inv(K0)
    dS = ForwardDiff.derivative(x -> collect(sym3_inv((K0 + x * E)[[1, 5, 9, 4, 7, 8]]...)), 0.0)
    @test dS ≈ vec(-inv(K0) * E * inv(K0))
    # Symbols: the cofactor inverse of a symbolic matrix times the matrix is the identity.
    Symbolics.@variables a b c d e f
    Ks = SMatrix{3, 3}(a, d, e, d, b, f, e, f, c)
    vals = Dict(a => 4.0, b => 3.0, c => 2.0, d => 1.0, e => 0.5, f => 0.2)
    num(x) = Symbolics.value(Symbolics.substitute(x, vals; fold = Val(true)))
    @test num.(inv3(Ks)) ≈ inv(K0)
    @test num.(collect(sym3_inv(a, b, c, d, e, f))) ≈ vec(inv(K0))
    @test (@inferred sym3_inv(K0[1, 1], K0[2, 2], K0[3, 3], K0[1, 2], K0[1, 3], K0[2, 3])) isa NTuple{9, Float64}
end
