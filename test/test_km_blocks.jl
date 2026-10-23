using TensND.StaticArrays: SMatrix, SVector
# Kelvin-Mandel block algebra about a normal (`src/km_blocks.jl`): the partition,
# the acoustic tensor on a triclinic stiffness (the only case that sees the
# (3, 2, 1) slot reversal), the plane pseudo-inverse against the SVD, the Schur
# inverse against `inv`, then duals and symbols.
@testsection "Kelvin-Mandel blocks about a normal" begin
    rng = Xoshiro(20261026)
    rand_km6(rng) = (A = randn(rng, 6, 6); SMatrix{6, 6}(A * transpose(A) + 8I))
    iso_km(k, μ) = SMatrix{6, 6}(KM(TensISO{3}(3k, 2μ)))

    @testset "partition, blocks and embeddings" begin
        @test KM_IP == SVector(1, 2, 6) && KM_OP == SVector(3, 4, 5)
        M = rand_km6(rng)
        @test km_op_block(km_op_embed(km_op_block(M))) == km_op_block(M)
        @test km_ip_block(km_ip_embed(km_ip_block(M))) == km_ip_block(M)
        @test iszero(km_ip_block(km_op_embed(km_op_block(M))))
        @test iszero(km_op_block(km_ip_embed(km_ip_block(M))))
        @test km_from_blocks(km_ip_block(M), M[KM_IP, KM_OP], M[KM_OP, KM_IP], km_op_block(M)) == M
    end

    @testset "acoustic tensor: the (3,2,1) reversal on a triclinic stiffness" begin
        M = rand_km6(rng)
        C = inv_KM(Matrix(M))
        arr = Array(get_array(C))
        K_naive = SMatrix{3, 3}(arr[3, :, 3, :])                   # Kₐᵦ = C₃ₐ₃ᵦ
        @test acoustic_tensor(M) ≈ K_naive atol = 1.0e-12
        @test !isapprox(K_naive, Diagonal(K_naive); atol = 1.0e-6)  # the reversal is exercised
        @test Array(get_array(acoustic_tensor(C, (0, 0, 1)))) ≈ K_naive atol = 1.0e-12
        # Any direction: the tensor form against the frame form.
        n = normalize(randn(rng, 3))
        ℬ = frame_from_axis(n)
        Kn = acoustic_tensor(SMatrix{6, 6}(KM(C, ℬ)))
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
        @test iszero(km_ip_block(P))
        M = rand_km6(rng)
        Mop = km_op_embed(km_op_block(M))
        @test plane_pinv(M) ≈ pinv(Matrix(Mop)) atol = 1.0e-10   # against the SVD
        @test plane_pinv(plane_pinv(M)) ≈ Mop atol = 1.0e-10
        Pm = plane_pinv(M)
        @test Pm * M * Pm ≈ Pm atol = 1.0e-10                     # ℙ:ℂ:ℙ = ℙ
        K3 = [4.0 0.3 0.2; 0.3 3.0 0.1; 0.2 0.1 2.0]
        @test plane_pinv2(K3) ≈ [0 0 0; 0 0 0; 0 0 1 / 2.0]
        # compliance_op_block undoes the acoustic map on the compliance side.
        X = [1.0 0.2 0.1; 0.2 2.0 0.3; 0.1 0.3 3.0]
        Sx = km_op_embed(compliance_op_block(X))
        Ka = acoustic_tensor(plane_pinv(Sx))
        @test Ka ≈ inv(X) atol = 1.0e-12
    end

    @testset "Schur inverse" begin
        for _ in 1:5
            M = rand_km6(rng)
            @test km_inv_schur(M) * M ≈ I atol = 1.0e-10
            @test km_inv_schur(M) ≈ inv(Matrix(M)) atol = 1.0e-10
        end
    end

    @testset "duals and symbols" begin
        M0 = rand_km6(rng)
        dM = rand_km6(rng)
        # d(M⁻¹) = −M⁻¹ dM M⁻¹, the analytic derivative of the inverse.
        Mi = inv(Matrix(M0))
        @test ForwardDiff.derivative(t -> km_inv_schur(M0 + t * dM), 0.0) ≈ -Mi * dM * Mi atol = 1.0e-10
        Bop = km_op_block(M0)
        @test ForwardDiff.derivative(t -> plane_pinv(M0 + t * dM), 0.0) ≈
            km_op_embed(-inv(Bop) * km_op_block(dM) * inv(Bop)) atol = 1.0e-10
        @syms k::positive μ::positive
        Ks = acoustic_tensor(SMatrix{6, 6}(KM(TensISO{3}(3k, 2μ))))
        @test simplify(Ks[3, 3] - (k + 4μ / 3)) == 0 && simplify(Ks[1, 1] - μ) == 0
        Ps = plane_pinv(SMatrix{6, 6}(KM(TensISO{3}(3k, 2μ))))
        @test simplify(Ps[4, 4] - 1 / (2μ)) == 0
        TensND.Symbolics.@variables a b
        Mn = km_from_blocks(SMatrix{3, 3}(a * I), zeros(3, 3), zeros(3, 3), SMatrix{3, 3}(b * I))
        sub(x) = Float64(TensND.Symbolics.value(TensND.Symbolics.substitute(x, Dict(a => 3.0, b => 2.0); fold = Val(true))))
        @test sub(km_inv_schur(Mn)[3, 3]) ≈ 1 / 2.0 && sub(km_inv_schur(Mn)[1, 1]) ≈ 1 / 3.0
    end
end
