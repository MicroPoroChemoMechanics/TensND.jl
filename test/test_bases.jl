@testsection "Bases" begin
    for T in (Sym, Float64)
        @testsection "type $T" begin
            sq2 = √T(2)
            v = T[1 0 0; 0 1 0; 0 1 1]
            bv = Basis(v)
            @test eltype(bv) == T
            @test vecbasis(bv) == T[1 0 0; 0 1 0; 0 1 1]
            @test vecbasis(bv, :cont) == T[1 0 0; 0 1 -1; 0 0 1]
            @test metric(bv) == T[1 0 0; 0 2 1; 0 1 1]
            @test metric(bv, :cont) == bv.gⁱʲ == T[1 0 0; 0 1 -1; 0 -1 2]
            @test isone(vecbasis(bv, :cont)' * vecbasis(bv, :cov))
            @test isone(metric(bv, :cont)' * metric(bv, :cov))
            @test !isorthogonal(bv)

            bw = Basis(vecbasis(bv, :cont), :cont)
            @test eltype(bw) == T
            @test vecbasis(bw) == T[1 0 0; 0 1 0; 0 1 1]
            @test vecbasis(bw, :cont) == T[1 0 0; 0 1 -1; 0 0 1]
            @test metric(bw) == T[1 0 0; 0 2 1; 0 1 1]
            @test metric(bw, :cont) == T[1 0 0; 0 1 -1; 0 -1 2]
            @test isone(vecbasis(bw, :cont)' * vecbasis(bw, :cov))
            @test isone(metric(bw, :cont)' * metric(bw, :cov))
            @test !isorthogonal(bw)

            nb = normalize(bv)
            @test opequal(vecbasis(nb), vecbasis(nb, :cov))
            @test opequal(vecbasis(nb), T[1 0 0; 0 sq2 / 2 0; 0 sq2 / 2 1])
            θ, ϕ, ψ = T == Sym ? symbols("θ, ϕ, ψ", real = true) : rand(T, 3)
            br = Basis(θ, ϕ, ψ)
            @test br isa RotatedBasis
            @test opequal(metric(br), 1I)
            @test opequal(vecbasis(br, :cont), vecbasis(br, :cov))
            if T == Sym
                @test angles(br) == (θ = θ, ϕ = ϕ, ψ = ψ)
            else
                angbr = angles(br)
                @test angbr.θ ≈ θ && angbr.ϕ ≈ ϕ && angbr.ψ ≈ ψ
                angbr = angles(Array(br))
                @test angbr.θ ≈ θ && angbr.ϕ ≈ ϕ && angbr.ψ ≈ ψ
            end
            λ = T == Sym ? [symbols("λ$i", real = true) for i in 1:get_dim(br)] : rand(get_dim(br))
            e = vecbasis(br, :cov) .* λ'
            bo = Basis(e)
            @test bo isa OrthogonalBasis
            @test bo isa AllOrthogonalBasis

            br = Basis(θ)
            @test opequal(metric(br), 1I)
            @test opequal(vecbasis(br, :cont), vecbasis(br, :cov))
            if T == Sym
                @test angles(br) == (θ = θ,)
            else
                @test angles(br).θ ≈ θ
                @test angles(Array(br)).θ ≈ θ
            end
            λ = T == Sym ? [symbols("λ$i", real = true) for i in 1:get_dim(br)] : rand(get_dim(br))
            e = vecbasis(br, :cov) .* λ'
            bo = Basis(e)
            @test bo isa OrthogonalBasis
            @test bo isa AllOrthogonalBasis

            b2 = Basis(zero(T))
            @test b2 isa CanonicalBasis
            @test get_dim(b2) == 2
            b3 = Basis(zero(T), zero(T), zero(T))
            @test b3 isa CanonicalBasis
            @test get_dim(b3) == 3

            if T ≠ Sym
                # General (non-orthogonal) basis built from a random invertible
                # matrix.  Note the duality relations are only satisfied up to
                # rounding here, unlike the exact integer basis `bv` above, so
                # compare with `≈` rather than `isone`.
                v = rand(T, 3, 3)
                while det(v) ≈ 0
                    v = rand(T, 3, 3)
                end
                bvr = Basis(v)
                @test eltype(bvr) == T
                @test !isorthogonal(bvr)
                @test vecbasis(bvr, :cont)' * vecbasis(bvr, :cov) ≈ 1I
                @test metric(bvr, :cont)' * metric(bvr, :cov) ≈ 1I
                @test Array(vecbasis(bvr, :cov)) ≈ v
            end

        end
    end
end

@testsection "Bases from integer angles" begin
    # A rotated basis took the element type of its angles, which integer
    # angles cannot hold — their cosines are floating-point — so building one
    # threw `InexactError`. It now takes the type of the rotation.
    @test RotatedBasis(1, 0, 0) isa RotatedBasis{3, Float64}
    @test vecbasis(RotatedBasis(1, 2, 3), :cov) == vecbasis(RotatedBasis(1.0, 2.0, 3.0), :cov)
    @test RotatedBasis(1) isa RotatedBasis{2, Float64}
    @test vecbasis(RotatedBasis(1), :cov) == vecbasis(RotatedBasis(1.0), :cov)
    @test RotatedBasis(0, 0, 0) isa CanonicalBasis
end

@testsection "A variance is :cov or :cont" begin
    ℬ = RotatedBasis(0.3, 0.2, 0.1)
    @test vecbasis(ℬ, :cov) == vecbasis(ℬ, Val(:cov))
    @test_throws ArgumentError vecbasis(ℬ, :covariant)
    @test_throws ArgumentError metric(ℬ, :foo)
    # Every function taking a variance as a `Symbol` reads it through the same
    # bridge, so a wrong one says so instead of failing to find a `Val` method.
    @test Basis(Matrix(ℬ), :cont) == Basis(Matrix(ℬ), Val(:cont))
    @test_throws ArgumentError Basis(ℬ, :foo)
    @test tens_basis(ℬ, 2, :cont) == tens_basis(ℬ, 2, Val(:cont))
    @test tens_basis(ℬ, :cont) == tens_basis(ℬ, Val(:cont))
    @test length(tens_basis(ℬ)) == 3
    @test_throws ArgumentError tens_basis(ℬ, :foo)
    @test TensND.strvecbasis(ℬ, 2, :cont) == "𝐞²"
    @test TensND.strvecbasis(ℬ, "θ", :cov; vec = "𝐚") == "𝐚" * TensND.subscriptchar("θ")
    @test_throws ArgumentError TensND.strvecbasis(ℬ, 1, :foo)
end
