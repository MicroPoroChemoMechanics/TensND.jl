# `KM`, `_static_tensor` and `change_tens` of the fourth-order structured tensors
# in closed form: the 36 Kelvin-Mandel entries, the congruence of the material
# frame, the class's own formula with its axis or frame written in the target
# basis. Each is checked against the route it replaced — the 81 canonical
# components, rotated by the generic change of basis — in random bases, right-
# and left-handed, on floats, duals and symbols.
@testsection "Structured tensors: Kelvin-Mandel and change of basis in closed form" begin
    rng = Xoshiro(20261022)
    _rot(rng) = (Q = Matrix(qr(randn(rng, 3, 3)).Q); det(Q) < 0 ? Q * Diagonal([1.0, 1.0, -1.0]) : Q)
    _unit(v) = Tuple(v / norm(v))
    _mirror = Diagonal([1.0, 1.0, -1.0])
    # The replaced routes, kept here as the reference.
    _old_KM(t) = tomandel(TensND._to_tensor(get_array(t), Val(3), Val(4)))
    _old_change(t, ℬ) = change_tens(TensND._as_canonical_tens(t), ℬ)
    _structured(rng) = (
        TensTI{4}(10.0, 7.0, 3.0, 4.0, 5.0, _unit(randn(rng, 3))),
        TensTI{4}(10.0, 7.0, 3.0, 2.5, 4.0, 5.0, _unit(randn(rng, 3))),
        TensND.TensTI{4, Float64, 8}((10.0, 7.0, 3.0, 2.5, 4.0, 5.0, 0.7, 0.3), _unit(randn(rng, 3))),
        TensOrtho(160.0, 120.0, 95.0, 55.0, 48.0, 42.0, 33.0, 29.0, 27.0, RotatedBasis(_rot(rng))),
        TensCubic(60.0, 16.0, 10.0, RotatedBasis(_rot(rng))),
    )

    @testset "against the 81 canonical components, random bases" begin
        for _ in 1:5, t in _structured(rng)
            K = KM(t)
            @test K isa Matrix{Float64}
            @test K ≈ _old_KM(t) rtol = 1.0e-14
            @test TensND._static_tensor(t) isa SymmetricTensor{4, 3, Float64}
            @test Array(TensND._static_tensor(t)) ≈ get_array(t) rtol = 1.0e-14
            R = _rot(rng)
            for ℬ in (RotatedBasis(R), RotatedBasis(R * _mirror), CanonicalBasis{3, Float64}())
                new, old = change_tens(t, ℬ), _old_change(t, ℬ)
                @test typeof(new) == typeof(old)
                @test TensND.get_var(new) == TensND.get_var(old)
                @test Array(get_array(new)) ≈ Array(get_array(old)) rtol = 1.0e-14
                # And the tensor it denotes is the one it was.
                @test KM(new, CanonicalBasis{3, Float64}()) ≈ K rtol = 1.0e-13
            end
        end
    end

    @testset "duals: coefficients, frame and basis" begin
        F = RotatedBasis(_rot(rng))
        R = _rot(rng) * _mirror
        c = (160.0, 120.0, 95.0, 55.0, 48.0, 42.0, 33.0, 29.0, 27.0)
        e₁ = ntuple(i -> Float64(i == 1), 9)
        # Linear in its coefficients: the derivative along C₁₁ is the tensor
        # with C₁₁ = 1 and the others zero.
        d = ForwardDiff.derivative(x -> KM(TensOrtho(x, c[2:end]..., F)), 160.0)
        @test d ≈ KM(TensOrtho(e₁..., F)) rtol = 1.0e-14
        ℬd(x) = RotatedBasis(R .+ zero(x))
        d = ForwardDiff.derivative(x -> Array(get_array(change_tens(TensOrtho(x, c[2:end]..., F), ℬd(x)))), 160.0)
        @test d ≈ Array(get_array(change_tens(TensOrtho(e₁..., F), RotatedBasis(R)))) rtol = 1.0e-14
        d = ForwardDiff.derivative(x -> Array(get_array(change_tens(TensCubic(x, 16.0, 10.0, F), ℬd(x)))), 60.0)
        @test d ≈ Array(get_array(change_tens(TensCubic(1.0, 0.0, 0.0, F), RotatedBasis(R)))) rtol = 1.0e-14
        d = ForwardDiff.derivative(x -> Array(get_array(change_tens(TensTI{4}(x, 7.0, 3.0, 4.0, 5.0, (0.6, 0.0, 0.8)), ℬd(x)))), 10.0)
        @test d ≈ Array(get_array(change_tens(TensTI{4}(1.0, 0.0, 0.0, 0.0, 0.0, (0.6, 0.0, 0.8)), RotatedBasis(R)))) rtol = 1.0e-14
        # Along the orientation of the frame, the axis and the target basis:
        # the replaced route, differentiated, is the reference.
        Fθ(θ) = RotatedBasis(θ, 0.4, 1.1)
        orth(θ) = TensOrtho((c .+ zero(θ))..., Fθ(θ))
        @test ForwardDiff.derivative(θ -> KM(orth(θ)), 0.7) ≈ ForwardDiff.derivative(θ -> _old_KM(orth(θ)), 0.7) rtol = 1.0e-12
        ti(θ) = TensTI{4}(10.0, 7.0, 3.0, 2.5, 4.0, 5.0, (sin(θ), zero(θ), cos(θ)))
        @test ForwardDiff.derivative(θ -> KM(ti(θ)), 0.7) ≈ ForwardDiff.derivative(θ -> _old_KM(ti(θ)), 0.7) rtol = 1.0e-12
        for t in (orth, ti)
            new(θ) = Array(get_array(change_tens(t(θ), Fθ(2θ))))
            old(θ) = Array(get_array(_old_change(t(θ), Fθ(2θ))))
            @test ForwardDiff.derivative(new, 0.7) ≈ ForwardDiff.derivative(old, 0.7) rtol = 1.0e-12
        end
    end

    @testset "symbols: the same tensor as with numbers" begin
        R = _rot(rng) * _mirror
        F = _rot(rng)
        n = (0.6, 0.0, 0.8)
        vals = (160.0, 120.0, 95.0, 55.0, 48.0, 42.0, 33.0, 29.0, 27.0)
        # Symbolics: a symbolic coefficient takes the entry route of `KM`, the
        # expressions it always had.
        TensND.Symbolics.@variables a b
        sub(x) = Float64.(TensND.Symbolics.value.(TensND.Symbolics.substitute.(x, Ref(Dict(a => 160.0, b => 7.0)); fold = Val(true))))
        tn = TensOrtho(a, vals[2:end]..., RotatedBasis(F))
        @test sub(KM(tn)) ≈ KM(TensOrtho(vals..., RotatedBasis(F))) rtol = 1.0e-14
        tin = TensTI{4}(10.0, b, 3.0, 4.0, 5.0, n)
        @test sub(KM(tin)) ≈ KM(TensTI{4}(10.0, 7.0, 3.0, 4.0, 5.0, n)) rtol = 1.0e-14
        ℬn = RotatedBasis(TensND.Symbolics.Num.(R))
        @test sub(Array(get_array(change_tens(tin, ℬn)))) ≈
            Array(get_array(change_tens(TensTI{4}(10.0, 7.0, 3.0, 4.0, 5.0, n), RotatedBasis(R)))) rtol = 1.0e-14
        # SymPy: a symbolic axis, written in a symbolic basis, exactly.
        @syms θ::real C::positive
        ts = TensTI{4}(C, 7, 3, 4, 5, (sin(θ), 0, cos(θ)))
        ℬs = RotatedBasis(Sym.(R))
        ex = Array(get_array(change_tens(ts, ℬs)))
        num = Array(get_array(change_tens(TensTI{4}(10.0, 7.0, 3.0, 4.0, 5.0, (sin(0.7), 0.0, cos(0.7))), RotatedBasis(R))))
        @test Float64.(subs.(ex, θ => 0.7, C => 10)) ≈ num rtol = 1.0e-12
        cs = TensCubic(C, 16, 10, ℬs)
        @test Float64.(subs.(KM(cs), C => 60)) ≈ KM(TensCubic(60.0, 16.0, 10.0, RotatedBasis(R))) rtol = 1.0e-12
    end
end

@testsection "Real and imaginary parts of a tensor" begin
    # Each part is checked against the parts of the canonical components, an
    # independent reading of the same tensor, and keeps the class.
    canon(t) = Array(components_canon(t))
    ℬ = RotatedBasis(0.7, -0.4, 1.1)
    n = (0.6, 0.0, 0.8)
    for t in (
            TensISO{3}(1.0 + 2im, 3.0 - 4im),
            TensTI{4}(10.0 + 1im, 7.0, 3.0 - 2im, 4.0, 5.0 + 0.5im, n),
            TensOrtho(160.0 + 1im, 120.0, 95.0, 55.0 - 2im, 48.0, 42.0, 33.0, 29.0 + 3im, 27.0, ℬ),
            TensCubic(60.0 + 2im, 16.0, 10.0 - 1im, ℬ),
            Tens(canon(TensISO{3}(1.0 + 2im, 3.0 - 4im)) .+ 0.1im, ℬ),
        )
        for f in (real, imag)
            p = f(t)
            @test eltype(p) == Float64
            @test nameof(typeof(p)) == nameof(typeof(t))
            @test canon(p) ≈ f.(canon(t))
        end
        @test canon(real(t)) + im * canon(imag(t)) ≈ canon(t)
    end
    @test axis(real(TensTI{4}(10.0 + 1im, 7.0, 3.0, 4.0, 5.0, n))) == n
    @test frame(real(TensOrtho(1.0 + 1im, 1, 1, 0, 0, 0, 1, 1, 1, ℬ))) == ℬ
    # In a complex basis the tensor comes back in the canonical one.
    B = Basis(ComplexF64[1 1im 0; 0 1 0; 0 0 1])
    tc = Tens(ComplexF64[1 + 2im 3; 4 5 - 1im][[1, 2, 2], [1, 2, 2]], B)
    @test get_basis(real(tc)) isa CanonicalBasis
    @test canon(real(tc)) ≈ real.(canon(tc))
    # Duals and symbols go through.
    @test ForwardDiff.derivative(x -> get_data(real(TensISO{3}(x + 2im, 3.0)))[1], 1.5) == 1
    @syms a::real b::real
    @test get_data(real(TensISO{3}(a + im * b, a))) == (a, a)
end

@testsection "Coefficients: rebuild, param_norms, flat vectors" begin
    ℬ = RotatedBasis(0.7, -0.4, 1.1)
    n = (0.6, 0.0, 0.8)
    structured = (
        TensISO{3}(3.0, 2.0), TensISO{3}(2.0), TensISO{2}(3.0, 2.0),
        TensTI{4}(10.0, 7.0, 3.0, 4.0, 5.0, n), TensTI{4}(10.0, 7.0, 3.0, 2.5, 4.0, 5.0, n),
        TensND.TensTI{4, Float64, 8}((10.0, 7.0, 3.0, 2.5, 4.0, 5.0, 0.7, 0.3), n),
        TensTI{2}(1.0, 2.0, n),
        TensOrtho(160.0, 120.0, 95.0, 55.0, 48.0, 42.0, 33.0, 29.0, 27.0, ℬ),
        TensCubic(60.0, 16.0, 10.0, ℬ),
    )
    for t in structured
        d = get_data(t)
        # The same tensor from its own coefficients, and the class kept.
        @test rebuild(t, d) == t
        @test rebuild(t, collect(d)) == t
        @test from_vector(t, to_vector(t)) == t
        @test to_vector(t) == collect(d)
        # The coefficient basis is orthogonal and `param_norms` makes it an
        # isometry: the norm of the tensor from its weighted coefficients.
        w = param_norms(t)
        @test norm(t) ≈ norm(collect(w) .* collect(d))
        q = collect(1.0:length(d))
        @test norm(rebuild(t, q)) ≈ norm(collect(w) .* q)
    end
    @test all(param_norms(TensISO{3}(3.0, 2.0)) .≈ (1, sqrt(5)))
    # Mixed coefficient types are promoted; a TI axis follows its coefficients.
    td = rebuild(TensTI{4}(10.0, 7.0, 3.0, 4.0, 5.0, n), (ForwardDiff.Dual(10.0, 1.0), 7.0, 3.0, 4.0, 5.0))
    @test eltype(td) <: ForwardDiff.Dual
    @test ForwardDiff.value.(axis(td)) == n
    @test ForwardDiff.derivative(x -> get_data(rebuild(TensISO{3}(3.0, 2.0), (x, 2x)))[2], 1.0) == 2
    # Generic tensors: canonical Kelvin-Mandel coordinates, and back.
    G = Tens(SymmetricTensor{4, 3}((i, j, k, l) -> Float64(i + j + k + l) + (i == k) * (j == l)), ℬ)
    @test length(to_vector(G)) == 36
    @test param_norms(G) == ones(36)
    @test norm(G) ≈ norm(param_norms(G) .* to_vector(G))
    @test to_vector(G) == vec(KM(G, CanonicalBasis{3, Float64}()))
    @test from_vector(G, to_vector(G)) ≈ G
    S2 = Tens(SymmetricTensor{2, 3}((i, j) -> Float64(i * j)), ℬ)
    @test length(to_vector(S2)) == 6
    @test from_vector(S2, to_vector(S2)) ≈ S2
    # Symbols.
    @syms a::positive b::positive
    @test rebuild(TensISO{3}(3.0, 2.0), (a, b)) == TensISO{3}(a, b)
end

@testsection "generic_tens and change_tens_canon of structured tensors" begin
    ℬ = RotatedBasis(0.7, -0.4, 1.1)
    # A generic tensor in a rotated basis: the same tensor, in the canonical basis.
    G = Tens(SymmetricTensor{4, 3}((i, j, k, l) -> Float64(i + 2j + k * l)), ℬ)
    @test get_basis(generic_tens(G)) isa CanonicalBasis
    @test Array(get_array(generic_tens(G))) ≈ Array(components_canon(G))
    for t in (
            TensISO{3}(3.0, 2.0), TensTI{4}(10.0, 7.0, 3.0, 4.0, 5.0, (0.6, 0.0, 0.8)),
            TensOrtho(160.0, 120.0, 95.0, 55.0, 48.0, 42.0, 33.0, 29.0, 27.0, ℬ), TensCubic(60.0, 16.0, 10.0, ℬ),
        )
        @test change_tens_canon(t) === t
        g = generic_tens(t)
        @test get_basis(g) isa CanonicalBasis
        @test Array(get_array(g)) ≈ Array(get_array(t))
        z = zero(g)
        @test iszero(Array(get_array(z))) && eltype(z) == eltype(t)
        # A zero that accepts a term outside the class of `t`.
        @test z + Tens(ones(3, 3, 3, 3)) ≈ Tens(ones(3, 3, 3, 3))
    end
end

@testsection "Transverse isotropy: TI6 projection, TI2 ⋅ vector, order-2 read-off" begin
    rng = Xoshiro(20261025)
    n = Tuple(normalize(randn(rng, 3)))
    # The TI6 projection is the orthogonal projection onto the six Walpole
    # tensors: a TI6 tensor is its own projection, the residual of any
    # minor-symmetric tensor is orthogonal to the six, and an eight-parameter
    # form loses exactly ℓ₇ and ℓ₈.
    A6 = TensTI{4}(10.0, 7.0, 3.0, 2.5, 4.0, 5.0, n)
    B, d, drel = proj_tens(Val(:TI6), change_tens(Tens(Array(get_array(A6))), RotatedBasis(0.3, 0.2, 0.1)), n)
    @test B isa TensTI{4, Float64, 6}
    @test all(get_data(B) .≈ get_data(A6)) && d < 1.0e-12 && drel < 1.0e-13
    X = Tens(SymmetricTensor{4, 3}((i, j, k, l) -> randn(rng)) ⊡ SymmetricTensor{4, 3}((i, j, k, l) -> randn(rng)))
    BX, dX, _ = proj_tens(:TI6, X, n)
    R = generic_tens(BX) - X
    for i in 1:6
        Wi = generic_tens(rebuild(BX, ntuple(j -> Float64(j == i), 6)))
        @test abs(qcontract(R, Wi)) < 1.0e-10 * norm(X)
    end
    @test dX ≈ norm(R)
    A8 = TensND.TensTI{4, Float64, 8}((10.0, 7.0, 3.0, 2.5, 4.0, 5.0, 0.7, 0.3), n)
    @test all(get_data(first(proj_tens(Val(:TI6), A8, n))) .≈ get_data(A8)[1:6])
    @test first(proj_tens(Val(:TI6), Array(get_array(A6)), n)) ≈ A6
    # TI2 ⋅ vector, against the matrix-vector product, for both forms.
    v = randn(rng, 3)
    for K in (TensTI{2}(2.0, 5.0, n), TensND.TensTI{2, Float64, 3}((2.0, 5.0, 0.7), n))
        @test K ⋅ v ≈ Array(get_array(K)) * v
        @test K ⋅ Vec{3}(Tuple(v)) ≈ Array(get_array(K)) * v
        @test (K ⋅ v) isa Vec{3}
    end
    @test_throws DimensionMismatch TensTI{2}(2.0, 5.0, n) ⋅ [1.0, 2.0]
    @test ForwardDiff.derivative(a -> (TensTI{2}(a, 5.0, n) ⋅ v)[1], 2.0) ≈ (v - (transpose(collect(n)) * v) * collect(n))[1]
    @syms a::real b::real
    Ks = TensTI{2}(a, b, (0, 0, 1))
    @test collect(Ks ⋅ Sym[1, 2, 3]) == [a, 2a, 3b]
    # The order-2 read-off in the frame of the axis.
    M = Array(get_array(TensND.TensTI{2, Float64, 3}((2.0, 5.0, 0.7), (0.0, 0.0, 1.0))))
    @test all(ti2_params_from_matrix(M) .≈ (2.0, 5.0, 0.7))
    @test all(ti2_params_from_matrix(M + M') .≈ (4.0, 10.0, 0.0))
end

@testsection "Projections about an axis or a frame: derivatives at e₃, symbols" begin
    # The axis or frame enters by its matrix, not by Euler angles: the
    # derivative with respect to the axis is finite at e₃ (it was NaN), and a
    # symbol not declared real is not conjugated. Reference: centered
    # differences with a step at the cube root of eps.
    A4 = Array(get_array(Tens(SymmetricTensor{4, 3}((i, j, k, l) -> 1.0 / (i + j + k + l)) + one(SymmetricTensor{4, 3}))))
    A2 = [3.0 0.4 0.1; 0.4 2.0 0.3; 0.1 0.3 1.0]
    h = cbrt(eps())
    cfd(f, x) = (f(x + h) - f(x - h)) / (2h)
    fti4(θ) = collect(get_data(first(proj_tens(Val(:TI), A4, (sin(θ), 0.0, cos(θ))))))
    fti2(θ) = collect(get_data(first(proj_tens(Val(:TI), A2, (sin(θ), 0.0, cos(θ))))))
    fortho(θ) = collect(get_data(first(proj_tens(Val(:ORTHO), A4, RotatedBasis(θ, 0.0, 0.0)))))
    fcub(θ) = collect(get_data(first(proj_tens(Val(:CUBIC), A4, RotatedBasis(θ, 0.0, 0.0)))))
    for f in (fti4, fti2, fortho, fcub), θ in (0.0, 0.4)
        d = ForwardDiff.derivative(f, θ)
        @test all(isfinite, d)
        @test d ≈ cfd(f, θ) rtol = 1.0e-6 atol = 1.0e-8
    end
    @syms φ
    B4 = first(proj_tens(Val(:TI), A4, (sin(φ), 0, cos(φ))))
    @test !occursin("conjugate", string(get_data(B4)))
    @test Float64.(subs.(collect(get_data(B4)), φ => 0.4)) ≈ fti4(0.4)
end
