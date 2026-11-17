# Lazy tensors: the zero tensor and the multiples of the identity. Every
# operation is checked against the same operation on the full arrays of
# canonical components, formed independently of the lazy methods.

import Symbolics

# Canonical contravariant components, as a plain array.
_lz_full(t) = Array(components(t, CanonicalBasis{get_dim(t), Float64}(), ntuple(_ -> :cont, Val(get_order(t)))))
_lz_full(t::Number) = t
_lz_δ(i, j) = i == j ? 1.0 : 0.0
_lz_id2(d) = [_lz_δ(i, j) for i in 1:d, j in 1:d]
_lz_id4(d) = [(_lz_δ(i, k) * _lz_δ(j, l) + _lz_δ(i, l) * _lz_δ(j, k)) / 2 for i in 1:d, j in 1:d, k in 1:d, l in 1:d]
_lz_id4f(d) = [_lz_δ(i, k) * _lz_δ(j, l) for i in 1:d, j in 1:d, k in 1:d, l in 1:d]

@testsection "Lazy tensors" begin
    rng = Random.MersenneTwister(4)
    𝟎₂, 𝟎₄ = TensZero{2, 3, Float64}(), TensZero{4, 3, Float64}()
    𝟏, 𝕀 = TensIdentity{2, 3}(2.0), TensIdentity{4, 3}(3.0)
    A = Tens(randn(rng, 3, 3))                      # no symmetry
    v = Tens(randn(rng, 3))
    C = Tens(randn(rng, 3, 3, 3, 3))
    S = Tens(SymmetricTensor{2, 3}((i, j) -> i + j + 0.5))

    @testsection "Interface" begin
        @test_throws ArgumentError TensIdentity{3, 3}(1.0)
        @test TensIdentity{2, 3}() isa TensIdentity{2, 3, Bool}
        for d in (2, 3)
            @test Array(TensIdentity{2, d}(2.0)) == 2 * _lz_id2(d)
            @test Array(TensIdentity{4, d}(2.0)) ≈ 2 * _lz_id4(d)
            @test Array(TensZero{3, d, Float64}()) == zeros(d, d, d)
            @test mandel(TensIdentity{4, d}(2.0)) ≈ mandel(TensISO(TensIdentity{4, d}(2.0)))
            @test mandel(TensIdentity{2, d}(2.0)) ≈ mandel(TensISO(TensIdentity{2, d}(2.0)))
            @test mandel(TensZero{4, d, Float64}()) == zeros(d * (d + 1) ÷ 2, d * (d + 1) ÷ 2)
            @test mandel(TensZero{2, d, Float64}()) == zeros(d * (d + 1) ÷ 2)
        end
        @test get_data(𝟎₂) == () && get_data(𝕀) == (3.0,)
        @test get_basis(𝕀) == CanonicalBasis{3, Float64}()
        @test get_var(𝕀) == (:cont, :cont, :cont, :cont) && get_var(𝟏, 2) == :cont
        @test components(𝟏) == 2 * _lz_id2(3)
        @test components(𝟏, (:cov, :cont)) == 2 * _lz_id2(3)
        ℬ = random_basis(rng, Float64)
        @test components(𝕀, ℬ, (:cont, :cov, :cont, :cov)) ≈ 3 * _lz_id4(3)
        @test _lz_full(express(𝕀, ℬ)) ≈ 3 * _lz_id4(3)
        @test _lz_full(express(𝟏, ℬ, (:cont, :cov))) ≈ 2 * _lz_id2(3)
        @test _lz_full(express(TensZero{3, 3, Float64}(), ℬ)) == zeros(3, 3, 3)
        @test TensISO(𝟏) == TensISO{3}(2.0) && TensISO(𝕀) == TensISO{3}(3.0, 3.0)
        @test TensISO(𝟎₂) == TensISO{3}(0.0) && TensISO(𝟎₄) == TensISO{3}(0.0, 0.0)
        @test isisotropic(𝕀) && !istransverselyisotropic(𝕀) && !isorthotropic(𝟎₂)
        @test issymmetric(𝟏) && issymmetric(𝟎₄) && isminorsymmetric(𝕀) && ismajorsymmetric(𝟎₄)
        @test iszero(𝟎₂) && !iszero(𝕀) && iszero(TensIdentity{2, 3}(0.0))
        @test sprint(show, 𝟎₂) == "𝟎" && sprint(show, 𝟏) == "(2.0) 𝟏" && sprint(show, 𝕀) == "(3.0) 𝕀"
        @test occursin("lazy: (3.0) 𝕀", sprint(show, MIME"text/plain"(), 𝕀))
    end

    @testsection "Scalars, inverse, powers, norms" begin
        @test -𝟎₂ === 𝟎₂ && -𝟏 == TensIdentity{2, 3}(-2.0)
        @test 2 * 𝟎₂ === 𝟎₂ && 𝟎₂ * 2 === 𝟎₂
        @test eltype(𝟎₂ / 2) === Float64 && eltype(TensZero{2, 3, Int}() / 2) === Float64
        @test 2 * 𝕀 == TensIdentity{4, 3}(6.0) && 𝕀 * 2 == 2 * 𝕀 && 𝕀 / 2 == TensIdentity{4, 3}(1.5)
        Symbolics.@variables a
        @test isequal((𝕀 / a).λ, 3.0 / a)
        @test eltype(𝟎₂ / a) === Symbolics.Num
        @test zero(𝕀) === 𝟎₄ && one(𝟎₂) == TensIdentity{2, 3}(1.0)
        @test transpose(𝕀) === 𝕀 && adjoint(𝟎₂) === 𝟎₂
        @test adjoint(TensIdentity{2, 3}(1.0 + 2im)) == TensIdentity{2, 3}(1.0 - 2im)
        @test inv(𝕀) == TensIdentity{4, 3}(1 / 3)
        @test_throws DomainError inv(𝟎₂)
        @test 𝕀^2 == TensIdentity{4, 3}(9.0) && 𝟏^-1 == TensIdentity{2, 3}(0.5) && 𝕀^-1 == inv(𝕀)
        @test 𝟎₂^2 === 𝟎₂ && 𝟎₄^0 == TensIdentity{4, 3}(1.0)
        @test_throws DomainError 𝟎₂^-1
        @test_throws DomainError 𝟎₄^-2
        @test tr(𝟎₂) == 0 && tr(𝟏) == 6.0
        @test norm(𝟎₄) == 0 && norm(𝟏) ≈ norm(_lz_full(𝟏)) && norm(𝕀) ≈ norm(_lz_full(𝕀))
        @test 𝟎₂ == TensZero{2, 3, Int}() && 𝕀 == TensIdentity{4, 3}(3) && 𝕀 != TensIdentity{4, 3}(2.0)
        @test +𝕀 === 𝕀
    end

    @testsection "Sums" begin
        @test 𝟎₂ + A === A && A + 𝟎₂ === A && A - 𝟎₂ === A
        @test _lz_full(𝟎₂ - A) == -_lz_full(A)
        @test 𝟎₂ + TensZero{2, 3, Int}() isa TensZero{2, 3, Float64} && 𝟎₂ - 𝟎₂ === 𝟎₂
        @test 𝟎₂ + 𝟏 == 𝟏 && 𝟎₂ - 𝟏 == -𝟏 && 𝟏 + 𝟎₂ == 𝟏 && 𝟏 - 𝟎₂ == 𝟏
        @test 𝕀 + 𝕀 == TensIdentity{4, 3}(6.0) && 𝕀 - 𝕀 == TensIdentity{4, 3}(0.0)
        for (L, B) in ((𝟏, A), (𝕀, C), (𝟏, TensISO{3}(5.0)), (𝕀, TensISO{3}(5.0, 2.0)), (𝟎₄, TensISO{3}(5.0, 2.0)))
            @test _lz_full(L + B) ≈ _lz_full(L) + _lz_full(B)
            @test _lz_full(B + L) ≈ _lz_full(L) + _lz_full(B)
            @test _lz_full(L - B) ≈ _lz_full(L) - _lz_full(B)
            @test _lz_full(B - L) ≈ _lz_full(B) - _lz_full(L)
        end
    end

    @testsection "Products and contractions of the zero tensor" begin
        partners = (
            A, C, v, TensISO{3}(2.0), TensISO{3}(2.0, 1.0), rand(rng, TensTI{2}), rand(rng, TensTI{4}),
            rand(rng, TensOrtho{2}), rand(rng, TensOrtho{4}), rand(rng, TensCubic),
        )
        for (op, n) in ((dot, 1), (dcontract, 2), (qcontract, 4), (otimes, 0), (Tensors.otimesu, 0), (TensND.otimesul, 0), (TensND.sotimes, 0))
            for B in partners, Z in (TensZero{get_order(B), 3, Float64}(), TensZero{4, 3, Float64}())
                o = get_order(Z) + get_order(B) - 2n
                if get_order(Z) ≥ n && get_order(B) ≥ n
                    @test op(Z, B) == (o == 0 ? 0.0 : TensZero{o, 3, Float64}())
                    @test op(B, Z) == (o == 0 ? 0.0 : TensZero{o, 3, Float64}())
                else
                    @test_throws DimensionMismatch op(Z, B)
                end
            end
            @test op(𝟎₄, 𝟎₄) == (8 - 2n == 0 ? 0.0 : TensZero{8 - 2n, 3, Float64}())
            # with no array formed
            f(Z, B) = op(Z, B)
            f(𝟎₄, C)
            @test (@allocated f(𝟎₄, C)) == 0
        end
    end

    @testsection "Contractions of the identities" begin
        λ₂, λ₄ = 𝟏.λ, 𝕀.λ
        # λ 𝟏 ⋅ 𝐁 = λ 𝐁, on either side, whatever the order of 𝐁
        for B in (v, A, C, TensISO{3}(2.0), rand(rng, TensTI{2}), rand(rng, TensOrtho{2}), TensISO{3}(2.0, 1.0), rand(rng, TensTI{4}), rand(rng, TensOrtho{4}), rand(rng, TensCubic))
            @test _lz_full(dot(𝟏, B)) ≈ λ₂ * _lz_full(B)
            @test _lz_full(dot(B, 𝟏)) ≈ λ₂ * _lz_full(B)
        end
        @test dot(𝟏, 𝟏) == TensIdentity{2, 3}(4.0)
        @test dot(𝟏, 𝟎₄) === 𝟎₄ && dot(𝟎₂, 𝟏) === 𝟎₂
        # λ 𝟏 ⊡ 𝐁 contracts the first two indices of 𝐁, 𝐁 ⊡ λ 𝟏 the last two
        c = _lz_full(C)
        @test dcontract(𝟏, A) ≈ λ₂ * tr(_lz_full(A))
        @test _lz_full(dcontract(𝟏, C)) ≈ λ₂ * [sum(c[i, i, k, l] for i in 1:3) for k in 1:3, l in 1:3]
        @test _lz_full(dcontract(C, 𝟏)) ≈ λ₂ * [sum(c[i, j, k, k] for k in 1:3) for i in 1:3, j in 1:3]
        # λ 𝕀 ⊡ 𝐁: the symmetric part in the contracted pair
        a = _lz_full(A)
        @test _lz_full(dcontract(𝕀, A)) ≈ λ₄ * (a + a') / 2
        @test _lz_full(dcontract(A, 𝕀)) ≈ λ₄ * (a + a') / 2
        @test _lz_full(dcontract(𝕀, C)) ≈ λ₄ * (c + permutedims(c, (2, 1, 3, 4))) / 2
        @test _lz_full(dcontract(C, 𝕀)) ≈ λ₄ * (c + permutedims(c, (1, 2, 4, 3))) / 2
        # a tensor symmetric by its type comes back multiplied, not formed
        for B in (S, TensISO{3}(2.0), TensISO{3}(2.0, 1.0), rand(rng, TensTI{2}), rand(rng, TensTI{4}), rand(rng, TensOrtho{2}), rand(rng, TensOrtho{4}), rand(rng, TensCubic))
            @test dcontract(𝕀, B) ≈ λ₄ * B
            @test dcontract(B, 𝕀) ≈ λ₄ * B
            @test _lz_full(dcontract(𝟏, B)) ≈ λ₂ * (get_order(B) == 2 ? tr(_lz_full(B)) : _lz_full(dcontract(TensISO{3}(1.0), B)))
            @test _lz_full(dcontract(B, 𝟏)) ≈ λ₂ * (get_order(B) == 2 ? tr(_lz_full(B)) : _lz_full(dcontract(B, TensISO{3}(1.0))))
        end
        g(L, B) = dcontract(L, B)
        Ciso = TensISO{3}(2.0, 1.0)
        g(𝕀, Ciso)
        @test (@allocated g(𝕀, Ciso)) == 0
        # vectors are refused by the double contraction
        @test_throws Exception dcontract(𝟏, v)
        @test_throws Exception dcontract(𝕀, v)
        # two identities
        @test dcontract(𝕀, 𝕀) == TensIdentity{4, 3}(9.0)
        @test dcontract(𝕀, 𝟏) == TensIdentity{2, 3}(6.0) && dcontract(𝟏, 𝕀) == TensIdentity{2, 3}(6.0)
        @test dcontract(𝟏, 𝟏) == 12.0
        @test dcontract(𝕀, 𝟎₂) === 𝟎₂ && dcontract(𝟎₄, 𝟏) === 𝟎₂
    end

    @testsection "Constructors from a tensor" begin
        Ciso = TensISO{3}(3.0, 4.0)
        @test TensZero(Ciso) === TensZero{4, 3, Float64}()
        @test TensZero(A) === TensZero{2, 3, Float64}()
        @test TensIdentity(Ciso) === TensIdentity{4, 3}(1.0)
        @test TensIdentity(TensISO{2}(1.0f0)) === TensIdentity{2, 2}(1.0f0)
        @test TensIdentity(Ciso, Val(:full)) === TensIdentity{4, 3}(1.0, Val(:full))
        @test_throws ArgumentError TensIdentity(v)
        # zero and one stay in the type of their argument
        @test zero(Ciso) isa TensISO && one(Ciso) isa TensISO
    end

    @testsection "Full identity" begin
        λ = 1.5
        𝕀ᶠ = TensIdentity{4, 3}(λ, Val(:full))
        @test 𝕀ᶠ isa TensIdentity{4, 3, Float64, :full}
        @test TensIdentity{2, 3}(1.0, Val(:full)) isa TensIdentity{2, 3, Float64, :sym}
        @test_throws ArgumentError TensIdentity{4, 3}(1.0, Val(:other))
        for d in (2, 3)
            I = TensIdentity{4, d}(λ, Val(:full))
            @test Array(I) == λ * _lz_id4f(d)
            @test mandel(I) == λ * Matrix(1.0LinearAlgebra.I, d^2, d^2)
            @test norm(I) ≈ norm(λ * _lz_id4f(d))
        end
        @test !isisotropic(𝕀ᶠ) && !issymmetric(𝕀ᶠ) && !Tensors.isminorsymmetric(𝕀ᶠ) && Tensors.ismajorsymmetric(𝕀ᶠ)
        @test sprint(show, 𝕀ᶠ) == "(1.5) 𝕀ᶠ"
        @test sprint(show, MIME("text/plain"), 𝕀ᶠ) == "TensIdentity{4, 3, Float64, :full}\n→ lazy: (1.5) 𝕀ᶠ"
        # the same in every orthonormal basis
        ℬ = random_basis(rng, Float64)
        @test _lz_full(express(𝕀ᶠ, ℬ)) ≈ λ * _lz_id4f(3)
        @test Array(components(𝕀ᶠ, ℬ, (:cont, :cont, :cont, :cont))) ≈ λ * _lz_id4f(3)
        # scalars keep the kind
        for t in (-𝕀ᶠ, 2𝕀ᶠ, 𝕀ᶠ * 2, 𝕀ᶠ / 2, inv(𝕀ᶠ), 𝕀ᶠ^2, one(𝕀ᶠ), adjoint(𝕀ᶠ), transpose(𝕀ᶠ))
            @test t isa TensIdentity{4, 3, Float64, :full}
        end
        @test get_data(inv(𝕀ᶠ)) == (1 / λ,) && get_data(𝕀ᶠ^2) == (λ^2,) && get_data(one(𝕀ᶠ)) == (1.0,)
        @test 𝕀ᶠ / Symbolics.Num(2) isa TensIdentity{4, 3, <:Any, :full}
        # equality: 𝕀 and 𝕀ᶠ differ unless both vanish
        @test 𝕀ᶠ == TensIdentity{4, 3}(λ, Val(:full))
        @test 𝕀ᶠ != TensIdentity{4, 3}(λ) && TensIdentity{4, 3}(λ) != 𝕀ᶠ
        @test 0𝕀ᶠ == TensIdentity{4, 3}(0.0) && TensIdentity{4, 3}(0.0) == 0𝕀ᶠ
        # sums: one kind adds the factors, otherwise the components
        @test 𝕀ᶠ + 𝕀ᶠ === TensIdentity{4, 3}(2λ, Val(:full))
        @test 𝕀ᶠ - 𝕀ᶠ === TensIdentity{4, 3}(0.0, Val(:full))
        @test _lz_full(𝕀ᶠ + 𝕀) ≈ λ * _lz_id4f(3) + 3 * _lz_id4(3)
        @test _lz_full(𝕀 - 𝕀ᶠ) ≈ 3 * _lz_id4(3) - λ * _lz_id4f(3)
        @test _lz_full(𝕀ᶠ + C) ≈ λ * _lz_id4f(3) + _lz_full(C)
        @test _lz_full(TensISO{3}(2.0, 1.0) - 𝕀ᶠ) ≈ _lz_full(TensISO{3}(2.0, 1.0)) - λ * _lz_id4f(3)
        @test 𝟎₄ + 𝕀ᶠ === 𝕀ᶠ && 𝕀ᶠ - 𝟎₄ === 𝕀ᶠ && 𝟎₄ - 𝕀ᶠ === -𝕀ᶠ
        # contractions: the identity map, on either side, whatever the symmetry
        @test _lz_full(dcontract(𝕀ᶠ, A)) ≈ λ * _lz_full(A)
        @test _lz_full(dcontract(A, 𝕀ᶠ)) ≈ λ * _lz_full(A)
        @test _lz_full(dcontract(𝕀ᶠ, C)) ≈ λ * _lz_full(C)
        @test _lz_full(dcontract(C, 𝕀ᶠ)) ≈ λ * _lz_full(C)
        for B in (TensISO{3}(2.0, 1.0), rand(rng, TensTI{4}), rand(rng, TensOrtho{4}), rand(rng, TensCubic))
            @test dcontract(𝕀ᶠ, B) ≈ λ * B && dcontract(B, 𝕀ᶠ) ≈ λ * B
        end
        # products of identities, against their components
        I4 = (𝕀, 𝕀ᶠ)
        for a in I4, b in I4
            p = dcontract(a, b)
            @test _lz_full(p) ≈ reshape(reshape(_lz_full(a), 9, 9) * reshape(_lz_full(b), 9, 9), 3, 3, 3, 3)
            @test p isa TensIdentity{4, 3, Float64, (a === b === 𝕀ᶠ ? :full : :sym)}
        end
        @test dcontract(𝕀ᶠ, 𝟏) === TensIdentity{2, 3}(2λ) && dcontract(𝟏, 𝕀ᶠ) === TensIdentity{2, 3}(2λ)
        @test dcontract(𝕀ᶠ, 𝟎₄) === 𝟎₄ && dcontract(𝟎₂, 𝕀ᶠ) === 𝟎₂
        Ciso = TensISO{3}(2.0, 1.0)
        g(L, B) = dcontract(L, B)
        g(𝕀ᶠ, Ciso)
        @test (@allocated g(𝕀ᶠ, Ciso)) == 0
        @test (@inferred dcontract(𝕀ᶠ, Ciso)) isa TensISO{4, 3, Float64}
        @test ForwardDiff.derivative(t -> _lz_full(dcontract(TensIdentity{4, 3}(t, Val(:full)), C)), 0.4) ≈ _lz_full(C)
    end

    @testsection "Derivatives and symbols" begin
        Ciso = TensISO{3}(2.0, 1.0)
        # d(λ 𝕀 ⊡ ℂ)/dλ = ℂ
        @test ForwardDiff.derivative(λ -> _lz_full(dcontract(TensIdentity{4, 3}(λ), Ciso)), 1.5) ≈ _lz_full(Ciso)
        @test ForwardDiff.derivative(λ -> _lz_full(TensIdentity{2, 3}(λ) + A), 0.7) ≈ _lz_full(TensIdentity{2, 3}(1.0))
        Symbolics.@variables μ
        Ls = TensIdentity{4, 3}(μ)
        @test isequal(get_data(dcontract(Ls, TensIdentity{4, 3}(2))), (2μ,))
        @test isequal(tr(TensIdentity{2, 3}(μ)), 3μ)
        @test (@inferred dcontract(𝕀, Ciso)) isa TensISO{4, 3, Float64}
        @test (@inferred dot(𝟏, v)) == 2 * v
    end
end
