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
