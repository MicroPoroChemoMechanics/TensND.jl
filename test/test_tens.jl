@testsection "Tens" begin
    sq2 = √Sym(2)
    v = Sym[0 1 1; 1 0 1; 1 1 0]
    b = Basis(v)
    bn = normalize(b)
    for i in 1:3
        @eval $(Symbol("v$i")) = symbols($"v$i", real = true)
    end
    V = Tens(Tensor{1, 3}(i -> eval(Symbol("v$i"))))
    @test tsimplify(components(V, b, (:cont,))) ==
        [(-v1 + v2 + v3) / 2, (v1 - v2 + v3) / 2, (v1 + v2 - v3) / 2]
    @test tsimplify(components(V, b, (:cov,))) == [v2 + v3, v1 + v3, v1 + v2]
    @test tsimplify(components(V, bn, (:cov,))) ==
        [sq2 * (v2 + v3) / 2, sq2 * (v1 + v3) / 2, sq2 * (v1 + v2) / 2]

    for i in 1:3, j in 1:3
        @eval $(Symbol("t$i$j")) = symbols($"t$i$j", real = true)
    end
    T = Tens(Tensor{2, 3}((i, j) -> eval(Symbol("t$i$j"))))
    @test tsimplify(components(T, b, (:cov, :cov))) == [
        t22 + t23 + t32 + t33 t21 + t23 + t31 + t33 t21 + t22 + t31 + t32
        t12 + t13 + t32 + t33 t11 + t13 + t31 + t33 t11 + t12 + t31 + t32
        t12 + t13 + t22 + t23 t11 + t13 + t21 + t23 t11 + t12 + t21 + t22
    ]
    @test tsimplify(components(T, b, (:cont, :cov))) ==
        [
        -t12 - t13 + t22 + t23 + t32 + t33 -t11 - t13 + t21 + t23 + t31 + t33 -t11 - t12 + t21 + t22 + t31 + t32
        t12 + t13 - t22 - t23 + t32 + t33 t11 + t13 - t21 - t23 + t31 + t33 t11 + t12 - t21 - t22 + t31 + t32
        t12 + t13 + t22 + t23 - t32 - t33 t11 + t13 + t21 + t23 - t31 - t33 t11 + t12 + t21 + t22 - t31 - t32
    ] / 2

    for i in 1:3
        @eval $(Symbol("a$i")) = symbols($"a$i", real = true)
        @eval $(Symbol("b$i")) = symbols($"b$i", real = true)
    end
    a = Tens(Vec{3}((i) -> eval(Symbol("a$i"))))
    b = Tens(Vec{3}((i) -> eval(Symbol("b$i"))))
    @test a ⊗ b == Sym[a1 * b1 a1 * b2 a1 * b3; a2 * b1 a2 * b2 a2 * b3; a3 * b1 a3 * b2 a3 * b3]
    @test a ⊗ˢ b == Sym[
        a1 * b1 a1 * b2 / 2 + a2 * b1 / 2 a1 * b3 / 2 + a3 * b1 / 2
        a1 * b2 / 2 + a2 * b1 / 2 a2 * b2 a2 * b3 / 2 + a3 * b2 / 2
        a1 * b3 / 2 + a3 * b1 / 2 a2 * b3 / 2 + a3 * b2 / 2 a3 * b3
    ]

    θ, ϕ, ψ = symbols("θ, ϕ, ψ", real = true)
    R = rot3(θ, ϕ, ψ)
    Λ = [symbols("λ$i", positive = true) for i in 1:3]
    bg = Basis(R .* Λ')
    bo = Basis(R)
    𝕀, 𝕁, 𝕂 = ISO()
    α, β = symbols("α β", real = true)
    𝕋 = α * 𝕁 + β * 𝕂
    @test change_tens(𝕋, bg) == 𝕋
    @test change_tens(𝕋, bo) == 𝕋
    A = components(𝕋, bg)
    𝕋₂ = Tens(A, bg)
    @test 𝕋₂ == 𝕋


end


@testsection "pprint" begin
    # `pprint` is the display entry point, renamed twice: `intrinsic` (through
    # v0.2.7) → `print_tensor` (v0.3.0) → `pprint` (v0.3.1). Both former names
    # must keep forwarding. `_pprint_string` is the formatting half, split out
    # so the layout can be checked without capturing `stdout`.
    S = coorsys_spherical()
    𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(S)
    θ, ϕ, r = getcoords(S)

    # Expanded on the basis; a unit coefficient carries no parentheses, and
    # zero components are dropped entirely.
    t = 𝐞ʳ ⊗ 𝐞ʳ + 2 * 𝐞ᶿ ⊗ 𝐞ᶿ
    @test TensND._pprint_string(t; coords = ("θ", "ϕ", "r")) == "(2)𝐞ᶿ⊗𝐞ᶿ + 𝐞ʳ⊗𝐞ʳ"

    # A vanishing tensor prints as `0`, not as an empty line.
    @test TensND._pprint_string(t - t) == "0"

    # `vec` renames the basis symbol, `coords` its indices. The index sits
    # *above* because the components of `𝐞ʳ` are contravariant.
    @test TensND._pprint_string(𝐞ʳ; vec = '𝐚', coords = (1, 2, 3)) == "𝐚³"

    # The chart supplies the coordinate names, so the derived field reads in
    # spherical notation rather than 𝐞₁, 𝐞₂, 𝐞₃.
    g = change_tens(GRAD(𝐞ʳ, S), normalized_basis(S))
    @test TensND._pprint_string(g; coords = ("θ", "ϕ", "r")) == "(1/r)𝐞ᶿ⊗𝐞ᶿ + (1/r)𝐞ᵠ⊗𝐞ᵠ"

    # Every entry point returns `nothing` and writes to `stdout`; the two
    # deprecated aliases resolve (a malformed `@deprecate` would be a
    # `MethodError` here) — they emit a deprecation warning on stderr.
    @test pprint(t, S) === nothing
    @test print_tensor(t, S) === nothing
    @test intrinsic(t, S) === nothing

    # Fallback on a non-tensor: a scalar field, via the rich `text/plain`
    # display rather than `println`.
    @test pprint(LAPLACE(1 / r, S)) === nothing
end

@testsection "Storage of a rotated tensor — ForwardDiff.Dual" begin
    # Writing a structured tensor about a NON-canonical axis leaves minor-
    # antisymmetric round-off of order 1e-16. `_store_symmetric` absorbs exactly
    # that residue — but its tolerant method used to be written
    # `T <: AbstractFloat`, and `ForwardDiff.Dual` is not one, so a `Dual`-valued
    # tensor fell through to the exact comparison and landed on the
    # 81-component `Tensor` where the same tensor in `Float64` landed on the
    # 36-component `SymmetricTensor`.
    #
    # The two are NOT interchangeable downstream: `KM` then returns a 9×9 matrix
    # instead of 6×6, and `inv` of the full form solves a 9×9 system whose three
    # minor-antisymmetric directions are null. Automatic differentiation through
    # any rotated structured tensor was therefore impossible.
    ℬʳ = RotatedBasis(0.3, 0.4, 0.1)
    d = ForwardDiff.Dual{:tag}(0.8, 1.0)

    for t in (
            TensISO{3}(3 * 2.0, 2 * d), TensTI{4}(
                20.0 * d, 25.0, 7.0, 6.0, 5.0,
                (1, 1, 1) ./ √3
            ),
        )
        @test size(KM(t, ℬʳ)) == (6, 6)
        @test size(KM(t)) == size(KM(t, ℬʳ))          # the frame changes nothing
    end
    @test size(KM(TensISO{3}(d), ℬʳ)) == (6,)          # order 2

    # `Float64` and `Complex` must be unaffected, and both stay on the tolerant
    # side (`eps` is defined on a `Dual` type but not on a `Complex` one, which
    # is why the tolerance goes through `TensND._approx_eps`).
    @test size(KM(TensISO{3}(6.0, 1.6), ℬʳ)) == (6, 6)
    z = ComplexF64(1.0, 0.5)
    @test size(KM(TensISO{3}(3z, 2z), ℬʳ)) == (6, 6)
    @test TensND._approx_eps(ComplexF64) === eps(Float64)
    @test TensND._approx_eps(typeof(d)) === eps(Float64)

    # The point of it all: a derivative that reaches through the rotation.
    f(x) = Matrix(KM(TensISO{3}(6.0, 2x), ℬʳ))[4, 4]
    h = 1.0e-7
    @test ForwardDiff.derivative(f, 0.8) ≈ (f(0.8 + h) - f(0.8 - h)) / (2h) rtol = 1.0e-5

    # …and the tolerance still refuses a genuine minor antisymmetry rather than
    # averaging it away.
    raw = zeros(3, 3, 3, 3)
    raw[1, 2, 1, 1] = 1.0
    raw[2, 1, 1, 1] = -1.0
    @test !(TensND.tensor_or_array(raw) isa SymmetricTensor)
    rawd = [ForwardDiff.Dual{:tag}(v, 0.0) for v in raw]
    @test !(TensND.tensor_or_array(rawd) isa SymmetricTensor)
end

@testsection "A rotation by a dual angle keeps its derivatives at zero" begin
    # A rotation by a dual angle θ = 0 is the identity in value, but its
    # derivative is the generator of the rotation. `isidentity` compared values
    # only, so `RotatedBasis` collapsed to `CanonicalBasis` and every derivative
    # with respect to the angle came out zero — for any |θ| below about 1e-8,
    # the tolerance of `≈`.
    @test RotatedBasis(ForwardDiff.Dual(0.0, 1.0), 0.0, 0.0) isa RotatedBasis
    @test RotatedBasis(ForwardDiff.Dual(0.0, 1.0)) isa RotatedBasis
    @test RotatedBasis(0.0, 0.0, 0.0) isa CanonicalBasis           # unchanged
    @test RotatedBasis(ForwardDiff.Dual(0.0, 0.0), 0.0, 0.0) isa CanonicalBasis

    # References in closed form: with ϕ = ψ = 0, the ZYZ rotation is a rotation
    # by θ about e₂, and a vector of components v in it has the canonical
    # components R(θ)·v.
    v = [1.0, 2.0, 3.0]
    f(θ) = components_canon(Tens(v .* one(θ), RotatedBasis(θ, zero(θ), zero(θ))))
    dR(θ) = [-sin(θ) 0 cos(θ); 0 0 0; -cos(θ) 0 -sin(θ)]
    d2R(θ) = [-cos(θ) 0 -sin(θ); 0 0 0; sin(θ) 0 -cos(θ)]
    for θ in (0.0, 0.7)
        @test ForwardDiff.derivative(f, θ) ≈ dR(θ) * v
        @test ForwardDiff.derivative(θ -> ForwardDiff.derivative(f, θ), θ) ≈ d2R(θ) * v
    end
    g(θ) = components_canon(Tens([1.0, 2.0] .* one(θ), RotatedBasis(θ)))
    @test ForwardDiff.derivative(g, 0.0) ≈ [0.0 -1.0; 1.0 0.0] * [1.0, 2.0]

    # A structured tensor in that frame: dℂ/dθ at 0 is the generator W acting
    # on each index, Σ W_ia ℂ_ajkl + W_jb ℂ_ibkl + W_kc ℂ_ijcl + W_ld ℂ_ijkd.
    C(θ) = TensOrtho(30.0, 25.0, 20.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, RotatedBasis(θ, zero(θ), zero(θ)))
    c = components_canon(C(0.0))
    W = dR(0.0)
    dc = [
        sum(W[i, a] * c[a, j, k, l] + W[j, a] * c[i, a, k, l] + W[k, a] * c[i, j, a, l] + W[l, a] * c[i, j, k, a] for a in 1:3)
            for i in 1:3, j in 1:3, k in 1:3, l in 1:3
    ]
    @test ForwardDiff.derivative(θ -> Matrix(KM(C(θ))), 0.0) ≈ Matrix(KM(Tens(dc))) atol = 1.0e-12
    @test maximum(abs, Matrix(KM(Tens(dc)))) > 1                    # not vacuous

    # The diagonal test, for the same reason, and on exact element types, where
    # its absolute `eps(T)` tolerance threw: an integer or rational basis could
    # not be built.
    z = ForwardDiff.Dual(0.0, 1.0)
    @test !TensND.isdiagonal([1.0 z; 0.0 2.0])
    @test !isidentity([1.0 z; 0.0 1.0])
    @test TensND.isdiagonal([1 0; 0 2]) && !TensND.isdiagonal([1 1; 0 2])
    @test TensND.isdiagonal([1 // 2 0; 0 2]) && !TensND.isdiagonal([1 // 2 1; 0 2])
    @test Basis([2 0 0; 0 1 0; 0 0 1]) isa OrthogonalBasis
    @test Basis([2 // 1 0 0; 0 1 0; 0 0 1]) isa OrthogonalBasis
end

@testsection "Plain-number tensors in a basis of duals or symbols" begin
    # The components of a tensor written in another basis were converted back
    # to the element type of the tensor. For a `Float64` tensor in a basis of
    # duals or symbols that conversion throws (`Float64(::Dual)`,
    # `Float64(::Sym)`): they now take the type that holds both. The canonical
    # basis imposes no type of its own — its entries are exact zeros and ones.
    v = [1.0, 2.0, 3.0]
    R(θ) = [cos(θ) 0 sin(θ); 0 1 0; -sin(θ) 0 cos(θ)]
    dR(θ) = [-sin(θ) 0 cos(θ); 0 0 0; -cos(θ) 0 -sin(θ)]
    f(θ) = components_canon(Tens(v, RotatedBasis(θ, zero(θ), zero(θ))))
    for θ in (0.0, 0.7)
        @test ForwardDiff.derivative(f, θ) ≈ dR(θ) * v
    end

    # Into a rotated basis, for an order-2 tensor: components Rᵀ·A·R.
    A = [2.0 0.5 0.1; 0.5 3.0 0.2; 0.1 0.2 4.0]
    h(θ) = Matrix(components(Tens(A), RotatedBasis(θ, zero(θ), zero(θ)), (:cont, :cont)))
    θ = 0.7
    @test ForwardDiff.derivative(h, θ) ≈ dR(θ)' * A * R(θ) + R(θ)' * A * dR(θ)

    # A contraction of a plain tensor with a dual one, through the generic
    # array route (a fourth-order operand is no Tensors.jl `dot` pair).
    𝔸 = Tens([0.1 * i + 0.2 * j + 0.3 * k + 0.05 * l for i in 1:3, j in 1:3, k in 1:3, l in 1:3])
    n(x) = [sin(x), 0.0, cos(x)]
    a(x) = Array(components_canon(𝔸 ⋅ Tens(n(x))))
    dn = [cos(0.4), 0.0, -sin(0.4)]
    @test ForwardDiff.derivative(a, 0.4) ≈ [sum(get_array(𝔸)[i, j, k, l] * dn[l] for l in 1:3) for i in 1:3, j in 1:3, k in 1:3]

    # An orthogonal basis whose scale factors are duals: x = Σ vⁱ χᵢ 𝐫ᵢ, so
    # ∂x/∂χ₁ = v¹ 𝐫₁.
    ℬʳ = RotatedBasis(0.3, 0.2, 0.1)
    s(χ) = components_canon(Tens(v, Basis(ℬʳ, [χ, 1.0, 1.0])))
    @test Basis(ℬʳ, [ForwardDiff.Dual(2.0, 1.0), 1.0, 1.0]) isa OrthogonalBasis
    @test ForwardDiff.derivative(s, 2.0) ≈ v[1] * Array(vecbasis(ℬʳ, :cov))[:, 1]

    # Symbolic: a plain vector in a basis rotated by a symbolic angle.
    @syms ϑ::real
    c = components_canon(Tens([1.0, 0.0], RotatedBasis(ϑ)))
    @test eltype(c) <: Sym
    @test simplify(c[1] - cos(ϑ)) == 0 && simplify(c[2] - sin(ϑ)) == 0

    # Plain tensors in plain bases keep their element type.
    @test eltype(components(Tens(A), RotatedBasis(0.3, 0.2, 0.1), (:cont, :cont))) == Float64
    @test eltype(components(Tens(A, RotatedBasis(0.3, 0.2, 0.1)), CanonicalBasis(), (:cont, :cont))) == Float64
end

@testsection "symbolic passes reach structured tensors and tuples" begin
    # REGRESSION. `structured_tens_ops.jl` implements the whole `tsimplify`
    # family as `_rebuild(A, OP(get_data(A)))`, and `get_data` returns an
    # `NTuple`. A `Tuple` is not an `AbstractArray`, so the elementwise methods
    # never matched and the generic identity fallback caught the call: every
    # one of these was a SILENT no-op on `TensISO`, `TensTI` and `TensOrtho`.
    @syms a::positive b::positive w::positive

    t = TensISO{3}(sin(a)^2 + cos(a)^2, 2b)
    @test isequal(get_data(tsimplify(t))[1], Sym(1))

    ti = TensTI{4}(sin(a)^2 + cos(a)^2, b, a, 2a, 3a, (0, 0, 1))
    @test isequal(get_data(tsimplify(ti))[1], Sym(1))
    @test isequal(get_data(tsubs(ti, b => 5))[2], Sym(5))

    # Tuples directly.
    @test isequal(tsimplify((sin(a)^2 + cos(a)^2, 2b))[1], Sym(1))
end

@testsection "tlimit — the limit passage keeps the symmetry class" begin
    @syms a::positive b::positive w::positive

    # Scalars, both sided and one-sided.
    @test isequal(tlimit(1 / (1 + w), w, 0), Sym(1))
    @test isequal(tlimit(1 / w, w, 0, "+"), oo)
    @test isequal(tlimit(a / (a + 1 / w), w, 0), Sym(0))

    # Identity on numeric types — a number has no free symbol.
    @test tlimit(3.0, :anything, 0) === 3.0
    @test tlimit([1.0, 2.0], :anything, 0) == [1.0, 2.0]

    # Structured tensors: the class survives, and only the coefficients move.
    ti = TensTI{4}(1 / (1 + w), b * w, a, 2a, 3a, (0, 0, 1))
    lim = tlimit(ti, w, 0)
    @test lim isa TensTI{4}
    @test axis(lim) == axis(ti)
    @test all(isequal.(collect(get_data(lim)), Sym[1, 0, a, 2a, 3a]))

    iso = TensISO{3}(3 / (1 + w), 2b * w)
    @test all(isequal.(collect(get_data(tlimit(iso, w, 0))), Sym[3, 0]))

    # An incompressible limit, the motivating case.
    @syms k::positive μ::positive
    @test isequal(tlimit(TensND.get_data(TensISO{3}(3k, 2μ))[2], k, oo), 2μ)

    # Symbolics has no limit and must say so rather than return its input.
    # `Symbolics` is a hard dependency of TensND, reached here through it since
    # the test environment does not import it directly.
    TensND.Symbolics.@variables z
    @test_throws ArgumentError tlimit(z, z, 0)
end
