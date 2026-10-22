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
    # `pprint` is the display entry point; `intrinsic` and `print_tensor`, its
    # former names, must keep forwarding. `_pprint_string` is the formatting
    # half, split out so the layout can be checked without capturing `stdout`.
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

    # A symbolic tensor built from `Sym[…]` (abstract element type `Sym`)
    # has components of the concrete `Sym{PyObject}`, as the conversions
    # through broadcasting always gave them.
    ℬs = Basis(Sym[1 0 0; 0 1 0; 0 1 1])
    ts = Tens(metric(ℬs, :cov), ℬs, (:cov, :cov))
    @test isconcretetype(eltype(components(ts, (:cont, :cov))))
    @test isconcretetype(eltype(components(ts, CanonicalBasis{3, Sym}(), (:cont, :cont))))

    # Plain tensors in plain bases keep their element type.
    @test eltype(components(Tens(A), RotatedBasis(0.3, 0.2, 0.1), (:cont, :cont))) == Float64
    @test eltype(components(Tens(A, RotatedBasis(0.3, 0.2, 0.1)), CanonicalBasis(), (:cont, :cont))) == Float64
end

@testsection "Equality, hashing, Kelvin-Mandel signature, display" begin
    # `==` compares tensors as mathematical objects, in whatever basis. The
    # same comparison stood for `isequal`, whose contract with `hash` it broke:
    # the same vector written in a permuted basis was `isequal` to itself in
    # the canonical one, with different components and so a different hash.
    P = Basis([0.0 1.0 0.0; 0.0 0.0 1.0; 1.0 0.0 0.0])
    t = Tens([1.0, 2.0, 3.0])
    tp = change_tens(t, P)
    @test t == tp
    tensors = (
        t, tp, Tens([1.0, 2.0, 3.0]), TensISO{3}(3.0, 2.0),
        Tens(get_array(TensISO{3}(3.0, 2.0))), Tens([1.0, 2.0, 3.0], (:cov,), P),
    )
    for a in tensors, b in tensors
        isequal(a, b) && @test hash(a) == hash(b)
    end
    @test isequal(t, Tens([1.0, 2.0, 3.0]))                          # structural
    @test isequal(TensISO{3}(3.0, 2.0), Tens(get_array(TensISO{3}(3.0, 2.0))))
    @test Dict(t => 1)[Tens([1.0, 2.0, 3.0])] == 1
    @test length(Set((t, tp, Tens([1.0, 2.0, 3.0])))) == 2

    # `KM(t, b, var)` took its keywords as trailing positional arguments, so
    # none could be passed.
    C = TensISO{3}(3.0, 2.0)
    ℬʳ = RotatedBasis(0.3, 0.2, 0.1)
    @test KM(C, ℬʳ, (:cont, :cont, :cont, :cont); order = nothing) == KM(C, ℬʳ, (:cont, :cont, :cont, :cont))

    # Displaying a Tensors.jl fourth-order tensor is left to Tensors.jl.
    @test parentmodule(which(show, Tuple{Tensors.SymmetricTensor{4, 3, Float64, 36}})) !== TensND

    # A basis is built from a square matrix, checked by an exception that
    # cannot be compiled out.
    @test_throws DimensionMismatch Basis([1.0 0.0 0.0; 0.0 1.0 0.0])
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

@testsection "Short cuts and less common shapes" begin
    # Same basis, same variances: the components, or the tensor, come back
    # as they are.
    ℬ = Basis([1.0 0.2 0.1; 0.0 1.0 0.3; 0.1 0.0 1.0])
    t = Tens([1.0 2.0 0.5; 0.3 4.0 1.0; 0.2 0.1 3.0], ℬ, (:cont, :cov))
    @test components(t, get_var(t)) === get_array(t)
    @test change_tens(t, get_var(t)) === t
    o = Tens([2.0 0.5 0.1; 0.5 3.0 0.2; 0.1 0.2 4.0])
    @test components(o, (:cov, :cont)) === get_array(o)
    @test components(o, get_basis(o)) === get_array(o)
    ℬo = Basis(RotatedBasis(0.3, 0.2, 0.1), [2.0, 1.0, 0.5])
    to = Tens([1.0, 2.0, 3.0], ℬo, (:cov,))
    @test components(to, get_basis(to), (:cont,)) == components(to, (:cont,))

    # The product of a tensor with itself takes a short cut; it agrees with
    # the general route, which the same tensor written in another basis takes.
    v = Tens([1.0, 2.0, 3.0], ℬ, (:cont,))
    v2 = change_tens(v, Basis([1.0 0.0 0.3; 0.2 1.0 0.0; 0.0 0.1 1.0]))
    @test !isequal(v, v2) && components_canon(v) ≈ components_canon(v2)
    @test components_canon(v ⊗ v) ≈ components_canon(v ⊗ v2)
    # A non-symmetric order-2 tensor, and a vector of dimension 4, have no
    # symmetric storage for their square: they used to raise a MethodError.
    # Beyond dimension 3, where Tensors.jl has no type, the components stay a
    # plain array; building such a tensor used to raise a MethodError too.
    @test get_array(Tens([1.0, 2.0, 3.0, 4.0])) isa Vector{Float64}
    @test get_array(Tens(reshape(collect(1.0:16.0), 4, 4))) isa Matrix{Float64}
    @test components_canon(t ⊗ t) ≈ components_canon(t ⊗ change_tens(t, RotatedBasis(0.3, 0.2, 0.1)))
    a = [1.0 2.0 0.5; 0.3 4.0 1.0; 0.2 0.1 3.0]
    @test get_array(Tens(a) ⊗ Tens(a)) == [a[i, j] * a[k, l] for i in 1:3, j in 1:3, k in 1:3, l in 1:3]
    w = Tens([1.0, 2.0, 3.0, 4.0])
    @test get_array(w ⊗ w) == [1.0, 2.0, 3.0, 4.0] * [1.0, 2.0, 3.0, 4.0]'
    u = Tens([1.0, 2.0, 3.0])
    @test get_array(u ⊗ˢ u) == get_array(u ⊗ u)
    # The Frobenius norm of an order-2 tensor, in an orthonormal basis.
    @test norm(o) ≈ norm(get_array(o))

    # Kelvin-Mandel matrices of every shape go back to the right tensor type;
    # a shape of none is refused.
    A = [0.1 * (i + 2j) * (k + l) + (i == k) for i in 1:3, j in 1:3, k in 1:3, l in 1:3]
    @test get_array(inv_KM(KM(Tens(A)))) ≈ A
    @test get_array(inv_KM(KM(Tens(A)))) isa Tensors.Tensor{4, 3}
    S2 = SymmetricTensor{4, 2}((i, j, k, l) -> 1.0 * (i + j) * (k + l))
    @test get_array(inv_KM(KM(S2))) ≈ S2
    T2 = Tensor{4, 2}((i, j, k, l) -> 1.0 * (i + 2j) * (k + l))
    @test get_array(inv_KM(KM(T2))) ≈ T2
    @test_throws DimensionMismatch inv_KM(rand(5, 5))
    @test_throws DimensionMismatch inv_KM(rand(5))
end

@testsection "zero of a vector, and tensors of dimension 1" begin
    # The zero of an order-2 or order-4 tensor is isotropic (`TensISO`); that
    # of a vector stays in its basis, with its variance.
    ℬ = Basis([1.0 0.5 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], :cov)
    z = zero(Tens([1.0, 2.0, 3.0], ℬ, (:cov,)))
    @test iszero(get_array(z))
    @test get_basis(z) === ℬ && get_var(z) == (:cov,)
    # Dimension 1 has Tensors.jl types too, and every such tensor is symmetric.
    @test get_array(Tens(fill(2.0, 1, 1))) == SymmetricTensor{2, 1}((2.0,))
    @test get_array(Tens(fill(2.0, 1, 1, 1, 1))) == SymmetricTensor{4, 1}((2.0,))
end

@testsection "Symmetry predicates on symbolic components" begin
    # A macro is resolved before the testset runs: name it by its module path.
    TensND.Symbolics.@variables a b
    S = TensND.Symbolics
    # TensND's own predicates (Tensors.jl's are not redefined on its types).
    @test TensND._issym(Tensor{2, 2, S.Num}((a, b, b, a)))
    @test !TensND._issym(Tensor{2, 2, S.Num}((a, b, a, a)))
    # Order 4: one component without its minor partner, then one without its
    # major partner.
    m = Tensor{4, 2, S.Num}((i, j, k, l) -> (i, j, k, l) == (1, 2, 1, 1) ? a : zero(a))
    @test !TensND._isminorsym(m)
    M = Tensor{4, 2, S.Num}((i, j, k, l) -> (i, j, k, l) == (1, 1, 2, 2) ? a : zero(a))
    @test TensND._isminorsym(M)
    @test !TensND._ismajorsym(M)
    # A symmetric symbolic array is stored symmetric, a tensor built from it too.
    @test get_array(Tens([a b; b a])) isa SymmetricTensor{2, 2}
    # Loading TensND leaves Tensors.jl's own predicates exact on floats.
    @test !issymmetric(Tensor{2, 2}((1.0, 2.0, 2.0 + 1.0e-15, 1.0)))
    S4 = SymmetricTensor{4, 2}((i, j, k, l) -> 10i + j + 0.1k + 0.01l)
    S4ᵀ = TensND._majortranspose(S4)
    @test S4ᵀ isa SymmetricTensor{4, 2}
    @test all(S4ᵀ[i, j, k, l] == S4[k, l, i, j] for i in 1:2, j in 1:2, k in 1:2, l in 1:2)
end

@testsection "Rich display goes to the io it is given" begin
    n = (0.0, 0.0, 1.0)
    for x in (
            TensISO{3}(3.0, 2.0), TensISO{3}(2.0), TensTI{4}(1.0, 2.0, 3.0, 4.0, 5.0, n),
            TensTI{2}(1.0, 2.0, n), TensOrtho(9.0, 8.0, 7.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, CanonicalBasis{3, Float64}()),
            TensCubic(9.0, 4.0, 3.0, CanonicalBasis{3, Float64}()),
        )
        rich = repr(MIME"text/plain"(), x)
        @test startswith(rich, string(typeof(x)))
        @test occursin("→ decomposition: " * sprint(show, x), rich)
        @test occursin("→ KM", rich) == (ndims(x) == 4)
        # `print` is the compact form, written where it is asked to.
        @test sprint(print, x) == sprint(show, x)
    end
    @test occursin("KM (material frame)", repr(MIME"text/plain"(), TensOrtho(9.0, 8.0, 7.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, CanonicalBasis{3, Float64}())))
    t = Tens([1.0 2.0; 3.0 4.0], Basis([1.0 0.5; 0.0 1.0]))
    rich = repr(MIME"text/plain"(), t)
    @test occursin("→ array: ", rich) && occursin("→ basis: ", rich) && occursin("→ var: ", rich)
    @test occursin("→ KM: ", repr(MIME"text/plain"(), Tens(get_array(TensISO{3}(3.0, 2.0)))))
    @test occursin("→ dual basis: ", repr(MIME"text/plain"(), Basis([1.0 0.5; 0.0 1.0])))
end
