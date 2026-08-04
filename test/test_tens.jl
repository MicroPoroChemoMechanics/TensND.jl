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
