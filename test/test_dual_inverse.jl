import Symbolics

# The inverse of a generic 4th-order tensor of duals: by the implicit-function
# theorem on the values, never by a factorization of the duals.

@testsection "Inverse of a 4th-order tensor of duals" begin
    rng = Xoshiro(20261024)
    # A(x) = A₀ + x B in Mandel form, A₀ well conditioned.
    G = randn(rng, 6, 6)
    A₀ = G * transpose(G) / 6 + 3I
    B = randn(rng, 6, 6); B = B + transpose(B)
    sym4(M) = frommandel(SymmetricTensor{4, 3}, M)
    tens_of(x) = Tens(sym4(A₀ + x * B))
    mandel_inv(x) = tomandel(get_array(inv(tens_of(x))))

    Ai = inv(A₀)
    # First derivative against d(A⁻¹)/dx = -A⁻¹ B A⁻¹.
    @test ForwardDiff.derivative(mandel_inv, 0.0) ≈ -Ai * B * Ai rtol = 1.0e-12
    # Second derivative, through nested duals, against 2 A⁻¹ B A⁻¹ B A⁻¹.
    d2 = ForwardDiff.derivative(x -> ForwardDiff.derivative(mandel_inv, x), 0.0)
    @test d2 ≈ 2 * Ai * B * Ai * B * Ai rtol = 1.0e-11
    # Several partials at once: the gradient along two directions.
    B₂ = randn(rng, 6, 6); B₂ = B₂ + transpose(B₂)
    J = ForwardDiff.jacobian(p -> vec(tomandel(get_array(inv(Tens(sym4(A₀ + p[1] * B + p[2] * B₂)))))), [0.0, 0.0])
    @test reshape(J[:, 1], 6, 6) ≈ -Ai * B * Ai rtol = 1.0e-12
    @test reshape(J[:, 2], 6, 6) ≈ -Ai * B₂ * Ai rtol = 1.0e-12

    # Without minor symmetries: the Voigt (9×9) form.
    H = randn(rng, 9, 9) + 4I
    K = randn(rng, 9, 9)
    full_inv(x) = tovoigt(get_array(inv(Tens(fromvoigt(Tensor{4, 3}, H + x * K)))))
    Hi = inv(H)
    @test ForwardDiff.derivative(full_inv, 0.0) ≈ -Hi * K * Hi rtol = 1.0e-12

    # The value and the type: as the plain inverse, and inferred.
    td = Tens(sym4(A₀ + ForwardDiff.Dual{Nothing}(0.0, 1.0) * B))
    @test ForwardDiff.value.(tomandel(get_array(inv(td)))) ≈ Ai rtol = 1.0e-13
    @test get_array(@inferred(inv(td))) isa SymmetricTensor{4, 3, <:ForwardDiff.Dual}
end

@testsection "inverse of a fourth-order tensor of symbolic components" begin
    # Tensors.jl factorizes with pivoting, which compares the components; the
    # Mandel (or Voigt) matrix is inverted instead, against the inverse of the
    # same tensor at numbers.
    Symbolics.@variables k μ c
    entries(a, b, d) = (i, j, m, n) -> (i == m && j == n ? b : 0) + (i == j && m == n ? a : 0) + (i == j == m == n == 3 ? d : 0)
    at = Dict(k => 2.0, μ => 1.5, c => 0.7)
    num(x) = Symbolics.value(Symbolics.substitute(x, at; fold = Val(true)))
    ref = Array(get_array(inv(Tens(SymmetricTensor{4, 3}(entries(2.0, 1.5, 0.7))))))
    Ci = inv(Tens(SymmetricTensor{4, 3}(entries(k, μ, c))))
    @test get_array(Ci) isa SymmetricTensor{4, 3, Symbolics.Num}
    @test num.(Array(get_array(Ci))) ≈ ref
    # Full storage, without the minor symmetries: the Voigt matrix.
    full(a, b) = Tensor{4, 3}((i, j, m, n) -> (i == m && j == n ? b : 0) + (i == j && m == n ? a : 0) + (i == 1 && j == 2 && m == 1 && n == 2 ? a / 7 : 0))
    Fi = inv(Tens(full(k, μ)))
    @test num.(Array(get_array(Fi))) ≈ Array(get_array(inv(Tens(full(2.0, 1.5)))))
    # SymPy.
    sk, sμ, sc = symbols("k μ c", positive = true)
    Si = inv(Tens(SymmetricTensor{4, 3}(entries(sk, sμ, sc))))
    @test [Float64(N(subs(x, sk => 2.0, sμ => 1.5, sc => 0.7))) for x in Array(get_array(Si))] ≈ ref
    # Numbers keep the inverse of Tensors.jl, to the bit.
    A = SymmetricTensor{4, 3}(entries(2.0, 1.5, 0.7))
    @test TensND._inv_array(A) === inv(A)
    F = full(2.0, 1.5)
    @test TensND._inv_full(F) === inv(F)
end
