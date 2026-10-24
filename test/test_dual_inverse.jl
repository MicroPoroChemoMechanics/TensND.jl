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
