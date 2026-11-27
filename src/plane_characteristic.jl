# ─────────────────────────────────────────────────────────────────────────────
#  plane_characteristic.jl — the characteristic roots of a stiffness of the
#  plane.
#
#  A displacement u = U f(x₁ + μ x₂) of the plane is in equilibrium in the
#  medium of stiffness ℂ when
#
#       [Q + μ (R + Rᵀ) + μ² T]·U = 0 ,
#       Q = e₁·ℂ·e₁ ,  R = e₁·ℂ·e₂ ,  T = e₂·ℂ·e₂ ,
#
#  the slowness problem of an interface (`waves_interface.jl`) without its
#  inertia term. The determinant is a quartic in μ with no real root for a
#  positive definite stiffness: two roots μ₁, μ₂ of positive imaginary part
#  and their conjugates. With the compliances of the plane in the engineering
#  form,
#
#       ε₁₁ = a₁₁σ₁₁ + a₁₂σ₂₂ + a₁₆σ₁₂ ,  ε₂₂ = a₁₂σ₁₁ + a₂₂σ₂₂ + a₂₆σ₁₂ ,
#       2ε₁₂ = a₁₆σ₁₁ + a₂₆σ₂₂ + a₆₆σ₁₂ ,
#
#  the quartic reads (Lekhnitskii)
#
#       a₁₁ μ⁴ − 2a₁₆ μ³ + (2a₁₂ + a₆₆) μ² − 2a₂₆ μ + a₂₂ = 0 .
#
#  An isotropic medium has the double root i. What is symmetric in (μ₁, μ₂)
#  stays smooth there: the sum s and the product p are what is computed, the
#  roots themselves being read from them.
# ─────────────────────────────────────────────────────────────────────────────

# z₁ / z₂ without the comparisons of the division of `Base`, which symbols do
# not take.
_cdiv(z₁, z₂) = (z₁ * conj(z₂)) / abs2(z₂)

# The compliances (a₁₁, a₂₂, a₁₂, a₆₆, a₁₆, a₂₆) read on the Kelvin-Mandel
# matrix `M` of the compliance of the plane.
function _engineering_compliances(M::AbstractMatrix)
    r = sqrt(2 * one(eltype(M)))
    return (M[1, 1], M[2, 2], M[1, 2], 2 * M[3, 3], r * M[1, 3], r * M[2, 3])
end

"""
    engineering_compliances(ℂ::AbstractTens{4,2}) → (a₁₁, a₂₂, a₁₂, a₆₆, a₁₆, a₂₆)

The compliances of the stiffness `ℂ` of the plane in the engineering form, in
the canonical axes:

```math
\\begin{gathered}
\\varepsilon_{11} = a_{11}\\sigma_{11} + a_{12}\\sigma_{22} + a_{16}\\sigma_{12}, \\\\
\\varepsilon_{22} = a_{12}\\sigma_{11} + a_{22}\\sigma_{22} + a_{26}\\sigma_{12}, \\\\
2\\varepsilon_{12} = a_{16}\\sigma_{11} + a_{26}\\sigma_{22} + a_{66}\\sigma_{12} .
\\end{gathered}
```
"""
engineering_compliances(ℂ::AbstractTens{4, 2}) = _engineering_compliances(mandel(inv(ℂ), CanonicalBasis{2, eltype(ℂ)}()))

# In the axes of an orthotropic medium a₁₆ = a₂₆ = 0 and the equation is
# biquadratic: p = −√(a₂₂/a₁₁), and s² = −(2a₁₂ + a₆₆)/a₁₁ + 2p, which is
# negative for a positive definite compliance.
function _material_pair(a₁₁, a₂₂, a₁₂, a₆₆)
    ρ = sqrt(a₂₂ / a₁₁)
    return complex(zero(ρ), sqrt((2 * a₁₂ + a₆₆) / a₁₁ + 2 * ρ)), complex(-ρ)
end

# The pair in axes from which the material ones are turned by the angle of
# cosine c and sine n: x₁′ + μ′x₂′ is proportional to x₁ + μ x₂ with
# μ = (μ′c + n)/(c − μ′n), whose sum and product are rational in (s′, p′).
function _turned_pair(s′, p′, c, n)
    d = c^2 - s′ * (c * n) + p′ * n^2
    return _cdiv(s′ * (c^2 - n^2) + 2 * (c * n) * (1 - p′), d), _cdiv(p′ * c^2 + s′ * (c * n) + n^2, d)
end

"""
    characteristic_pair(ℂ::AbstractTens{4,2}) → (s, p)

Sum ``s = \\mu_1 + \\mu_2`` and product ``p = \\mu_1\\mu_2`` of the two
characteristic roots of positive imaginary part of the stiffness `ℂ` of the
plane, in the canonical axes: the roots of

```math
a_{11}\\mu^4 - 2a_{16}\\mu^3 + (2a_{12} + a_{66})\\mu^2 - 2a_{26}\\mu + a_{22} = 0
```

with the engineering compliances ``a_{IJ}`` of the plane
([Characteristic roots](@ref th-characteristic-roots)).

An isotropic tensor has the double root ``\\mathrm{i}``, hence
``(2\\mathrm{i}, -1)``; an orthotropic or a square one has the pair of its own
axes in closed form, turned to the canonical ones: these three accept symbols.
Any other tensor has its equation solved numerically (numbers and dual
numbers), and must be major-symmetric. See [`characteristic_roots`](@ref).
"""
characteristic_pair(ℂ::TensISO{4, 2}) = (complex(zero(eltype(ℂ)), 2 * one(eltype(ℂ))), complex(-one(eltype(ℂ))))
function characteristic_pair(ℂ::Union{TensOrtho{4, 2, <:Any, 4}, TensCubic{2}})
    M = mandel_material(inv(ℂ))
    e = vecbasis(frame(ℂ))
    return _turned_pair(_material_pair(M[1, 1], M[2, 2], M[1, 2], 2 * M[3, 3])..., e[1, 1], e[2, 1])
end
characteristic_pair(ℂ::AbstractTens{4, 2}) = _characteristic_pair(Val(is_hard_numeric(eltype(ℂ))), ℂ)
_characteristic_pair(::Val{false}, ℂ) = throw(
    ArgumentError(
        "the characteristic equation of a symbolic $(nameof(typeof(ℂ))) of the plane is not solved in closed form; " *
            "give the tensor as a TensISO, a TensCubic or a TensOrtho"
    )
)
function _characteristic_pair(::Val{true}, ℂ)
    # On the values: the norm of a null difference has no derivative.
    M = primal.(mandel(ℂ, CanonicalBasis{2, eltype(ℂ)}()))
    asym = maximum(abs, M - transpose(M)) / maximum(abs, M)
    asym ≤ 1.0e-10 || throw(
        ArgumentError("the characteristic roots are those of a major-symmetric stiffness; got a relative asymmetry $asym")
    )
    return _characteristic_pair(engineering_compliances(ℂ))
end

# The quartic divided by a₁₁ is (μ² − sμ + p)(μ² − s̄μ + p̄): its coefficients
# give Re s and three equations on (Im s, Re p, Im p). The roots of the
# companion matrix, computed on the values, start a Newton iteration on these,
# which brings the pair to the precision of its type and carries the
# derivatives. The roots are sorted by a comparison: numbers only.
function _characteristic_pair(a::NTuple{6, T}) where {T}
    a₁₁, a₂₂, a₁₂, a₆₆, a₁₆, a₂₆ = a
    c₃, c₂, c₁, c₀ = -2 * a₁₆ / a₁₁, (2 * a₁₂ + a₆₆) / a₁₁, -2 * a₂₆ / a₁₁, a₂₂ / a₁₁
    f₃, f₂, f₁, f₀ = Float64.(primal.((c₃, c₂, c₁, c₀)))
    roots = filter(μ -> imag(μ) > 0, eigvals([0.0 0.0 0.0 -f₀; 1.0 0.0 0.0 -f₁; 0.0 1.0 0.0 -f₂; 0.0 0.0 1.0 -f₃]))
    length(roots) == 2 || throw(
        DomainError(
            (f₃, f₂, f₁, f₀),
            "the characteristic equation has real roots: the stiffness is not positive definite"
        )
    )
    sᵣ, c = -c₃ / 2, (c₂, c₁, c₀)
    x = SVector{3, T}(imag(roots[1] + roots[2]), real(roots[1] * roots[2]), imag(roots[1] * roots[2]))
    # The start has half the digits at worst (the companion matrix of a double
    # root) and each step doubles them: the values reach the precision of
    # `Float64` in two steps, and each further step makes one more order of
    # derivatives exact.
    for _ in 1:5
        x -= _pair_correction(x, sᵣ, c)
    end
    return complex(sᵣ, x[1]), complex(x[2], x[3])
end

# Newton correction of x = (Im s, Re p, Im p) on the three equations.
function _pair_correction(x::SVector{3, T}, sᵣ, (c₂, c₁, c₀)) where {T}
    sᵢ, pᵣ, pᵢ = x
    F = SVector{3, T}(sᵣ^2 + sᵢ^2 + 2 * pᵣ - c₂, 2 * (sᵣ * pᵣ + sᵢ * pᵢ) + c₁, pᵣ^2 + pᵢ^2 - c₀)
    J = SMatrix{3, 3, T}(2 * sᵢ, 2 * pᵢ, zero(T), 2 * one(T), 2 * sᵣ, 2 * pᵣ, zero(T), 2 * sᵢ, 2 * pᵢ)
    return J \ F
end

"""
    characteristic_roots(ℂ::AbstractTens{4,2}) → (μ₁, μ₂)

The two characteristic roots of positive imaginary part of the stiffness `ℂ`
of the plane, from their sum and product ([`characteristic_pair`](@ref)):
``\\mu_{1,2} = \\tfrac12\\big(s \\pm \\sqrt{s^2 - 4p}\\big)``. Near an
isotropic tensor the two roots merge and are not differentiable functions of
the stiffness, while their sum and product are: a computation that must stay
smooth is written with the pair.
"""
function characteristic_roots(ℂ::AbstractTens{4, 2})
    s, p = characteristic_pair(ℂ)
    r = sqrt(s^2 - 4 * p)
    return (s + r) / 2, (s - r) / 2
end

export engineering_compliances, characteristic_pair, characteristic_roots
