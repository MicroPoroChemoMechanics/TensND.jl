# ─────────────────────────────────────────────────────────────────────────────
#  waves_weak_anisotropy.jl — plane waves to first order in the deviation
#  from isotropy, for a stiffness of any symmetry.
#
#  Let ℂ = ℂ⁰ + δℂ with ℂ⁰ isotropic. The acoustic tensor of ℂ⁰ along n has
#  the eigenvalue λ + 2μ on n and the double eigenvalue μ on the plane normal
#  to n. To first order in δℂ, perturbation theory keeps these two subspaces:
#  the P eigenvalue becomes n·Γ·n, and the two S eigenvalues those of the
#  restriction of Γ to the plane normal to n, a 2×2 symmetric problem solved
#  in closed form (the degenerate case: the perturbation chooses the two
#  polarizations). Γ = Γ⁰ + δΓ being the acoustic tensor of ℂ itself, the
#  result does not depend on the isotropic tensor ℂ⁰ chosen:
#
#       ρV_P² ≈ n·Γ(n)·n = (n⊗n):ℂ:(n⊗n) ,
#       ρV_S² ≈ ½ [Γ₁₁ + Γ₂₂ ± √((Γ₁₁ − Γ₂₂)² + 4Γ₁₂²)] ,  Γ_ab = e_a·Γ(n)·e_b ,
#
#  with (e₁, e₂) any orthonormal pair normal to n. No velocity is compared
#  with another: symbols go through. The error is of second order in δℂ.
# ─────────────────────────────────────────────────────────────────────────────

"""
    weak_plane_waves(ℂ::AbstractTens{4,3}, ρ, n) → (qP = …, qS1 = …, qS2 = …)

The three plane waves along `n` to first order in the deviation of `ℂ` from
isotropy, each as a pair `(velocity, polarization)`:

```math
\\rho V_P^2 \\approx \\underline{n}\\cdot\\boldsymbol{\\Gamma}\\cdot\\underline{n},
\\qquad
\\rho V_{S1,S2}^2 \\approx \\tfrac12\\Bigl[\\Gamma_{11} + \\Gamma_{22} \\pm \\sqrt{(\\Gamma_{11} - \\Gamma_{22})^2 + 4\\,\\Gamma_{12}^2}\\Bigr],
```

where ``\\Gamma_{ab} = \\underline{e}_a\\cdot\\boldsymbol{\\Gamma}(\\underline{n})\\cdot\\underline{e}_b``
on an orthonormal pair normal to ``\\underline{n}``. The P wave is polarized
along ``\\underline{n}``, the S waves along the eigenvectors of the 2×2 problem
([Weak anisotropy of any symmetry](@ref th-weak-anisotropy)). The result does
not depend on a reference medium. Numbers, dual numbers and symbols;
[`plane_waves`](@ref) gives the exact waves.
"""
function weak_plane_waves(ℂ::AbstractTens{4, 3}, ρ, n)
    C = canonical_components(ℂ)
    ℬ = frame_from_axis(n)
    e = vecbasis(ℬ)
    e₁, e₂, e₃ = SVector{3}(e[:, 1]), SVector{3}(e[:, 2]), SVector{3}(e[:, 3])   # e₃ along n
    Γ = _christoffel(C, e₃)
    Γ₁₁, Γ₂₂, Γ₁₂ = _quadratic(Γ, e₁, e₁), _quadratic(Γ, e₂, e₂), _quadratic(Γ, e₁, e₂)
    r = sqrt((Γ₁₁ - Γ₂₂)^2 + 4 * Γ₁₂^2)
    λ₁, λ₂ = (Γ₁₁ + Γ₂₂ + r) / 2, (Γ₁₁ + Γ₂₂ - r) / 2
    # Eigenvectors of the 2×2 problem, not normalized: (Γ₁₂, λ − Γ₁₁) on (e₁, e₂).
    return (
        qP = (sqrt(_quadratic(Γ, e₃, e₃) / ρ), e₃),
        qS1 = (sqrt(λ₁ / ρ), Γ₁₂ * e₁ + (λ₁ - Γ₁₁) * e₂),
        qS2 = (sqrt(λ₂ / ρ), Γ₁₂ * e₁ + (λ₂ - Γ₁₁) * e₂),
    )
end

# a·Γ·b without the conjugation of `dot`, which symbols do not take.
_quadratic(Γ, a, b) = sum(a[i] * Γ[i, j] * b[j] for i in 1:3, j in 1:3)

"""
    generalized_thomsen_parameters(ℂ::AbstractTens{4,3}) → NamedTuple

The dimensionless parameters that describe, to first order, the velocities of
a stiffness of any symmetry about the vertical ``\\underline{e}_3``, in the
two-index notation of the canonical components
([`voigt_stiffness`](@ref)) and with the reference moduli of the vertical
waves, ``M = C_{33}`` for P and ``C_{44}``, ``C_{55}`` for the two S:

```math
\\varepsilon_x = \\frac{C_{22} - M}{2M},\\quad
\\varepsilon_y = \\frac{C_{11} - M}{2M},\\quad
\\delta_x = \\frac{C_{23} + 2C_{44} - M}{M},\\quad
\\delta_y = \\frac{C_{13} + 2C_{55} - M}{M},\\quad
\\delta_z = \\frac{C_{12} + 2C_{66} - M}{M},
```

```math
\\chi_x = \\frac{C_{14} + 2C_{56}}{M},\\quad
\\chi_y = \\frac{C_{25} + 2C_{46}}{M},\\quad
\\chi_z = \\frac{C_{36} + 2C_{45}}{M},\\quad
\\varepsilon_{IJ} = \\frac{C_{IJ}}{M}\\ \\ (IJ = 15, 16, 24, 26, 34, 35),
```

```math
\\gamma_x = \\frac{C_{44} - C_{55}}{2C_{55}},\\quad
\\gamma_y = \\frac{C_{66} - C_{44}}{2C_{44}},\\quad
\\gamma_z = \\frac{C_{55} - C_{44}}{2C_{44}} .
```

The fifteen first ones are the combinations of the stiffness that the P wave
sees to first order; the index names the plane normal to that axis. For a
transversely isotropic medium of axis ``\\underline{e}_3`` they reduce to
Thomsen's: ``\\varepsilon_x = \\varepsilon_y = \\varepsilon``,
``\\gamma_y = \\gamma``, and ``\\delta_x = \\delta_y`` equals Thomsen's
``\\delta`` to first order
([Weak anisotropy of any symmetry](@ref th-weak-anisotropy)). Numbers, dual
numbers and symbols.
"""
function generalized_thomsen_parameters(ℂ::AbstractTens{4, 3})
    C = voigt_stiffness(ℂ)
    M, C44, C55, C66 = C[3, 3], C[4, 4], C[5, 5], C[6, 6]
    return (
        εx = (C[2, 2] - M) / (2 * M), εy = (C[1, 1] - M) / (2 * M),
        δx = (C[2, 3] + 2 * C44 - M) / M, δy = (C[1, 3] + 2 * C55 - M) / M, δz = (C[1, 2] + 2 * C66 - M) / M,
        χx = (C[1, 4] + 2 * C[5, 6]) / M, χy = (C[2, 5] + 2 * C[4, 6]) / M, χz = (C[3, 6] + 2 * C[4, 5]) / M,
        ε15 = C[1, 5] / M, ε16 = C[1, 6] / M, ε24 = C[2, 4] / M, ε26 = C[2, 6] / M, ε34 = C[3, 4] / M, ε35 = C[3, 5] / M,
        γx = (C44 - C55) / (2 * C55), γy = (C66 - C44) / (2 * C44), γz = (C55 - C44) / (2 * C44),
    )
end

export weak_plane_waves, generalized_thomsen_parameters
