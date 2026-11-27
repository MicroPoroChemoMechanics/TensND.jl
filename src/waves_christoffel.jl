# ─────────────────────────────────────────────────────────────────────────────
#  waves_christoffel.jl — plane waves of a homogeneous elastic medium.
#
#  A plane wave u = U f(n·x − V t) of unit direction n, phase velocity V and
#  polarization U solves, in a medium of stiffness ℂ and mass density ρ,
#
#       (Γ(n) − ρV² 𝟏)·U = 0 ,     Γ(n) = n·ℂ·n ,
#
#  the Christoffel problem: ρV² is an eigenvalue of the acoustic tensor
#  (`acoustic_tensor`, `km_blocks.jl`) and U its eigenvector. Three waves go
#  in each direction, the fastest quasi-longitudinal (qP), the two others
#  quasi-transverse (qS1 faster than qS2).
#
#  The energy travels at the group velocity g = ∂ω/∂k, which is not along n in
#  an anisotropic medium. With a unit polarization,
#
#       g = (U·ℂ·U)·n / (ρV) ,
#
#  that is the acoustic tensor of the direction U applied to n; and g·n = V.
#
#  - `plane_waves` solves the problem for any stiffness, on numbers and dual
#    numbers: the three waves sorted by decreasing velocity. Sorting compares
#    values, which symbols do not take.
#  - `wave_modes` gives the closed forms of the classes that have them,
#    isotropic and transversely isotropic, each mode under its name and
#    without any comparison: symbols go through.
#  - `christoffel_invariants` gives the coefficients of the characteristic
#    polynomial of Γ/ρ, polynomial in the stiffness and the direction.
# ─────────────────────────────────────────────────────────────────────────────

"""
    PlaneWaves{T}

The three plane waves of a direction `direction` (unit), sorted by decreasing
phase velocity (qP, qS1, qS2): `velocity[m]` the phase velocity of the wave
`m`, `polarization[:, m]` its unit polarization and `group[:, m]` its group
velocity. Built by [`plane_waves`](@ref).
"""
struct PlaneWaves{T}
    direction::SVector{3, T}
    velocity::SVector{3, T}
    polarization::SMatrix{3, 3, T, 9}
    group::SMatrix{3, 3, T, 9}
end

function Base.show(io::IO, w::PlaneWaves)
    print(io, "PlaneWaves along ", Tuple(w.direction), ": V = ", Tuple(w.velocity))
    return nothing
end

# Unit direction at the element type `T`; a zero direction is refused.
function _unit_direction(n, ::Type{T}) where {T}
    v = _direction_components(n)
    nn = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
    _refuse_zero(nn, "the direction of propagation is zero")
    return SVector{3, T}(v[1] / nn, v[2] / nn, v[3] / nn)
end

# The three components of a direction of space. The length is checked first:
# `_extract_vec` reads the first three components of whatever it is given.
function _direction_components(n)
    length(n) == 3 || throw(DimensionMismatch("a direction of space has three components; got $(length(n))"))
    return _extract_vec(n)
end

_sym3(K11, K22, K33, K12, K13, K23) = SMatrix{3, 3}(K11, K12, K13, K12, K22, K23, K13, K23, K33)

# Γ(ζ) of an array of canonical components.
_christoffel(C::AbstractArray{<:Any, 4}, ζ) = _sym3(acoustic_components(C, ζ)...)

# Eigenvalues, ascending, and eigenvectors of a symmetric 3×3 matrix: the
# closed form of StaticArrays on floats, the differentiated decomposition of
# `_sym_eigen` (`tens_projection.jl`) on dual numbers.
_eigen3(Γ::SMatrix{3, 3, T}) where {T <: AbstractFloat} = (e = eigen(Symmetric(Γ)); (SVector{3, T}(e.values), SMatrix{3, 3, T}(e.vectors)))
function _eigen3(Γ::SMatrix{3, 3, T}) where {T <: ForwardDiff.Dual}
    λ, V = _sym_eigen(Γ)
    return SVector{3, T}(λ), SMatrix{3, 3, T}(V)
end
_eigen3(::SMatrix{3, 3, T}) where {T} = throw(
    ArgumentError(
        "plane_waves sorts the three waves by their velocities, which a $T does not allow; " *
            "wave_modes gives the closed forms of the isotropic and transversely isotropic classes, " *
            "christoffel_invariants the characteristic polynomial of any tensor"
    )
)

@noinline _throw_unstable(λ) = throw(
    DomainError(
        primal(λ),
        "the acoustic tensor has an eigenvalue that is not positive: no wave travels in that direction " *
            "(the stiffness is not positive definite, or the density is not positive)"
    )
)

"""
    plane_waves(ℂ::AbstractTens{4,3}, ρ, n) → PlaneWaves

The three plane waves that travel along the direction `n` (any nonzero
vector) in a medium of stiffness `ℂ` and mass density `ρ`: solutions of the
Christoffel problem

```math
\\bigl(\\boldsymbol{\\Gamma}(\\underline{n}) - \\rho V^2\\,\\boldsymbol{1}\\bigr)\\cdot\\underline{U} = \\underline{0},
\\qquad \\boldsymbol{\\Gamma}(\\underline{n}) = \\underline{n}\\cdot\\mathbb{C}\\cdot\\underline{n} ,
```

sorted by decreasing phase velocity (qP, qS1, qS2), with their unit
polarizations and their group velocities
``\\underline{g} = (\\underline{U}\\cdot\\mathbb{C}\\cdot\\underline{U})\\cdot\\underline{n}/(\\rho V)``
([Plane waves](@ref th-plane-waves)). The components are canonical, whatever
the basis `ℂ` is written in.

Numbers and dual numbers; the derivative of a polarization is left out
between two waves of the same velocity, where it is not defined. A direction
in which an eigenvalue is not positive throws a `DomainError`. For symbols see
[`wave_modes`](@ref) and [`christoffel_invariants`](@ref).

# Examples
```jldoctest
julia> w = plane_waves(iso_stiffness(5.0, 2.0), 2.0, (0, 0, 1));

julia> w.velocity ≈ [sqrt((5 + 4 * 2 / 3) / 2), 1, 1]
true
```
"""
function plane_waves(ℂ::AbstractTens{4, 3}, ρ, n)
    C = canonical_components(ℂ)
    T = float_like(promote_type(eltype(C), typeof(ρ), map(typeof, _direction_components(n))...))
    return _plane_waves(C, convert(T, ρ), _unit_direction(n, T))
end

function _plane_waves(C::AbstractArray{<:Any, 4}, ρ::T, n::SVector{3, T}) where {T}
    λ, U = _eigen3(SMatrix{3, 3, T}(_christoffel(C, n)) / ρ)
    primal(λ[1]) > 0 || _throw_unstable(λ[1])
    # ascending eigenvalues: the waves are taken from the last to the first
    V = SVector{3, T}(sqrt(λ[3]), sqrt(λ[2]), sqrt(λ[1]))
    P = SMatrix{3, 3, T}(U[:, 3]..., U[:, 2]..., U[:, 1]...)
    g = map(m -> _group_velocity(C, ρ, n, V[m], P[:, m]), SVector(1, 2, 3))
    return PlaneWaves{T}(n, V, P, SMatrix{3, 3, T}(g[1]..., g[2]..., g[3]...))
end

# g = (U·ℂ·U)·n / (ρV), the polarization being a unit vector.
_group_velocity(C, ρ, n, V, U) = (SMatrix{3, 3}(_christoffel(C, U)) * n) / (ρ * V)

"""
    phase_velocities(ℂ, ρ, n) → SVector{3}

The phase velocities of the three plane waves along `n`, in decreasing order
(qP, qS1, qS2): the field `velocity` of [`plane_waves`](@ref).
"""
phase_velocities(ℂ::AbstractTens{4, 3}, ρ, n) = plane_waves(ℂ, ρ, n).velocity

"""
    slowness(w::PlaneWaves) → SMatrix{3,3}

The slowness vectors ``\\underline{n}/V`` of the three waves, one per column.
"""
slowness(w::PlaneWaves) = w.direction * transpose(inv.(w.velocity))

"""
    christoffel_invariants(ℂ::AbstractTens{4,3}, ρ, n) → (I₁, I₂, I₃)

The invariants of ``\\boldsymbol{\\Gamma}(\\underline{n})/\\rho``: its trace, the
trace of its adjugate and its determinant. The squared phase velocities
``v = V^2`` along the unit direction `n` are the roots of

```math
v^3 - I_1\\,v^2 + I_2\\,v - I_3 = 0 .
```

Polynomial in the components of `ℂ` and of `n`: numbers, dual numbers and
symbols.
"""
function christoffel_invariants(ℂ::AbstractTens{4, 3}, ρ, n)
    K = map(k -> k / ρ, acoustic_components(canonical_components(ℂ), _extract_vec(n)))
    A11, A22, A33, _, _, _ = sym3_adj(K...)
    K11, K22, K33, K12, K13, K23 = K
    det = K11 * A11 + K12 * (K13 * K23 - K12 * K33) + K13 * (K12 * K23 - K13 * K22)
    return K11 + K22 + K33, A11 + A22 + A33, det
end

# ── Closed forms, mode by mode ───────────────────────────────────────────────

"""
    wave_modes(ℂ::TensISO{4,3}, ρ, n) → (P = (velocity, polarization), S = …)
    wave_modes(ℂ::TensTI{4}, ρ, n) → (qP = …, qSV = …, SH = …)

The plane waves of an isotropic or of a transversely isotropic medium along
`n`, in closed form and under the name of their mode, each as a pair
`(velocity, polarization)`. No velocity is compared with another: symbols go
through, and a mode keeps its name where it crosses another.

Isotropic: ``\\rho V_P^2 = \\lambda + 2\\mu`` along ``\\underline{n}``, and
``\\rho V_S^2 = \\mu`` for every polarization normal to it (`polarization` is
then `nothing`).

Transversely isotropic of axis ``\\underline{a}``, with ``\\theta`` the angle
between ``\\underline{n}`` and ``\\underline{a}`` and the components
``C_{11}, C_{33}, C_{13}, C_{44}, C_{66}`` in the usual two-index notation of
the axis 3:

```math
\\rho V_{SH}^2 = C_{66}\\sin^2\\theta + C_{44}\\cos^2\\theta ,
\\qquad
2\\rho V_{qP,\\,qSV}^2 = C_{11}\\sin^2\\theta + C_{33}\\cos^2\\theta + C_{44} \\pm \\sqrt{D} ,
```

```math
D = \\bigl[(C_{11} - C_{44})\\sin^2\\theta - (C_{33} - C_{44})\\cos^2\\theta\\bigr]^2 + 4\\,(C_{13} + C_{44})^2\\sin^2\\theta\\cos^2\\theta .
```

SH is polarized along ``\\underline{a}\\times\\underline{n}``; qP and qSV in the
plane of ``\\underline{a}`` and ``\\underline{n}``, given by their components
`(along the transverse direction, along the axis)` in that plane, for a
direction on the side of the axis, ``\\underline{n}\\cdot\\underline{a} \\geq 0``
(the direction ``-\\underline{n}`` carries the same waves)
([Plane waves](@ref th-plane-waves)).
"""
function wave_modes(ℂ::TensISO{4, 3}, ρ, n)
    α, β = get_data(ℂ)                         # 3k and 2μ
    return (P = (sqrt((α + 2 * β) / (3 * ρ)), _extract_vec(n)), S = (sqrt(β / (2 * ρ)), nothing))
end

function wave_modes(ℂ::TensTI{4}, ρ, n)
    C11, C12, C13, C33, C44 = _ti_components(ℂ)
    C66 = (C11 - C12) / 2
    v, a = _extract_vec(n), axis(ℂ)
    n² = v[1]^2 + v[2]^2 + v[3]^2
    c = (v[1] * a[1] + v[2] * a[2] + v[3] * a[3])        # |n| cos θ
    c², s² = c^2 / n², 1 - c^2 / n²
    A, B = C11 * s² + C44 * c², C44 * s² + C33 * c²      # the 2×2 problem in the plane (transverse, axis)
    D = (A - B)^2 + 4 * (C13 + C44)^2 * s² * c²
    rD = sqrt(D)
    H = (C13 + C44) * sqrt(s² * c²)                      # its coupling, up to the sign of sin θ cos θ
    qP, qSV = (A + B + rD) / 2, (A + B - rD) / 2
    SH = cross(SVector{3}(a), SVector{3}(v))
    return (
        qP = (sqrt(qP / ρ), (qP - B, H)), qSV = (sqrt(qSV / ρ), (qSV - B, H)),
        SH = (sqrt((C66 * s² + C44 * c²) / ρ), SH),
    )
end

export PlaneWaves, plane_waves, phase_velocities, slowness, christoffel_invariants, wave_modes
