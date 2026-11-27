# ─────────────────────────────────────────────────────────────────────────────
#  waves_thomsen.jl — Thomsen's parameters of a transversely isotropic medium
#  and the velocities of weak anisotropy.
#
#  With the components of the stiffness in the two-index notation of the axis
#  3 and the density ρ, Thomsen (1986) describes the medium by the two
#  velocities along its axis and three dimensionless numbers,
#
#       V_P0 = √(C₃₃/ρ) ,  V_S0 = √(C₄₄/ρ) ,
#       ε = (C₁₁ − C₃₃)/(2C₃₃) ,  γ = (C₆₆ − C₄₄)/(2C₄₄) ,
#       δ = [(C₁₃ + C₄₄)² − (C₃₃ − C₄₄)²] / [2C₃₃(C₃₃ − C₄₄)] ,
#
#  exact definitions, whatever the strength of the anisotropy. They are a
#  `form` of `ti_constants` / `ti_from_constants` (`tens_anisotropic.jl`),
#  with the density as a keyword.
#
#  To first order in (ε, δ, γ) the exact velocities of `wave_modes`
#  (`waves_christoffel.jl`) become, θ being the angle to the axis,
#
#       V_P  ≈ V_P0 (1 + δ sin²θ cos²θ + ε sin⁴θ) ,
#       V_SV ≈ V_S0 [1 + (V_P0/V_S0)² (ε − δ) sin²θ cos²θ] ,
#       V_SH ≈ V_S0 (1 + γ sin²θ) ,
#
#  the expansions that the tests recover as series of the exact forms.
# ─────────────────────────────────────────────────────────────────────────────

# (V_P0, V_S0, ε, δ, γ) of the components (C₁₁, C₁₂, C₁₃, C₃₃, C₄₄) and ρ.
function _thomsen_of_components(C11, C12, C13, C33, C44, ρ)
    C66 = (C11 - C12) / 2
    δ = ((C13 + C44)^2 - (C33 - C44)^2) / (2 * C33 * (C33 - C44))
    return (sqrt(C33 / ρ), sqrt(C44 / ρ), (C11 - C33) / (2 * C33), δ, (C66 - C44) / (2 * C44))
end

# The converse. δ fixes (C₁₃ + C₄₄)² only: the root kept is the positive one,
# C₁₃ + C₄₄ > 0, which every stable medium met in practice satisfies.
function _components_of_thomsen(V_P0, V_S0, ε, δ, γ, ρ)
    C33, C44 = ρ * V_P0^2, ρ * V_S0^2
    C11, C66 = C33 * (1 + 2 * ε), C44 * (1 + 2 * γ)
    C13 = sqrt(2 * δ * C33 * (C33 - C44) + (C33 - C44)^2) - C44
    return (C11, C11 - 2 * C66, C13, C33, C44)
end

@noinline _throw_no_density() = throw(
    ArgumentError("Thomsen's parameters hold two velocities: give the mass density, `form = :thomsen, density = ρ`")
)
_thomsen_density(::Nothing) = _throw_no_density()
_thomsen_density(ρ) = ρ

# The stiffness written as a transversely isotropic one, which Thomsen's
# parameters need: a TensTI as it is; an isotropic one about any axis (its ε, δ
# and γ vanish), e₃ when none is given; any other when it belongs to the class,
# about the axis given or, numbers only, about the axis found.
_as_ti(t::TensTI{4}, ::Nothing) = t
_as_ti(t::TensTI{4}, n) = isparallel(Tuple(n), axis(t)) ? t : _throw_not_ti(t, n)
_as_ti(t::TensISO{4, 3}, n) = iso_to_ti(t, something(n, (0, 0, 1)))
function _as_ti(t::AbstractTens{4, 3}, n)
    ti = n === nothing ? project_exact(t, :TI) : project_exact(t, :TI, n)
    return ti === nothing ? _throw_not_ti(t, n) : ti
end
@noinline _throw_not_ti(t, axis) = throw(
    ArgumentError(
        "Thomsen's parameters are those of a transversely isotropic medium, and this " *
            "$(nameof(typeof(t))) is not one" * (axis === nothing ? "" : " about the axis $(Tuple(axis))") *
            ". For any symmetry, use generalized_thomsen_parameters (Mensch and Rasolofosaon)"
    )
)

"""
    thomsen_parameters(t::TensTI{4}, ρ) → (V_P0, V_S0, ε, δ, γ)
    thomsen_parameters(t::AbstractTens{4,3}, ρ; axis = nothing)

Thomsen's parameters of the transversely isotropic stiffness `t` of a medium
of density `ρ`, as a named tuple: the velocities of the P and S waves along
the axis, and

```math
\\varepsilon = \\frac{C_{11} - C_{33}}{2C_{33}},
\\qquad
\\gamma = \\frac{C_{66} - C_{44}}{2C_{44}},
\\qquad
\\delta = \\frac{(C_{13} + C_{44})^2 - (C_{33} - C_{44})^2}{2C_{33}\\,(C_{33} - C_{44})} ,
```

in the two-index notation of the axis 3. The definitions are exact; the
parameters are small for a weakly anisotropic medium
([Thomsen's parameters](@ref th-thomsen)). It is
`ti_constants(t; form = :thomsen, density = ρ)` with names; the converse is
`ti_from_constants(V_P0, V_S0, ε, δ, γ, axis; form = :thomsen, density = ρ)`.

Any other stiffness is accepted when it is transversely isotropic: an
isotropic one about any `axis` (``e_3`` by default), whose ``\\varepsilon``,
``\\delta`` and ``\\gamma`` vanish; a tensor of another type about the `axis` given
or, for numbers, about the axis found ([`project_exact`](@ref)). Otherwise the
call is refused: the generalized parameters of
[`generalized_thomsen_parameters`](@ref) hold for any symmetry.
"""
thomsen_parameters(t::AbstractTens{4, 3}, ρ; axis = nothing) = _thomsen_parameters(_as_ti(t, axis), ρ)
_thomsen_parameters(t::TensTI{4}, ρ) =
    NamedTuple{(:V_P0, :V_S0, :ε, :δ, :γ)}(_thomsen_of_components(_ti_components(t)..., ρ))

"""
    anellipticity(t::TensTI{4}) → η
    anellipticity(t::AbstractTens{4,3}; axis = nothing)

The anellipticity ``\\eta = (\\varepsilon - \\delta)/(1 + 2\\delta)`` of a
transversely isotropic stiffness: zero when the slowness curve of the P wave
is an ellipse. It does not depend on the density. Another stiffness is taken
as [`thomsen_parameters`](@ref) takes it.
"""
anellipticity(t::AbstractTens{4, 3}; axis = nothing) = _anellipticity(_as_ti(t, axis))
function _anellipticity(t::TensTI{4})
    _, _, ε, δ, _ = _thomsen_of_components(_ti_components(t)..., one(eltype(t)))
    return (ε - δ) / (1 + 2 * δ)
end

"""
    weak_wave_modes(t::TensTI{4}, ρ, n) → (qP = V, qSV = V, SH = V)
    weak_wave_modes(t::AbstractTens{4,3}, ρ, n; axis = nothing)

The phase velocities of the three modes along `n` to first order in Thomsen's
parameters,

```math
\\begin{gathered}
V_{P} \\approx V_{P0}\\,(1 + \\delta\\sin^2\\theta\\cos^2\\theta + \\varepsilon\\sin^4\\theta),
 \\\\
V_{SV} \\approx V_{S0}\\Bigl[1 + \\frac{V_{P0}^2}{V_{S0}^2}(\\varepsilon - \\delta)\\sin^2\\theta\\cos^2\\theta\\Bigr],
 \\\\
V_{SH} \\approx V_{S0}\\,(1 + \\gamma\\sin^2\\theta) ,
\\end{gathered}
```

``\\theta`` being the angle between `n` and the axis. [`wave_modes`](@ref)
gives the exact velocities, of which these are the first-order series.
Another stiffness is taken as [`thomsen_parameters`](@ref) takes it.
"""
weak_wave_modes(t::AbstractTens{4, 3}, ρ, n; axis = nothing) = _weak_wave_modes(_as_ti(t, axis), ρ, n)
function _weak_wave_modes(t::TensTI{4}, ρ, n)
    V_P0, V_S0, ε, δ, γ = _thomsen_of_components(_ti_components(t)..., ρ)
    v, a = _extract_vec(n), axis(t)
    c² = (v[1] * a[1] + v[2] * a[2] + v[3] * a[3])^2 / (v[1]^2 + v[2]^2 + v[3]^2)
    s² = 1 - c²
    return (
        qP = V_P0 * (1 + δ * s² * c² + ε * s²^2),
        qSV = V_S0 * (1 + (V_P0 / V_S0)^2 * (ε - δ) * s² * c²),
        SH = V_S0 * (1 + γ * s²),
    )
end

"""
    squared_moveout_velocities(t::TensTI{4}, ρ) → (qP = V², qSV = V², SH = V²)
    squared_moveout_velocities(t::AbstractTens{4,3}, ρ; axis = nothing)

The squares of the short-spread moveout velocities of a horizontal reflector
in a medium of vertical axis, exact for any strength of the anisotropy:
``V_{P0}^2(1 + 2\\delta)``,
``V_{S0}^2\\,[1 + 2(V_{P0}/V_{S0})^2(\\varepsilon - \\delta)]`` and
``V_{S0}^2(1 + 2\\gamma)``, that is the squared velocity of each mode along its
axis times one plus the second derivative, at the axis, of the logarithm of
its phase velocity ([Thomsen's parameters](@ref th-thomsen)). The squares are
returned because the qSV one is negative in many rocks
(``\\varepsilon - \\delta`` large and negative): the travel time of the
reflected qSV wave then decreases with the offset, and no moveout velocity
exists. Another stiffness is taken as [`thomsen_parameters`](@ref) takes it.
"""
squared_moveout_velocities(t::AbstractTens{4, 3}, ρ; axis = nothing) = _squared_moveout_velocities(_as_ti(t, axis), ρ)
function _squared_moveout_velocities(t::TensTI{4}, ρ)
    V_P0, V_S0, ε, δ, γ = _thomsen_of_components(_ti_components(t)..., ρ)
    return (qP = V_P0^2 * (1 + 2 * δ), qSV = V_S0^2 + 2 * V_P0^2 * (ε - δ), SH = V_S0^2 * (1 + 2 * γ))
end

export thomsen_parameters, anellipticity, weak_wave_modes, squared_moveout_velocities
