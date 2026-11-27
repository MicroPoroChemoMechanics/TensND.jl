# ─────────────────────────────────────────────────────────────────────────────
#  waves_interface.jl — plane waves at a plane interface: the law of
#  Snell-Descartes and the reflection and transmission coefficients.
#
#  A wave u = a U exp[iω(s·x − t)] has the slowness vector s. On a plane
#  interface of unit normal ν, all the waves that an incident one excites
#  share its tangential slowness s_t (the phases must match along the
#  interface: Snell-Descartes). In each medium the normal slowness ξ of
#  s = s_t + ξ ν solves the Christoffel problem
#
#       [Q + ξ (R + Rᵀ) + ξ² T − ρ 𝟏]·U = 0 ,
#       Q = s_t·ℂ·s_t ,  R = s_t·ℂ·ν ,  T = ν·ℂ·ν ,
#
#  a quadratic eigenproblem with six roots, real or in conjugate pairs. A real
#  root is a homogeneous wave, which carries energy along its group velocity;
#  a complex one an evanescent wave, bound to the interface. Three of the six
#  go (or decay) towards +ν, three towards −ν.
#
#  The traction on the interface is iω a τ with τ = (ν·ℂ·s)·U = (Rᵀ + ξ T)·U,
#  and the mean energy flux through it ½ω²|a|² Re(τ·Ū). A welded interface
#  keeps the displacement and the traction: six equations for the three
#  reflected and the three transmitted amplitudes.
#
#  The roots are those of a 6×6 matrix, found numerically and sorted by
#  comparisons: numbers, and dual numbers through the derivative of each
#  simple root; not symbols.
# ─────────────────────────────────────────────────────────────────────────────

"""
    InterfaceWave{T}

One of the six plane waves of a medium for a given tangential slowness:
`normal_slowness` (ξ, complex), `slowness` (``\\underline{s}_t + \\xi\\,\\underline{\\nu}``),
`polarization` (unit, complex), `traction` (``\\underline{\\tau} = (\\underline{\\nu}\\cdot\\mathbb{C}\\cdot\\underline{s})\\cdot\\underline{U}``),
`flux` (``\\operatorname{Re}(\\underline{\\tau}\\cdot\\overline{\\underline{U}})``, the mean
energy flux through the interface per unit ``\\tfrac12\\omega^2|a|^2``) and
`homogeneous` (true for a real ξ).
"""
struct InterfaceWave{T}
    normal_slowness::Complex{T}
    slowness::SVector{3, Complex{T}}
    polarization::SVector{3, Complex{T}}
    traction::SVector{3, Complex{T}}
    flux::T
    homogeneous::Bool
end

# Q, R, T of the header from canonical components.
function _interface_blocks(C::AbstractArray{<:Any, 4}, sₜ, ν)
    blk(a, b) = SMatrix{3, 3}(ntuple(q -> ((i, k) = (mod1(q, 3), cld(q, 3)); sum(a[j] * C[j, i, k, l] * b[l] for j in 1:3, l in 1:3)), Val(9)))
    return blk(sₜ, sₜ), blk(sₜ, ν), blk(ν, ν)
end

"""
    interface_waves(ℂ::AbstractTens{4,3}, ρ, ν, sₜ) → (forward, backward)

The six plane waves of the medium of stiffness `ℂ` and density `ρ` whose
slowness has the tangential part `sₜ` on the plane of unit normal `ν` (the
part of `sₜ` along `ν` is removed): three `forward` waves, which carry energy
or decay towards ``+\\underline{\\nu}``, and three `backward` ones, as tuples of
[`InterfaceWave`](@ref) sorted by increasing normal slowness in absolute
value (the fastest wave, quasi-longitudinal, first).

A homogeneous wave (real normal slowness) is classed by the sign of the normal
component of its energy flux; an evanescent one by the sign of the imaginary
part of its normal slowness, positive when it decays towards
``+\\underline{\\nu}``. The angle beyond which a wave turns evanescent is its
critical angle ([Waves at an interface](@ref th-interface-waves)).

Two homogeneous waves of the same normal slowness (the two shear waves of an
isotropic medium) share a plane of polarizations, of which any orthonormal pair
solves the problem: the pair returned is the one whose energy fluxes add. Their
individual amplitudes in [`reflection_transmission`](@ref) depend on that
choice; the field and the energy they carry together do not.

Numbers and dual numbers. An isotropic stiffness with symbols (`Symbolics.Num`)
gives the six waves in closed form, below the critical angles and for a
nonzero tangential slowness: no root is compared there, the forward waves are
those of positive normal slowness.
"""
function interface_waves(ℂ::AbstractTens{4, 3}, ρ, ν, sₜ)
    C = canonical_components(ℂ)
    T = float_like(promote_type(eltype(C), typeof(ρ), eltype(_extract_vec(ν)), eltype(_extract_vec(sₜ))))
    T <: _EigenType || throw(
        ArgumentError("interface_waves sorts the six roots of the slowness problem by comparisons: it takes numbers and dual numbers, not $T")
    )
    n = _unit_direction(ν, T)
    v = SVector{3, T}(_extract_vec(sₜ)...)
    s = v - (transpose(v) * n) * n
    return _interface_waves(C, convert(T, ρ), n, s)
end

# An isotropic medium with symbols: the six waves in closed form, the root
# being chosen without comparison. Below the critical angles the normal
# slownesses ξ_P = √(ρ/(λ+2μ) − p²) and ξ_S = √(ρ/μ − p²) are real, p = |sₜ|,
# and the forward waves are those of +ξ; the polarizations are along the
# slowness (P), in the plane of incidence normal to it (SV) and along ν × t
# (SH), t = sₜ/p. Numbers keep the general route of `interface_waves`.
interface_waves(ℂ::TensISO{4, 3}, ρ, ν, sₜ) =
    _iso_interface_waves(Val(is_hard_numeric(promote_type(eltype(ℂ), typeof(ρ), map(typeof, _extract_vec(ν))..., map(typeof, _extract_vec(sₜ))...))), ℂ, ρ, ν, sₜ)
_iso_interface_waves(::Val{true}, ℂ, ρ, ν, sₜ) = invoke(interface_waves, Tuple{AbstractTens{4, 3}, Any, Any, Any}, ℂ, ρ, ν, sₜ)
function _iso_interface_waves(::Val{false}, ℂ, ρ, ν, sₜ)
    α, β = get_data(ℂ)                              # 3k and 2μ
    M, μ = (α + 2 * β) / 3, β / 2                   # λ + 2μ and μ
    # The waves carry complex numbers, which need a real element type: the
    # symbols of Symbolics are, those of SymPy are not.
    M isa Real || throw(
        ArgumentError("interface_waves takes the symbols of Symbolics (`Num`), whose complex numbers exist, not $(typeof(M))")
    )
    n = _unit_direction(ν, typeof(M))
    v = SVector{3}(_extract_vec(sₜ)...)
    s = v - (transpose(v) * n) * n
    p = sqrt(transpose(s) * s)
    t = s * inv(p)                                  # `s / p` is a symbolic array for symbols
    ξP, ξS = sqrt(ρ / M - p^2), sqrt(ρ / μ - p^2)
    _, R, Tm = _interface_blocks(canonical_components(ℂ), s, n)
    c(x) = map(complex, x)
    wave(ξ, U) = _interface_wave(complex(ξ), c(U), R, Tm, s, n, true)
    waves(σ) = (
        wave(σ * ξP, (s + σ * ξP * n) * sqrt(M / ρ)),
        wave(σ * ξS, (σ * ξS * t - p * n) * sqrt(μ / ρ)),
        wave(σ * ξS, cross(n, t)),
    )
    return waves(1), waves(-1)
end

# One wave from its normal slowness and its unit polarization.
function _interface_wave(ξ::Complex{T}, U::SVector{3, Complex{T}}, R, Tm, sₜ, ν, homogeneous::Bool) where {T}
    τ = (transpose(R) + ξ * Tm) * U
    return InterfaceWave{T}(ξ, sₜ + ξ * ν, U, τ, real(transpose(τ) * conj.(U)), homogeneous)
end

function _interface_waves(C, ρ::T, ν::SVector{3, T}, sₜ::SVector{3, T}) where {T <: AbstractFloat}
    Q, R, Tm = _interface_blocks(C, sₜ, ν)
    Q̃ = Q - ρ * one(Q)
    Ti = inv(Tm)
    # ξ U as the second block: [0 𝟏; −T⁻¹Q̃ −T⁻¹(R + Rᵀ)]
    A = [zeros(T, 3, 3) Matrix{T}(I, 3, 3); -Matrix(Ti * Q̃) -Matrix(Ti * (R + transpose(R)))]
    e = eigen(A)
    scale = sqrt(ρ / maximum(abs, Tm))                       # a slowness of the medium
    tol = 1.0e-8 * scale
    ξ = Complex{T}.(e.values)
    hom = map(z -> abs(imag(z)) ≤ tol / 10, ξ)
    U = [SVector{3, Complex{T}}(e.vectors[1, m], e.vectors[2, m], e.vectors[3, m]) for m in 1:6]
    S = R + transpose(R)
    for m in 1:6
        group = [q for q in 1:6 if hom[m] && hom[q] && abs(ξ[q] - ξ[m]) ≤ tol]
        (length(group) > 1 && first(group) == m) && _split_coinciding!(ξ, U, group, S, Tm)
    end
    waves = map(1:6) do m
        _interface_wave(hom[m] ? complex(real(ξ[m])) : ξ[m], U[m] / sqrt(real(sum(abs2, U[m]))), R, Tm, sₜ, ν, hom[m])
    end
    isforward(w) = w.homogeneous ? w.flux > 0 : imag(w.normal_slowness) > 0
    forward, backward = filter(isforward, waves), filter(!isforward, waves)
    length(forward) == 3 || throw(
        DomainError(length(forward), "the slowness problem did not split into three waves each way: the stiffness is not positive definite, or a wave grazes the interface")
    )
    byspeed(w) = (abs(real(w.normal_slowness)), abs(imag(w.normal_slowness)))
    three(ws) = (v = sort(ws; by = byspeed); ntuple(i -> v[i], Val(3)))       # three each way, checked above
    return three(forward), three(backward)
end

# Real roots that coincide (the two shear waves of an isotropic medium) share
# an eigenspace, of which any basis solves the problem. The energy fluxes of
# two of its waves add only if τ_a·U_b + U_a·τ_b = U_a·P′·U_b vanishes, with
# P′ = R + Rᵀ + 2ξT: the polarizations kept are the eigenvectors of P′
# restricted to the eigenspace, orthonormal and orthogonal for the flux.
function _split_coinciding!(ξ, U, group, S, Tm)
    g = length(group)
    ξ̄ = sum(real, view(ξ, group)) / g
    V = svd(hcat((Vector(real.(U[q])) for q in group)..., (Vector(imag.(U[q])) for q in group)...)).U[:, 1:g]
    W = V * eigen(Symmetric(transpose(V) * Matrix(S + 2 * ξ̄ * Tm) * V)).vectors
    for (j, q) in enumerate(group)
        ξ[q] = ξ̄
        U[q] = SVector{3}(complex.(view(W, :, j)))
    end
    return nothing
end

# Dual numbers: the waves of the values, then the derivative of each root and
# of its polarization from the pencil P(ξ) = Q̃ + ξ(R + Rᵀ) + ξ²T, which is
# symmetric:  P·U̇ + ξ̇ P′·U = −Ṗ·U.  Multiplied on the left by U, whose
# transpose P cancels, it gives  ξ̇ = −(U·Ṗ·U)/(U·P′·U) ; U̇ then solves the
# system with Ū·U̇ = 0 (unit polarization). One nesting level at a time, as
# `_sym_eigen` does, and with its convention for the roots that coincide (the
# two shear waves of an isotropic medium): within their plane no polarization
# is preferred, so U̇ is taken orthogonal to all of them, by least squares.
function _interface_waves(C, ρ::D, ν::SVector{3, D}, sₜ::SVector{3, D}) where {D <: ForwardDiff.Dual}
    val = ForwardDiff.value
    forward₀, backward₀ = _interface_waves(val.(Array(C)), val(ρ), val.(ν), val.(sₜ))
    Q, R, Tm = _interface_blocks(C, sₜ, ν)
    Q̃ = Q - ρ * one(Q)
    S = R + transpose(R)
    tol = 1.0e-8 * sqrt(primal(ρ) / maximum(abs ∘ primal, Tm))            # against a slowness of the medium
    all₀ = (forward₀..., backward₀...)
    function lift(w)
        group = filter(v -> abs(v.normal_slowness - w.normal_slowness) ≤ tol, all₀)
        return _lift_interface_wave(w, map(v -> v.polarization, group), Q̃, S, R, Tm, sₜ, ν, D)
    end
    return map(lift, forward₀), map(lift, backward₀)
end

function _lift_interface_wave(w::InterfaceWave{V}, group, Q̃, S, R, Tm, sₜ, ν, ::Type{D}) where {V, D <: ForwardDiff.Dual}
    val = ForwardDiff.value
    ξ₀, U₀ = w.normal_slowness, w.polarization
    P = val.(Q̃) + ξ₀ * val.(S) + ξ₀^2 * val.(Tm)
    P′U = (val.(S) + 2 * ξ₀ * val.(Tm)) * U₀
    A = vcat(Matrix(P), map(u -> adjoint(Vector(u)), group)...)
    ẋ = ntuple(Val(ForwardDiff.npartials(D))) do k
        ∂(M) = ForwardDiff.partials.(M, k)
        ṖU = (∂(Q̃) + ξ₀ * ∂(S) + ξ₀^2 * ∂(Tm)) * U₀
        ξ̇ = -(transpose(U₀) * ṖU) / (transpose(U₀) * P′U)
        U̇ = A \ vcat(-Vector(ṖU + ξ̇ * P′U), zeros(Complex{V}, length(group)))
        (U̇[1], U̇[2], U̇[3], ξ̇)
    end
    dual(z₀, ż) = complex(D(real(z₀), ForwardDiff.Partials(map(real, ż))), D(imag(z₀), ForwardDiff.Partials(map(imag, ż))))
    U = SVector{3}(ntuple(a -> dual(U₀[a], map(x -> x[a], ẋ)), Val(3)))
    ξ = dual(ξ₀, map(x -> x[4], ẋ))
    return _interface_wave(w.homogeneous ? complex(real(ξ)) : ξ, U, R, Tm, sₜ, ν, w.homogeneous)
end

"""
    reflection_transmission(ℂ₁, ρ₁, ℂ₂, ρ₂, ν, incident::InterfaceWave; slip = nothing)
        → (reflected, transmitted, R, T, energy)

The waves that the `incident` wave, one of the forward waves of the medium 1
([`interface_waves`](@ref), on the side ``\\underline{\\nu}\\cdot\\underline{x} < 0``),
excites at the plane interface with the medium 2: the three `reflected` waves
of the medium 1 and the three `transmitted` waves of the medium 2, with their
complex amplitudes `R` and `T` per unit incident amplitude (displacement
amplitudes, unit polarizations) and `energy`, the fractions of the incident
energy flux carried by each of the six waves (zero for an evanescent wave):
they sum to one.

The interface is welded, displacement and traction continuous. With `slip`, a
3×3 matrix ``\\omega\\boldsymbol{Z}`` (angular frequency times a compliance of
the interface), the displacement jumps by ``\\boldsymbol{Z}`` times the
traction (linear slip); the traction stays continuous
([Waves at an interface](@ref th-interface-waves)). Numbers and dual numbers.
"""
function reflection_transmission(ℂ₁::AbstractTens{4, 3}, ρ₁, ℂ₂::AbstractTens{4, 3}, ρ₂, ν, incident::InterfaceWave{T}; slip = nothing) where {T}
    sₜ = real.(incident.slowness)
    _, reflected = interface_waves(ℂ₁, ρ₁, ν, sₜ)
    transmitted, _ = interface_waves(ℂ₂, ρ₂, ν, sₜ)
    Z = _slip_matrix(slip, T)
    # u₂ − u₁ = iωZ·τ and τ₂ = τ₁, with the common factor iω of the tractions kept in Z
    col(w, side) = vcat(side * w.polarization - _slip_term(Z, w, side), side * w.traction)
    A = hcat(map(w -> col(w, -1), reflected)..., map(w -> col(w, 1), transmitted)...)
    b = vcat(incident.polarization, incident.traction)
    x = Matrix(A) \ Vector(b)
    R, Tt = SVector{3}(x[1:3]), SVector{3}(x[4:6])
    energy = (
        map((a, w) -> -abs2(a) * w.flux / incident.flux, R, SVector(reflected)),
        map((a, w) -> abs2(a) * w.flux / incident.flux, Tt, SVector(transmitted)),
    )
    return (reflected = reflected, transmitted = transmitted, R = R, T = Tt, energy = energy)
end

"""
    free_surface_reflection(ℂ::AbstractTens{4,3}, ρ, ν, incident::InterfaceWave) → (reflected, R, energy)

The waves that the `incident` wave, one of the forward waves of the medium
([`interface_waves`](@ref), on the side ``\\underline{\\nu}\\cdot\\underline{x} < 0``),
excites at the free surface ``\\underline{\\nu}\\cdot\\underline{x} = 0``: the three
`reflected` waves, their complex amplitudes `R` per unit incident amplitude
(unit polarizations) and `energy`, the fractions of the incident energy flux
they carry, which sum to one. The traction vanishes on the surface,

```math
\\underline{\\tau}^{\\mathrm{inc}} + \\sum_{m=1}^{3} R_m\\,\\underline{\\tau}_m^{\\mathrm{r}} = \\underline{0}
```

([Waves at an interface](@ref th-interface-waves)). Numbers and dual numbers.
"""
function free_surface_reflection(ℂ::AbstractTens{4, 3}, ρ, ν, incident::InterfaceWave)
    _, reflected = interface_waves(ℂ, ρ, ν, real.(incident.slowness))
    A = hcat(map(w -> w.traction, reflected)...)
    R = SVector{3}(Matrix(A) \ Vector(-incident.traction))
    energy = map((a, w) -> -abs2(a) * w.flux / incident.flux, R, SVector(reflected))
    return (reflected = reflected, R = R, energy = energy)
end

"""
    weak_contrast_reflection(ℂ, ρ, δℂ, δρ, ν, incident::InterfaceWave) → (reflected, R)

The reflection coefficients, to first order in the contrast, at the welded
interface between the medium ``(\\mathbb{C}, \\rho)`` of the `incident` wave and
the medium ``(\\mathbb{C} + \\delta\\mathbb{C}, \\rho + \\delta\\rho)``: with the
waves of the first medium only, ``i`` incident and ``n`` reflected,

```math
R_n = \\frac{\\langle\\underline{\\eta}_n, \\delta\\underline{\\eta}_i\\rangle}{\\langle\\underline{\\eta}_n, \\underline{\\eta}_n\\rangle},
\\qquad
\\langle\\underline{\\eta}_a, \\underline{\\eta}_b\\rangle = \\underline{U}_a\\cdot\\underline{\\tau}_b + \\underline{\\tau}_a\\cdot\\underline{U}_b ,
```

``\\delta\\underline{\\eta}_i`` being the change of the displacement and traction
of the incident wave when its medium becomes the second one, written without
its polarization change ([Waves at an interface](@ref th-weak-contrast)). Valid
for any anisotropy, the incident wave a simple root. Numbers, dual numbers,
and symbols where the waves are symbolic.
"""
function weak_contrast_reflection(ℂ::AbstractTens{4, 3}, ρ, δℂ::AbstractTens{4, 3}, δρ, ν, incident::InterfaceWave)
    sₜ = real.(incident.slowness)
    _, reflected = interface_waves(ℂ, ρ, ν, sₜ)
    n = _unit_direction(ν, real(eltype(incident.slowness)))
    s = sₜ - (transpose(sₜ) * n) * n
    Q, R, Tm = _interface_blocks(canonical_components(ℂ), s, n)
    δQ, δR, δT = _interface_blocks(canonical_components(δℂ), s, n)
    ξᵢ, Uᵢ = incident.normal_slowness, incident.polarization
    P′(ξ) = R + transpose(R) + 2 * ξ * Tm
    δP = δQ - δρ * one(δQ) + ξᵢ * (δR + transpose(δR)) + ξᵢ^2 * δT
    δξ = -(transpose(Uᵢ) * δP * Uᵢ) / (transpose(Uᵢ) * P′(ξᵢ) * Uᵢ)
    coefficient(w) = (
        transpose(w.polarization) * ((transpose(δR) + ξᵢ * δT + δξ * Tm) * Uᵢ) -
            transpose(w.polarization) * ((δξ * P′(ξᵢ) + δP) * Uᵢ) / (ξᵢ - w.normal_slowness)
    ) / (transpose(w.polarization) * P′(w.normal_slowness) * w.polarization)
    return (reflected = reflected, R = map(coefficient, SVector(reflected)))
end

_slip_matrix(::Nothing, ::Type{T}) where {T} = nothing
_slip_matrix(Z::AbstractMatrix, ::Type{T}) where {T} = SMatrix{3, 3, Complex{T}}(Z)
# The jump is carried by the transmitted side: u₂ = u₁ + iωZ·τ₂.
_slip_term(::Nothing, w, side) = zero(w.polarization)
_slip_term(Z::SMatrix, w, side) = side == 1 ? im * (Z * w.traction) : zero(w.polarization)

export InterfaceWave, interface_waves, reflection_transmission, free_surface_reflection, weak_contrast_reflection
