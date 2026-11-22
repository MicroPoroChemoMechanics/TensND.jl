# ── Distances between elasticity tensors ───────────────────────────────────────
#
# A stiffness ℂ and its compliance 𝕊 = ℂ⁻¹ describe the same material, and the
# Euclidean (Frobenius) distance tells them apart: ‖ℂ₁ − ℂ₂‖ and ‖𝕊₁ − 𝕊₂‖ have
# neither the same value nor the same unit, and the tensor of a symmetry class
# closest to ℂ is not the inverse of the one closest to 𝕊. The distances of this
# file are built on a function of the tensor, f(ℂ), taken on its eigenvalues:
#
#   - `EuclideanDistance`:      ‖ℂ₁ − ℂ₂‖, the distance of `project`;
#   - `LogEuclideanDistance`:   ‖ln ℂ₁ − ln ℂ₂‖;
#   - `FunctionEuclideanDistance(f, f⁻¹, E₀)`: ‖f(ℂ₁/E₀) − f(ℂ₂/E₀)‖, with a gauge
#     E₀ that makes the tensors dimensionless (`ArctanEuclideanDistance`,
#     `PowerEuclideanDistance`);
#   - `RiemannianDistance`:     ‖ln(ℂ₁^{-1/2} ℂ₂ ℂ₁^{-1/2})‖.
#
# All but the first are invariant under inversion when f(x) + f(1/x) is a
# constant. A symmetric fourth-order tensor (minor and major symmetries) is a
# symmetric matrix in the Kelvin-Mandel notation, whose eigenvalues and
# projectors define f(ℂ); a symmetric second-order tensor is its own matrix.

# ── A function of a symmetric tensor ───────────────────────────────────────────

# Relative closeness under which two eigenvalues are one for the divided
# difference of a function: below it (f(λᵢ) − f(λⱼ))/(λᵢ − λⱼ) loses its digits,
# and the derivative at their mean differs from it by a term of second order.
const _DIVDIFF_RTOL = 1.0e-6

# f of a real symmetric matrix, V f(Λ) Vᵀ.
function _matrix_function(f::F, M::AbstractMatrix{<:AbstractFloat}) where {F}
    λ, V = _sym_eigen(M)
    return V * Diagonal(f.(λ)) * transpose(V)
end

# With dual entries, one level at a time: the value is that of the value
# matrix, and each derivative is given by the formula of Daleckii and Krein,
#   Df(M)[Ṁ] = V (Δ ∘ (Vᵀ Ṁ V)) Vᵀ,   Δᵢⱼ = (f(λᵢ) − f(λⱼ))/(λᵢ − λⱼ),  Δᵢᵢ = f′(λᵢ),
# which, unlike the derivatives of the eigenvectors, has a limit when two
# eigenvalues meet.
function _matrix_function(f::F, M::AbstractMatrix{D}) where {F, D <: ForwardDiff.Dual}
    M₀ = ForwardDiff.value.(M)
    λ, V = _sym_eigen(M₀)
    fλ = f.(λ)
    f′λ = map(x -> ForwardDiff.derivative(f, x), λ)
    n = length(λ)
    scale = max(abs(λ[n]), abs(λ[1]))
    Δ = [
        abs(λ[i] - λ[j]) ≤ _DIVDIFF_RTOL * scale ? (f′λ[i] + f′λ[j]) / 2 : (fλ[i] - fλ[j]) / (λ[i] - λ[j])
            for i in 1:n, j in 1:n
    ]
    Vᵀ = transpose(V)
    F₀ = V * Diagonal(fλ) * Vᵀ
    Ḟ = ntuple(k -> V * (Δ .* (Vᵀ * ForwardDiff.partials.(M, k) * V)) * Vᵀ, Val(ForwardDiff.npartials(D)))
    return [D(F₀[a, b], ForwardDiff.Partials(map(Ḟₖ -> Ḟₖ[a, b], Ḟ))) for a in 1:n, b in 1:n]
end

# Symbols: no eigendecomposition. The structured classes below are functions of
# their coefficients and need none.
_matrix_function(f, M::AbstractMatrix) = throw(
    ArgumentError(
        "the function of a general tensor is taken on its eigenvalues, which needs numbers " *
            "(or duals); with element type $(eltype(M)), use a structured tensor (`TensISO`, " *
            "`TensCubic`, `TensTI`), whose function is that of its coefficients"
    )
)

# The symmetric matrix that stands for a symmetric tensor: the tensor itself at
# order two, its Kelvin-Mandel matrix at order four. Refused when it is not
# symmetric, rather than symmetrized.
const _SYMMETRY_RTOL = 1.0e-10
function _symmetric_matrix(t::AbstractTens{order, dim}) where {order, dim}
    order == 2 || order == 4 || throw(
        ArgumentError("the function of a tensor is defined at orders 2 and 4; got order $order")
    )
    M = Matrix(order == 2 ? canonical_components(t) : mandel(t, CanonicalBasis{dim, eltype(t)}()))
    Mᵀ = transpose(M)
    _is_symmetric_matrix(Val(is_hard_numeric(eltype(M))), M, Mᵀ) || throw(
        ArgumentError(
            "the function of a tensor is defined for a symmetric one (with the major symmetry " *
                "at order 4); take its symmetric part first if that is what is meant"
        )
    )
    return (M + Mᵀ) / 2
end
_is_symmetric_matrix(::Val{true}, M, Mᵀ) = _frobenius(M - Mᵀ) ≤ _SYMMETRY_RTOL * _frobenius(M)
_is_symmetric_matrix(::Val{false}, M, Mᵀ) = all(iszero ∘ _simplified, M - Mᵀ)

_tensor_of_matrix(F, ::AbstractTens{2}) = Tens(F)
_tensor_of_matrix(F, ::AbstractTens{4}) = from_mandel(F)

"""
    tensor_function(f, t::AbstractTens) → AbstractTens

The function `f` of the symmetric tensor `t`: the tensor that has the
eigenprojectors of `t` and the eigenvalues ``f(\\lambda_i)``,

```math
f(\\mathbb{C}) = \\sum_i f(\\lambda_i)\\,\\mathbb{P}_i
\\qquad\\text{for}\\qquad \\mathbb{C} = \\sum_i \\lambda_i\\,\\mathbb{P}_i .
```

At order 4 (with the minor and the major symmetries) the eigenvalues are those
of the Kelvin-Mandel matrix; at order 2 those of the tensor. `log`, `exp`,
`sqrt`, `inv`, `atan` or any scalar function may be given: `tensor_function(log, ℂ)`
is the logarithm of a positive-definite stiffness.

A structured tensor is a combination of orthogonal projectors, so that `f` is
taken on its coefficients and the class is kept, with symbolic coefficients
too: ``f(\\alpha\\,\\mathbb{J} + \\beta\\,\\mathbb{K}) = f(\\alpha)\\,\\mathbb{J} + f(\\beta)\\,\\mathbb{K}``
for a `TensISO`, the three coefficients of a `TensCubic`, the block of the two
coupled projectors and the two shear coefficients of a `TensTI`. A general
tensor goes through the eigendecomposition of its matrix, with numbers or
`ForwardDiff` duals, and comes back as a general tensor in the canonical basis.

# Examples
```jldoctest
julia> ℂ = iso_stiffness(10.0, 4.0);                 # 3k = 30, 2μ = 8

julia> get_data(tensor_function(log, ℂ)) == (log(30.0), log(8.0))
true

julia> get_array(tensor_function(sqrt, Tens(get_array(ℂ))) ⊡ tensor_function(sqrt, Tens(get_array(ℂ)))) ≈ get_array(ℂ)
true
```
"""
tensor_function(f::F, t::AbstractTens) where {F} = _tensor_of_matrix(_matrix_function(f, _symmetric_matrix(t)), t)

tensor_function(f::F, t::TensISO{order, dim}) where {F, order, dim} = TensISO{dim}(map(f, get_data(t))...)
tensor_function(f::F, t::TensCubic) where {F} = TensCubic(map(f, get_data(t))..., frame(t))

# A major-symmetric TI tensor: the projectors of ℓ₅ and ℓ₆ are orthogonal to the
# others, and (ℓ₁, ℓ₂, ℓ₃) is the symmetric 2×2 matrix [ℓ₁ ℓ₃; ℓ₃ ℓ₂] on the two
# coupled ones. With its eigenvalues λ±, f of that matrix is
#   (f(λ₊) + f(λ₋))/2 𝟏 + (f(λ₊) − f(λ₋))/(λ₊ − λ₋) (M − tr(M)/2 𝟏),
# whose second term has the limit f′ when the eigenvalues meet: a comparison,
# made on numbers; a symbolic matrix is taken with distinct eigenvalues.
function tensor_function(f::F, t::TensTI{4, T, 5}) where {F, T}
    ℓ₁, ℓ₂, ℓ₃, _, ℓ₅, ℓ₆ = walpole_coeffs(t)
    m, δ = (ℓ₁ + ℓ₂) / 2, (ℓ₁ - ℓ₂) / 2
    r = sqrt(δ^2 + ℓ₃^2)
    f₊, f₋ = f(m + r), f(m - r)
    slope = _block_slope(Val(is_hard_numeric(typeof(r))), f, f₊, f₋, m, r)
    return TensTI{4}((f₊ + f₋) / 2 + slope * δ, (f₊ + f₋) / 2 - slope * δ, slope * ℓ₃, f(ℓ₅), f(ℓ₆), axis(t))
end
_block_slope(::Val{false}, f, f₊, f₋, m, r) = (f₊ - f₋) / (2 * r)
_block_slope(::Val{true}, f::F, f₊, f₋, m, r) where {F} =
    abs(r) ≤ _DIVDIFF_RTOL * abs(m) ? ForwardDiff.derivative(f, m) : (f₊ - f₋) / (2 * r)

# ── The distances ─────────────────────────────────────────────────────────────

"""
    AbstractTensorDistance

Supertype of the distances between two symmetric tensors of the same order and
dimension, typically two stiffnesses: [`EuclideanDistance`](@ref),
[`LogEuclideanDistance`](@ref), [`RiemannianDistance`](@ref),
[`FunctionEuclideanDistance`](@ref). A distance `d` is evaluated by
[`tensor_distance`](@ref)`(A, B, d)`, and [`closest_tensor`](@ref) projects on
a symmetry class in its sense.
"""
abstract type AbstractTensorDistance end

"""
    EuclideanDistance()

The Euclidean distance ``\\lVert\\mathbb{A} - \\mathbb{B}\\rVert``, with
``\\lVert\\mathbb{A}\\rVert^2 = A_{ijkl}A_{ijkl}``: that of [`project`](@ref). It
is not invariant under inversion: the distance between two stiffnesses is not
that between their compliances.
"""
struct EuclideanDistance <: AbstractTensorDistance end

"""
    LogEuclideanDistance()

The log-Euclidean distance
``\\lVert\\ln\\mathbb{A} - \\ln\\mathbb{B}\\rVert`` between two positive-definite
tensors [arsigny2006](@cite). Since ``\\ln\\mathbb{A}^{-1} = -\\ln\\mathbb{A}``, it
is invariant under inversion, and it does not depend on the unit of the
tensors: multiplying both by a number changes nothing.
"""
struct LogEuclideanDistance <: AbstractTensorDistance end

"""
    RiemannianDistance()

The Riemannian (affine-invariant) distance between two positive-definite
tensors [moakher2006spd, moakher2006](@cite),

```math
d(\\mathbb{A}, \\mathbb{B}) = \\bigl\\lVert\\ln\\bigl(\\mathbb{A}^{-1/2} : \\mathbb{B} : \\mathbb{A}^{-1/2}\\bigr)\\bigr\\rVert
= \\Bigl(\\sum_i \\ln^2\\mu_i\\Bigr)^{1/2} ,
```

with ``\\mu_i`` the eigenvalues of ``\\mathbb{A}^{-1} : \\mathbb{B}``. It is
invariant under inversion and under any change ``\\mathbb{A} \\mapsto
\\mathbb{G} : \\mathbb{A} : \\mathbb{G}^{\\mathsf T}`` applied to both tensors.
"""
struct RiemannianDistance <: AbstractTensorDistance end

"""
    FunctionEuclideanDistance(f, f⁻¹, gauge)

The distance ``\\lVert f(\\mathbb{A}/E_0) - f(\\mathbb{B}/E_0)\\rVert`` built on a
strictly monotone function `f` of the positive reals, with the gauge
``E_0`` = `gauge`, a modulus that makes the tensors dimensionless
[morin2020](@cite). It is invariant under inversion when ``f(x) + f(1/x)`` is a
constant. `f⁻¹` is the inverse of `f`, which [`closest_tensor`](@ref) needs.

[`ArctanEuclideanDistance`](@ref) and [`PowerEuclideanDistance`](@ref) are two
such distances; [`LogEuclideanDistance`](@ref) is the one of `f = log`, which
no gauge enters.
"""
struct FunctionEuclideanDistance{F, G, T <: Number} <: AbstractTensorDistance
    f::F
    inverse::G
    gauge::T
end

"""
    ArctanEuclideanDistance(gauge)

The [`FunctionEuclideanDistance`](@ref) of ``f = \\arctan``, invariant under
inversion since ``\\arctan x + \\arctan(1/x) = \\pi/2`` [morin2020](@cite). For
two moduli it is the angle between the two straight lines of slopes
``E_1/E_0`` and ``E_2/E_0``.
"""
ArctanEuclideanDistance(gauge::Number) = FunctionEuclideanDistance(atan, tan, gauge)

# The function x ↦ xⁿ − x⁻ⁿ of the power-Euclidean distance and its inverse,
# y ↦ ((y + √(y² + 4))/2)^(1/n) for n > 0: xⁿ is the positive root of
# u² − y u − 1 = 0.
struct _PowerDifference{N <: Real}
    n::N
end
(p::_PowerDifference)(x) = x^p.n - x^(-p.n)
struct _InversePowerDifference{N <: Real}
    n::N
end
(q::_InversePowerDifference)(y) = ((y + sqrt(y^2 + 4)) / 2)^inv(q.n)

"""
    PowerEuclideanDistance(n, gauge)

The [`FunctionEuclideanDistance`](@ref) of ``f(x) = x^n - x^{-n}``, ``n > 0``,
for which ``f(x) + f(1/x) = 0`` [morin2020](@cite).
"""
function PowerEuclideanDistance(n::Real, gauge::Number)
    n > 0 || throw(ArgumentError("the exponent of a power-Euclidean distance is positive; got $n"))
    return FunctionEuclideanDistance(_PowerDifference(n), _InversePowerDifference(n), gauge)
end

# Frobenius norm of the difference of two tensors, on their canonical components.
_difference_norm(A::AbstractTens, B::AbstractTens) = _frobenius(canonical_components(A) - canonical_components(B))
# Two isotropic tensors: 𝕁 ∷ 𝕁 = 1 and 𝕂 ∷ 𝕂 = d(d+1)/2 − 1, with no array (and with symbols).
function _difference_norm(A::TensISO{4, dim}, B::TensISO{4, dim}) where {dim}
    (α₁, β₁), (α₂, β₂) = get_data(A), get_data(B)
    return sqrt((α₁ - α₂)^2 + (dim * (dim + 1) ÷ 2 - 1) * (β₁ - β₂)^2)
end

# The tensor a distance compares: f(A/E₀).
_transformed(::EuclideanDistance, A) = A
_transformed(::LogEuclideanDistance, A) = tensor_function(log, A)
_transformed(d::FunctionEuclideanDistance, A) = tensor_function(d.f, A / d.gauge)

"""
    tensor_distance(A, B, d::AbstractTensorDistance = EuclideanDistance())

The distance `d` between the symmetric tensors `A` and `B`, of the same order
and dimension. For every distance but the Euclidean one the tensors are
positive-definite.

# Examples
The log-Euclidean distance between two stiffnesses is that between their
compliances; the Euclidean one is not.
```jldoctest
julia> ℂ₁, ℂ₂ = iso_stiffness(10.0, 4.0), iso_stiffness(14.0, 3.0);

julia> tensor_distance(ℂ₁, ℂ₂, LogEuclideanDistance()) ≈ tensor_distance(inv(ℂ₁), inv(ℂ₂), LogEuclideanDistance())
true

julia> tensor_distance(ℂ₁, ℂ₂) ≈ tensor_distance(inv(ℂ₁), inv(ℂ₂))
false
```

See also [`closest_tensor`](@ref), [`tensor_function`](@ref).
"""
tensor_distance(A::AbstractTens, B::AbstractTens, d::AbstractTensorDistance = EuclideanDistance()) =
    _difference_norm(_transformed(d, A), _transformed(d, B))

# ln(A^{-1/2} B A^{-1/2}): its eigenvalues are the logarithms of those of A⁻¹B.
function _riemannian_logarithm(A::AbstractTens, B::AbstractTens)
    L = _matrix_function(x -> inv(sqrt(x)), _symmetric_matrix(A))
    M = L * _symmetric_matrix(B) * L
    return _matrix_function(log, (M + transpose(M)) / 2)
end
tensor_distance(A::AbstractTens, B::AbstractTens, ::RiemannianDistance) = _frobenius(_riemannian_logarithm(A, B))
# Two isotropic tensors commute: ln(A⁻¹B) = ln(α₂/α₁) 𝕁 + ln(β₂/β₁) 𝕂.
function tensor_distance(A::TensISO{4, dim}, B::TensISO{4, dim}, ::RiemannianDistance) where {dim}
    (α₁, β₁), (α₂, β₂) = get_data(A), get_data(B)
    return sqrt(log(α₂ / α₁)^2 + (dim * (dim + 1) ÷ 2 - 1) * log(β₂ / β₁)^2)
end

# ── The closest tensor of a symmetry class ────────────────────────────────────

# A symmetry class with its axis or its frame, as `project` takes them:
# `(Val(:ISO),)`, `(Val(:TI), axis)`, `(Val(:ORTHO), frame)`, `(Val(:CUBIC), frame)`.
_euclidean_projection(A, class...) = first(project(A, class...))
# The isotropic projection with no comparison, so that symbols go through.
_euclidean_projection(A, ::Val{:ISO}) = isotropify(A)

"""
    closest_tensor(A, class, [axis or frame]; distance = EuclideanDistance()) → (B, d)

The tensor `B` of the symmetry class `class` closest to the symmetric tensor
`A` in the sense of `distance`, and that distance `d`. The class is `:ISO`,
`:TI` with its axis, `:ORTHO` or `:CUBIC` with its frame, as in
[`project`](@ref), whose result this is for the Euclidean distance.

For a distance built on a function ``f``
([`LogEuclideanDistance`](@ref), [`FunctionEuclideanDistance`](@ref)), a
symmetry class with a given axis or frame is a linear space that ``f`` maps
into itself, so that the closest tensor is

```math
\\mathbb{B} = E_0\\,f^{-1}\\Bigl(\\Pi\\bigl(f(\\mathbb{A}/E_0)\\bigr)\\Bigr) ,
```

with ``\\Pi`` the Euclidean projection on the class. For the closest isotropic
tensor this gives ``f(3k/E_0) = \\mathbb{J} :: f(\\mathbb{A}/E_0)`` and
``f(2\\mu/E_0) = \\mathbb{K} :: f(\\mathbb{A}/E_0)/5`` [morin2020](@cite). The
closest tensor in the sense of the [`RiemannianDistance`](@ref) has no such
form: it is computed by the method of Gauss and Newton on the logarithm of
`B`, from the log-Euclidean one, for numbers only.

With a distance invariant under inversion, the tensor closest to a stiffness
is the inverse of the one closest to its compliance.

# Examples
```jldoctest
julia> ℂ = TensCubic(30.0, 8.0, 12.0, CanonicalBasis{3, Float64}());

julia> B, d = closest_tensor(ℂ, :ISO; distance = LogEuclideanDistance());

julia> get_data(B)[2] ≈ exp((2 * log(8.0) + 3 * log(12.0)) / 5)
true

julia> S, _ = closest_tensor(inv(ℂ), :ISO; distance = LogEuclideanDistance());

julia> all(get_data(inv(S)) .≈ get_data(B))
true
```

See also [`tensor_distance`](@ref), [`project`](@ref), [`best_fit_iso`](@ref).
"""
closest_tensor(A::AbstractTens, class::Symbol, args...; distance::AbstractTensorDistance = EuclideanDistance()) =
    _closest_tensor(distance, A, Val(class), args...)

function _closest_tensor(d::AbstractTensorDistance, A, class...)
    B = _closest(d, A, class...)
    return B, tensor_distance(A, B, d)
end

# The function of a tensor of the class is in the class; projecting it once
# more gives it the structured type of the class, whatever the route it took.
_closest(::EuclideanDistance, A, class...) = _euclidean_projection(A, class...)
_closest(::LogEuclideanDistance, A, class...) =
    _euclidean_projection(tensor_function(exp, _euclidean_projection(tensor_function(log, A), class...)), class...)
_closest(d::FunctionEuclideanDistance, A, class...) = _euclidean_projection(
    d.gauge * tensor_function(d.inverse, _euclidean_projection(tensor_function(d.f, A / d.gauge), class...)), class...
)

# The Riemannian case. The class is parameterized by the logarithm X of B, a
# point of the linear space of the class, written on the vector of its
# coefficients (`to_vector`), so that B = exp(X) is positive-definite whatever
# the coefficients. The squared distance is the sum of the squares of the
# entries of ln(A^{-1/2} exp(X) A^{-1/2}): a least-squares problem, solved by
# the method of Gauss and Newton with a step halved until the distance
# decreases. The Jacobian is exact, by the formula of Daleckii and Krein.
const _RIEMANN_MAXITER = 100
const _RIEMANN_RTOL = 1.0e-13

function _closest(::RiemannianDistance, A, class...)
    is_hard_numeric(eltype(A)) && !(eltype(A) <: ForwardDiff.Dual) || throw(
        ArgumentError(
            "the tensor closest in the Riemannian sense is computed by an iteration, for " *
                "floating-point numbers; got the element type $(eltype(A)). The log-Euclidean " *
                "distance has a closed form, with duals and structured symbolic tensors"
        )
    )
    X = _euclidean_projection(tensor_function(log, A), class...)
    L = _matrix_function(x -> inv(sqrt(x)), _symmetric_matrix(A))
    exponential(p) = tensor_function(exp, from_vector(X, p))
    residual(p) = vec(_matrix_function(log, _congruence(L, _symmetric_matrix(exponential(p)))))
    p = collect(float.(to_vector(X)))
    r = residual(p)
    for _ in 1:_RIEMANN_MAXITER
        J = ForwardDiff.jacobian(residual, p)
        δ = -(J \ r)
        step, r₊ = _halved_step(residual, p, δ, r)
        p += step
        converged = sum(abs2, step) ≤ _RIEMANN_RTOL^2 * max(sum(abs2, p), one(eltype(p)))
        r = r₊
        converged && break
    end
    return _euclidean_projection(exponential(p), class...)
end
_congruence(L, M) = (N = L * M * L; (N + transpose(N)) / 2)
# The Gauss-Newton step, halved until the residual decreases (or no longer can).
function _halved_step(residual, p, δ, r)
    t, r₀ = one(eltype(p)), sum(abs2, r)
    r₊ = residual(p + δ)
    while sum(abs2, r₊) > r₀ && t > eps(eltype(p))
        t /= 2
        r₊ = residual(p + t * δ)
    end
    return t * δ, r₊
end

export AbstractTensorDistance, EuclideanDistance, LogEuclideanDistance, RiemannianDistance
export FunctionEuclideanDistance, ArctanEuclideanDistance, PowerEuclideanDistance
export tensor_function, tensor_distance, closest_tensor
