# ─────────────────────────────────────────────────────────────────────────────
#  tens_plane_projection.jl — projections of a tensor of the plane onto its
#  symmetry classes, and its average over the rotations of the plane.
#
#  The classes of the plane are the isotropic one, the square one and the
#  orthotropic one (`tens_plane_classes.jl`); the plane has no transversely
#  isotropic class. With a frame given, the projections are those of space: the
#  Kelvin-Mandel matrix written in the frame, then the entries the class does
#  not have set to zero and those it identifies averaged.
#
#  With the frame free, the plane needs no search. Under a rotation of the frame
#  by θ a major- and minor-symmetric tensor of order 4 splits into two
#  invariants and two harmonics, of orders 2 and 4 in θ:
#
#      d₁₂ = C₁₁₁₂ + C₂₂₁₂ ,  δ = (C₁₁₁₁ − C₂₂₂₂)/2 ,
#      a = C₁₁₁₂ − C₂₂₁₂ ,    b = (C₁₁₁₁ + C₂₂₂₂ − 2C₁₁₂₂ − 4C₁₂₁₂)/4 ,
#
#      d₁₂′ = d₁₂ cos 2θ − δ sin 2θ ,    δ′ = δ cos 2θ + d₁₂ sin 2θ ,
#      a′  = a cos 4θ − b sin 4θ ,       b′ = b cos 4θ + a sin 4θ .
#
#  A tensor is orthotropic in the frame when d₁₂′ = a′ = 0 and square when
#  moreover δ′ = 0, so that the squared distances to the two classes in the
#  frame at θ are
#
#      square:       2 (d₁₂² + δ²) + 2 a′² ,
#      orthotropic:  2 (d₁₂′² + a′²) .
#
#  The first is least where a′ = 0, θ = atan(a, b)/4, in closed form. The
#  second is a trigonometric polynomial of z = e^{4iθ},
#  const + Re(c₁ z + c₂ z²) with c₁ = (d₁₂ + iδ)², c₂ = (a + ib)², whose
#  stationary points are the roots of modulus one of
#
#      2c₂ z⁴ + c₁ z³ − c̄₁ z − 2c̄₂ = 0 :
#
#  the least of them is taken, an exact global minimum.
# ─────────────────────────────────────────────────────────────────────────────

# ── With the frame given ─────────────────────────────────────────────────────

# The Kelvin-Mandel matrix of the minor-symmetric part of `A` in the frame.
function _plane_km_in_frame(A::AbstractArray{T, 4}, frame::OrthonormalBasis{2}) where {T}
    Q = mandel_rotation(_frame_matrix(frame))
    # The element type is said: what holds the array is only known at run time.
    return transpose(Q) * SMatrix{3, 3, T}(_KM_of_array(A)) * Q
end

"""
    project(A::AbstractArray{T,4}, ::Val{:ORTHO}, frame::OrthonormalBasis{2}) → (TensOrtho{4,2}, d, drel)
    project(A::AbstractArray{T,4}, ::Val{:CUBIC}, frame::OrthonormalBasis{2}) → (TensCubic{2}, d, drel)
    project(A::AbstractArray{T,2}, ::Val{:ORTHO}, frame::OrthonormalBasis{2}) → (TensOrtho{2,2}, d, drel)

The closest orthotropic or square tensor **of the plane** in the frame `frame`
(`Basis(θ)`, `CanonicalBasis{2, T}()`), with its distance `d` to `A` and the
relative one `drel`: in the frame, the entries `C₁₁₁₂`, `C₂₂₁₂` of the
Kelvin-Mandel matrix are dropped, and `C₁₁₁₁`, `C₂₂₂₂` averaged for the square
class. The orthotropic projection is major-symmetric (four constants). At
order 2 the orthotropic projection is the diagonal in the frame, and the
square one the isotropic one.

Without a frame, `project(t, :ORTHO)` and `project(t, :CUBIC)` find the best one
in closed form: see [Projections in the plane](@ref th-projection-plane).
"""
function project(A::AbstractArray{T, 4}, ::Val{:ORTHO}, frame::OrthonormalBasis{2}) where {T}
    nA = _frobenius(A)
    if nA ≈ zero(T)
        z = zero(T)
        return TensOrtho(z, z, z, z, frame), z, z
    end
    M = _plane_km_in_frame(A, frame)
    B = TensOrtho(M[1, 1], M[2, 2], (M[1, 2] + M[2, 1]) / 2, M[3, 3] / 2, frame)
    d = _frobenius(get_array(B) - A)
    return B, d, d / nA
end
function project(A::AbstractArray{T, 4}, ::Val{:CUBIC}, frame::OrthonormalBasis{2}) where {T}
    nA = _frobenius(A)
    if nA ≈ zero(T)
        z = zero(T)
        return TensCubic(z, z, z, frame), z, z
    end
    M = _plane_km_in_frame(A, frame)
    B = cubic_from_constants((M[1, 1] + M[2, 2]) / 2, (M[1, 2] + M[2, 1]) / 2, M[3, 3] / 2, frame)
    d = _frobenius(get_array(B) - A)
    return B, d, d / nA
end
function project(A::AbstractArray{T, 2}, ::Val{:ORTHO}, frame::OrthonormalBasis{2}) where {T}
    nA = _frobenius(A)
    if nA ≈ zero(T)
        z = zero(T)
        return TensOrtho{2}(z, z, frame), z, z
    end
    R = _frame_matrix(frame)
    M = transpose(R) * SMatrix{2, 2}(A) * R
    B = TensOrtho{2}(M[1, 1], M[2, 2], frame)
    d = _frobenius(get_array(B) - A)
    return B, d, d / nA
end
project(A::AbstractArray{T, 2}, ::Val{:CUBIC}, ::OrthonormalBasis{2}) where {T} = project(A, Val(:ISO))

# ── With the frame free ──────────────────────────────────────────────────────

# The harmonics (d₁₂, δ, a, b) of the major-symmetric part, from the
# Kelvin-Mandel matrix in the canonical basis.
function _plane_harmonics(M)
    s2 = sqrt(2 * one(eltype(M)))
    C₁₁₁₂, C₂₂₁₂ = (M[1, 3] + M[3, 1]) / (2 * s2), (M[2, 3] + M[3, 2]) / (2 * s2)
    C₁₁₂₂ = (M[1, 2] + M[2, 1]) / 2
    return C₁₁₁₂ + C₂₂₁₂, (M[1, 1] - M[2, 2]) / 2, C₁₁₁₂ - C₂₂₁₂, (M[1, 1] + M[2, 2] - 2 * C₁₁₂₂ - 2 * M[3, 3]) / 4
end

# The angle of the frame of the closest square tensor: a′ = 0.
function _plane_square_angle(A::AbstractArray{T, 4}) where {T}
    _, _, a, b = _plane_harmonics(_KM_of_array(A))
    return atan(a, b) / 4
end

# The angle of the frame of the closest orthotropic tensor: the least of
# Re(c₁ z + c₂ z²) on the unit circle, z = e^{4iθ}, among the roots of its
# stationarity quartic. A comparison of values and an eigenvalue problem: for
# floating-point numbers.
function _plane_ortho_angle(A::AbstractArray{T, 4}) where {T <: AbstractFloat}
    d₁₂, δ, a, b = _plane_harmonics(_KM_of_array(A))
    c₁, c₂ = complex(d₁₂, δ)^2, complex(a, b)^2
    h(z) = real(c₁ * z + c₂ * z^2)
    # Without the harmonic of order 4 (or of order 2) the least is in closed
    # form: both are kept as candidates, the quartic being singular there.
    best = _unit(-conj(c₁))
    for z in (_unit(sqrt(-conj(c₂))), -_unit(sqrt(-conj(c₂))))
        h(z) < h(best) && (best = z)
    end
    if abs(c₂) > eps(T) * (abs(c₁) + abs(c₂))
        # Companion matrix of z⁴ + (c₁/2c₂) z³ − (c̄₁/2c₂) z − c̄₂/c₂.
        p₃, p₁, p₀ = c₁ / (2 * c₂), -conj(c₁) / (2 * c₂), -conj(c₂) / c₂
        o, e = zero(p₀), one(p₀)
        for r in eigvals([o o o -p₀; e o o -p₁; o e o o; o o e -p₃])
            z = _unit(r)
            h(z) < h(best) && (best = z)
        end
    end
    return angle(best) / 4
end
_unit(z) = iszero(z) ? one(z) : z / abs(z)

# The frame of the plane at the angle θ, a rotated basis whatever the angle:
# the type of a projection does not depend on the frame found.
function _plane_frame(θ)
    c, s = cos(θ), sin(θ)
    return RotatedBasis(SMatrix{2, 2}(c, s, -s, c), Val(2), Val(:typed))
end

# The free projections of the plane, to which the searches of space hand an
# array of dimension 2 (`_proj_tens_opt`). A null tensor is of every class in
# every frame: `atan(0, 0) = 0` gives it the canonical one.
_proj_plane_opt(v::Val{:ORTHO}, A::AbstractArray{T, 4}) where {T <: AbstractFloat} =
    project(A, v, _plane_frame(_plane_ortho_angle(A)))
_proj_plane_opt(v::Val{:CUBIC}, A::AbstractArray{T, 4}) where {T <: AbstractFloat} =
    project(A, v, _plane_frame(_plane_square_angle(A)))
# Order 2: the eigenframe of the symmetric part, tan 2θ = 2A₁₂/(A₁₁ − A₂₂).
_proj_plane_opt(v::Val{:ORTHO}, A::AbstractArray{T, 2}) where {T <: AbstractFloat} =
    project(A, v, _plane_frame(atan(A[1, 2] + A[2, 1], A[1, 1] - A[2, 2]) / 2))

# ── Symmetry search ──────────────────────────────────────────────────────────

_no_plane_ti() = throw(
    ArgumentError(
        "the plane has no transversely isotropic class: a tensor of the plane is isotropic, " *
            "square (`:CUBIC`), orthotropic (`:ORTHO`) or without symmetry"
    )
)

for order in (2, 4)
    @eval begin
        # In the plane the best frame of each class is exact and costs one
        # projection: there is no cheaper candidate to try first, so
        # `optimize_angles` changes nothing.
        function best_sym_tens(
                t::AbstractTens{$order, 2, T}; proj = (:ISO, :CUBIC, :ORTHO), ε = 1.0e-6, optimize_angles::Bool = false
            ) where {T}
            :TI in proj && _no_plane_ti()
            basis = relevant_OrthonormalBasis(get_basis(t))
            newt = express(t, basis)
            A = Array(get_array(newt))
            Q = _rotation_to_canon(basis)
            return _best_sym_loop(newt, proj, ε, sym -> _from_frame(Q, project(A, sym)))
        end
        function best_sym_tens(
                t::AbstractTens{$order, 2, T}, frame::OrthonormalBasis{2}; proj = (:ISO, :CUBIC, :ORTHO), ε = 1.0e-6
            ) where {T}
            :TI in proj && _no_plane_ti()
            basis = relevant_OrthonormalBasis(get_basis(t))
            newt = express(t, basis)
            A = Array(get_array(newt))
            Q = _rotation_to_canon(basis)
            local_frame = _to_frame(Q, frame)
            return _best_sym_loop(
                newt, proj, ε, sym -> sym == :ISO ? project(A, sym) : _from_frame(Q, project(A, sym, local_frame))
            )
        end
    end
end

# ── Average over the rotations of the plane ──────────────────────────────────

"""
    transverse_isotropify(t::AbstractTens{4,2}) → Tens{4,2}
    transverse_isotropify(t::AbstractTens{2,2}) → Tens{2,2}

The exact average of a tensor **of the plane** over the rotations of the plane,
the axis being the normal to the plane. It is not always isotropic: the
rotations alone, without the reflections, leave three tensors of order 4
invariant, and two of order 2,

```math
\\alpha\\,\\mathbb{J} + \\beta\\,\\mathbb{K} + \\omega\\,\\mathbb{W} , \\qquad a\\,\\boldsymbol{1} + c\\,\\boldsymbol{\\epsilon} ,
```

where ``\\mathbb{W}`` turns the deviators by a quarter turn (it has no major
symmetry) and ``\\boldsymbol{\\epsilon}`` is the rotation by a quarter turn of
the plane. The result is a general tensor; [`isotropify`](@ref) gives its
isotropic part ``\\alpha\\,\\mathbb{J} + \\beta\\,\\mathbb{K}`` (or ``a\\,\\boldsymbol{1}``),
the average over the rotations *and* the reflections, and the two coincide for
a tensor with the major symmetry (a symmetric one at order 2).
"""
function transverse_isotropify(t::AbstractTens{4, 2})
    M = _mandel_minor(canonical_components(t))
    r = inv(sqrt(2 * one(eltype(M))))
    # In the basis (spherical, normal deviator, shear) of the Kelvin-Mandel space.
    e₀, e₁, e₂ = SVector(r, r, 0), SVector(r, -r, 0), SVector(0, 0, 1)
    α = transpose(e₀) * M * e₀
    β = (transpose(e₁) * M * e₁ + transpose(e₂) * M * e₂) / 2
    ω = (transpose(e₁) * M * e₂ - transpose(e₂) * M * e₁) / 2
    return from_mandel(
        α * e₀ * transpose(e₀) + β * (e₁ * transpose(e₁) + e₂ * transpose(e₂)) + ω * (e₁ * transpose(e₂) - e₂ * transpose(e₁))
    )
end
function transverse_isotropify(t::AbstractTens{2, 2})
    A = canonical_components(t)
    a, c = (A[1, 1] + A[2, 2]) / 2, (A[1, 2] - A[2, 1]) / 2
    return Tens(Tensor{2, 2}((a, -c, c, a)))
end
