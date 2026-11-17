# ── Solid harmonics, vector harmonics and the solutions of the Lamé equation ──
#
# With (r, θ, ϕ) the spherical coordinates of x about 𝐞₃ and
# χ_n^m = (-1)^m P_n^m(cos θ) e^{imϕ} the surface harmonic of degree n and order
# 0 ≤ m ≤ n (P_n^m of `legendre_P`, P_1^1(p) = -√(1-p²), so χ_1^1 = sin θ e^{iϕ}),
# the complex solid harmonics are normalized as in kushch2013, Eq. (2.9):
#
#     y_n^m(x) = rⁿ χ_n^m / (n+m)! ,     Y_n^m(x) = (n-m)! χ_n^m / r^{n+1} ,
#
# the normalization in which a decaying harmonic translates to another center
# with decaying harmonics as coefficients and no other factor,
#
#     Y_{n'}^{m'}(r + R) = Σ_{n,m} (-1)^{n+m} Y_{n'+n}^{m'-m}(R) y_n^m(r) ,  ‖r‖ < ‖R‖ .
#
# They are the meridian harmonics of `interior_harmonics` and
# `exterior_harmonics` in the spherical chart, times e^{imϕ}, (-1)^m and
# 1/(n+m)! or (n-m)!. Both families are evaluated by recurrences on the
# Cartesian components, with w = x₁ + i x₂,
#
#     y_m^m     = y_{m-1}^{m-1} w / (2m) ,
#     y_{n+1}^m = [(2n+1) x₃ y_n^m - r² y_{n-1}^m] / [(n-m+1)(n+m+1)] ,
#     Y_m^m     = Y_{m-1}^{m-1} (2m-1) w / r² ,
#     Y_{n+1}^m = [(2n+1) x₃ Y_n^m - (n-m)(n+m) Y_{n-1}^m] / r² ,
#
# with neither an angle nor a factorial, valid at the origin for y. Negative
# orders follow from h_n^{-m} = (-1)^m conj(h_n^m).

"""
    SolidHarmonicTable{C}

Values of a family of complex solid harmonics at one point, for the degrees
``0 \\le n \\le N``: `h[n, m]` for ``|m| \\le n``. Only the orders ``m \\ge 0`` are
stored, the others following from ``h_n^{-m} = (-1)^m\\,\\overline{h_n^m}``; an
order beyond the degree reads as zero, as the translation formulas need. The
fields are `nmax` (``N``) and `values`, which holds ``h_n^m``, ``m \\ge 0``, at
the position [`solid_harmonic_index`](@ref)`(n, m)`: tables of one degree add,
scale and accumulate on their `values`.

See also [`regular_solid_harmonics`](@ref), [`irregular_solid_harmonics`](@ref).
"""
struct SolidHarmonicTable{C <: Complex}
    nmax::Int
    values::Vector{C}
end

"""
    SolidHarmonicTable(f, N)

The table of the values `f(n, m)` for ``0 \\le m \\le n \\le N``, the negative
orders following from the conjugates: a family of harmonics that shares the
symmetry ``h_n^{-m} = (-1)^m\\,\\overline{h_n^m}`` of the solid harmonics, read
with the same indexing and by the same vector harmonics.

# Examples
```jldoctest
julia> h = SolidHarmonicTable((n, m) -> complex(n, m), 2);

julia> h[2, -1]
-2 + 1im
```
"""
SolidHarmonicTable(f::F, N::Integer) where {F} = SolidHarmonicTable(Int(N), [complex(f(n, m)) for n in 0:N for m in 0:n])

# The position of (n, m ≥ 0) in the packed storage.
@inline _packed(n::Integer, m::Integer) = (n * (n + 1)) >> 1 + m + 1

"""
    solid_harmonic_index(n, m) -> Int

The position of ``h_n^m`` (``0 \\le m \\le n``) in the `values` of a
[`SolidHarmonicTable`](@ref): the degrees in increasing order and, within a
degree, the orders ``0, \\dots, n``; a table of degree ``N`` holds
`solid_harmonic_index(N, N)` values.
"""
solid_harmonic_index(n::Integer, m::Integer) =
    0 ≤ m ≤ n ? _packed(n, m) : throw(DomainError((n, m), "a stored harmonic has 0 ≤ m ≤ n"))

function Base.getindex(h::SolidHarmonicTable{C}, n::Integer, m::Integer) where {C}
    0 ≤ n ≤ h.nmax || throw(BoundsError(h, (n, m)))
    abs(m) > n && return zero(C)
    v = h.values[_packed(n, abs(m))]
    return m ≥ 0 ? v : (isodd(m) ? -conj(v) : conj(v))
end

# A harmonic of a table, zero beyond the degrees it holds.
@inline _entry(h::SolidHarmonicTable{C}, n::Integer, m::Integer) where {C} = 0 ≤ n ≤ h.nmax ? h[n, m] : zero(C)

# The sectoral terms (n = m), then each order upward in the degree: recurrences.
function _solid_table(x::AbstractVector, nmax::Integer, sectoral::F, step::G) where {F, G}
    length(x) == 3 || throw(DimensionMismatch("solid harmonics are defined in 3D; got a $(length(x))-component point"))
    nmax ≥ 0 || throw(DomainError(nmax, "the highest degree of the solid harmonics must be non-negative"))
    T = float(eltype(x))
    w = complex(T(x[1]), T(x[2]))
    z = T(x[3])
    r² = abs2(w) + z * z
    v = Vector{Complex{T}}(undef, _packed(nmax, nmax))
    v[1] = sectoral(0, zero(Complex{T}), w, r²)
    for m in 0:nmax
        m > 0 && (v[_packed(m, m)] = sectoral(m, v[_packed(m - 1, m - 1)], w, r²))
        prev = zero(Complex{T})
        for n in m:(nmax - 1)
            cur = v[_packed(n, m)]
            v[_packed(n + 1, m)] = step(n, m, cur, prev, z, r²)
            prev = cur
        end
    end
    return SolidHarmonicTable(Int(nmax), v)
end

"""
    regular_solid_harmonics(x, N) -> SolidHarmonicTable

Regular solid harmonics
``y_n^m(\\underline{x}) = (-1)^m\\,r^n P_n^m(\\cos\\theta)\\,e^{\\mathrm{i}m\\varphi}/(n+m)!``
at the point `x`, up to the degree `N` [kushch2013](@cite). They are
homogeneous polynomials of the coordinates, defined at the origin too, and the
coordinates may be numbers, duals or symbols.

# Examples
```jldoctest
julia> y = regular_solid_harmonics([0.0, 0.0, 2.0], 2);

julia> y[2, 0] ≈ 2.0^2 / 2
true
```

See also [`irregular_solid_harmonics`](@ref), [`interior_harmonics`](@ref).
"""
regular_solid_harmonics(x::AbstractVector, N::Integer) = _solid_table(
    x, N,
    (m, prev, w, r²) -> m == 0 ? one(w) : prev * w / (2m),
    (n, m, cur, prev, z, r²) -> ((2n + 1) * z * cur - r² * prev) / ((n - m + 1) * (n + m + 1))
)

"""
    irregular_solid_harmonics(x, N) -> SolidHarmonicTable

Decaying solid harmonics
``Y_n^m(\\underline{x}) = (-1)^m\\,(n-m)!\\,P_n^m(\\cos\\theta)\\,e^{\\mathrm{i}m\\varphi}/r^{n+1}``
at the point `x ≠ 0`, up to the degree `N` [kushch2013](@cite).

# Examples
```jldoctest
julia> Y = irregular_solid_harmonics([0.0, 0.0, 2.0], 1);

julia> Y[1, 0] ≈ 1 / 2.0^2
true
```

See also [`regular_solid_harmonics`](@ref), [`exterior_harmonics`](@ref).
"""
function irregular_solid_harmonics(x::AbstractVector, N::Integer)
    all(iszero, x) && throw(DomainError(x, "the decaying solid harmonics are singular at the origin"))
    return _solid_table(
        x, N,
        (m, prev, w, r²) -> m == 0 ? complex(inv(sqrt(r²))) : prev * (2m - 1) * w / r²,
        (n, m, cur, prev, z, r²) -> ((2n + 1) * z * cur - (n - m) * (n + m) * prev) / r²
    )
end

"""
    solid_harmonic_norm(n, m) -> Float64

``\\sqrt{(n-m)!\\,(n+m)!}``, the factor that brings the coefficients of an
expansion in the solid harmonics of [`regular_solid_harmonics`](@ref) to the
same order at every degree [kushch2013; Eq. (2.52)](@cite); in floating point,
from the logarithms of the factorials, beyond the range of the factorials.
"""
solid_harmonic_norm(n::Integer, m::Integer) = exp((logfactorial(n - m) + logfactorial(n + m)) / 2)

# ── Vector harmonics ─────────────────────────────────────────────────────────
#
# In the complex Cartesian basis 𝐞₊ = (𝐢₁ + i 𝐢₂)/2, 𝐞₋ = conj(𝐞₊), 𝐞₃ = 𝐢₃, the
# three families of harmonic vectors (kushch2013, Eqs. (4.14), (4.15)) are
#
#   Y⁽¹⁾_nm = 𝐞₊ Y_{n+1}^{m-1} - 𝐞₋ Y_{n+1}^{m+1} - 𝐞₃ Y_{n+1}^m              (= grad Y_n^m)
#   Y⁽²⁾_nm = (i/n) [𝐞₊ (n+m) Y_n^{m-1} + 𝐞₋ (n-m) Y_n^{m+1} - 𝐞₃ m Y_n^m]
#   Y⁽³⁾_nm = [𝐞₊ (n+m-1)(n+m) Y_{n-1}^{m-1} - 𝐞₋ (n-m-1)(n-m) Y_{n-1}^{m+1}
#             + 𝐞₃ (n-m)(n+m) Y_{n-1}^m] / [n(2n-1)]
#
#   y⁽¹⁾_nm = 𝐞₊ y_{n-1}^{m-1} - 𝐞₋ y_{n-1}^{m+1} + 𝐞₃ y_{n-1}^m              (= grad y_n^m)
#   y⁽²⁾_nm = (i/(n+1)) [-𝐞₊ (n-m+1) y_n^{m-1} - 𝐞₋ (n+m+1) y_n^{m+1} + 𝐞₃ m y_n^m]
#   y⁽³⁾_nm = [𝐞₊ (n-m+1)(n-m+2) y_{n+1}^{m-1} - 𝐞₋ (n+m+1)(n+m+2) y_{n+1}^{m+1}
#             - 𝐞₃ (n-m+1)(n+m+1) y_{n+1}^m] / [(n+1)(2n+3)]
#
# The first is a gradient, the second a rotation field (Y⁽²⁾_nm = grad Y_n^m × x/n),
# both of zero divergence; the third has the divergence -Y_n^m (-y_n^m).

# 𝐞₊ a + 𝐞₋ b + 𝐞₃ c in Cartesian components.
@inline _complex_cartesian(a, b, c) = SVector((a + b) / 2, im * (a - b) / 2, c)

"""
    irregular_vector_harmonic(i, n, m, Y::SolidHarmonicTable) -> SVector{3}

Cartesian components of the decaying vector solid harmonic
``\\underline{Y}^{(i)}_{nm}`` (``i = 1, 2, 3``; ``n \\ge 1`` for ``i > 1``) of
[kushch2013](@citet), Eqs. (4.14)–(4.15), from the table `Y` of
[`irregular_solid_harmonics`](@ref) at the point, which reaches the degree
``n + 1``: ``\\underline{Y}^{(1)}_{nm} = \\underline{\\operatorname{grad}}\\,Y_n^m``,
``\\underline{Y}^{(2)}_{nm}`` tangent to the spheres, and
``\\underline{Y}^{(3)}_{nm}`` of divergence ``-Y_n^m``.

See also [`regular_vector_harmonic`](@ref), [`irregular_lame_solution`](@ref).
"""
function irregular_vector_harmonic(i::Integer, n::Integer, m::Integer, Y::SolidHarmonicTable)
    i == 1 && return _complex_cartesian(_entry(Y, n + 1, m - 1), -_entry(Y, n + 1, m + 1), -_entry(Y, n + 1, m))
    i == 2 && return (im / n) * _complex_cartesian(
        (n + m) * _entry(Y, n, m - 1), (n - m) * _entry(Y, n, m + 1), -m * _entry(Y, n, m)
    )
    i == 3 || throw(DomainError(i, "the family of a vector solid harmonic is 1, 2 or 3"))
    return _complex_cartesian(
        (n + m - 1) * (n + m) * _entry(Y, n - 1, m - 1), -(n - m - 1) * (n - m) * _entry(Y, n - 1, m + 1),
        (n - m) * (n + m) * _entry(Y, n - 1, m)
    ) / (n * (2n - 1))
end

"""
    regular_vector_harmonic(i, n, m, y::SolidHarmonicTable) -> SVector{3}

Cartesian components of the regular vector solid harmonic
``\\underline{y}^{(i)}_{nm}`` (``i = 1, 2, 3``) of [kushch2013](@citet), from the
table `y` of [`regular_solid_harmonics`](@ref) at the point, which reaches the
degree ``n + 1``.

See also [`irregular_vector_harmonic`](@ref), [`regular_lame_solution`](@ref).
"""
function regular_vector_harmonic(i::Integer, n::Integer, m::Integer, y::SolidHarmonicTable)
    i == 1 && return _complex_cartesian(_entry(y, n - 1, m - 1), -_entry(y, n - 1, m + 1), _entry(y, n - 1, m))
    i == 2 && return (im / (n + 1)) * _complex_cartesian(
        -(n - m + 1) * _entry(y, n, m - 1), -(n + m + 1) * _entry(y, n, m + 1), m * _entry(y, n, m)
    )
    i == 3 || throw(DomainError(i, "the family of a vector solid harmonic is 1, 2 or 3"))
    return _complex_cartesian(
        (n - m + 1) * (n - m + 2) * _entry(y, n + 1, m - 1), -(n + m + 1) * (n + m + 2) * _entry(y, n + 1, m + 1),
        -(n - m + 1) * (n + m + 1) * _entry(y, n + 1, m)
    ) / ((n + 1) * (2n + 3))
end

# ── Solutions of the Lamé equation ───────────────────────────────────────────
#
# The displacements of an isotropic solid of Poisson ratio ν with no body force
# are the first two families as they are and the Papkovich–Neuber field of the
# third harmonic vector,
#
#   U⁽¹⁾ = Y⁽¹⁾ ,  U⁽²⁾ = Y⁽²⁾ ,  U⁽³⁾ = 4(1-ν) Y⁽³⁾ - grad(x·Y⁽³⁾) ,
#
# with x·Y⁽³⁾_nm = r² Y_n^m/(2n-1) and x·y⁽³⁾_nm = -r² y_n^m/(2n+3) (kushch2013,
# Eq. (4.27), written without a biharmonic function).

"""
    irregular_lame_solution(i, n, m, x, ν[, Y]) -> SVector{3}

Decaying solution ``\\underline{U}^{(i)}_{nm}`` of the Lamé equation of an
isotropic solid of Poisson ratio `ν` at the point `x`: the vector harmonics
``\\underline{Y}^{(1)}_{nm}``, ``\\underline{Y}^{(2)}_{nm}`` and
``\\underline{U}^{(3)}_{nm} = 4(1-\\nu)\\,\\underline{Y}^{(3)}_{nm}
- \\underline{\\operatorname{grad}}\\,(\\underline{x}\\cdot\\underline{Y}^{(3)}_{nm})``
[kushch2013](@cite). `Y`, the table of [`irregular_solid_harmonics`](@ref) at
`x` up to a degree ``n + 1`` or more, is computed when not given.

See also [`regular_lame_solution`](@ref), [`lame_solution_translation`](@ref).
"""
irregular_lame_solution(i::Integer, n::Integer, m::Integer, x::AbstractVector, ν) =
    irregular_lame_solution(i, n, m, x, ν, irregular_solid_harmonics(x, n + 1))
function irregular_lame_solution(i::Integer, n::Integer, m::Integer, x::AbstractVector, ν, Y::SolidHarmonicTable)
    i < 3 && return irregular_vector_harmonic(i, n, m, Y)
    r² = sum(abs2, x)
    # grad(x·Y⁽³⁾) = grad(r² Y_n^m)/(2n-1) = (2 x Y_n^m + r² Y⁽¹⁾_nm)/(2n-1)
    grad = (2 * Y[n, m] * SVector{3}(x) + r² * irregular_vector_harmonic(1, n, m, Y)) / (2n - 1)
    return 4 * (1 - ν) * irregular_vector_harmonic(3, n, m, Y) - grad
end

"""
    regular_lame_solution(i, n, m, x, ν[, y]) -> SVector{3}

Regular solution ``\\underline{u}^{(i)}_{nm}`` of the Lamé equation at the point
`x`, for the Poisson ratio `ν`: the counterpart of
[`irregular_lame_solution`](@ref) on the regular vector harmonics. The
solutions of the lowest degrees are the dilatation ``\\underline{u}^{(3)}_{00}``,
the translations ``\\underline{u}^{(1)}_{1m}``, the rotations
``\\underline{u}^{(2)}_{1m}`` and the uniform deviatoric strains
``\\underline{u}^{(1)}_{2m}``. `y`, the table of
[`regular_solid_harmonics`](@ref) at `x` up to a degree ``n + 1`` or more, is
computed when not given.
"""
regular_lame_solution(i::Integer, n::Integer, m::Integer, x::AbstractVector, ν) =
    regular_lame_solution(i, n, m, x, ν, regular_solid_harmonics(x, n + 1))
function regular_lame_solution(i::Integer, n::Integer, m::Integer, x::AbstractVector, ν, y::SolidHarmonicTable)
    i < 3 && return regular_vector_harmonic(i, n, m, y)
    r² = sum(abs2, x)
    # grad(x·y⁽³⁾) = -grad(r² y_n^m)/(2n+3) = -(2 x y_n^m + r² y⁽¹⁾_nm)/(2n+3)
    grad = -(2 * y[n, m] * SVector{3}(x) + r² * regular_vector_harmonic(1, n, m, y)) / (2n + 3)
    return 4 * (1 - ν) * regular_vector_harmonic(3, n, m, y) - grad
end

# ── Translation of the Lamé solutions ────────────────────────────────────────
#
# About a center at R from the first, a decaying vector harmonic expands on the
# regular ones (kushch2013, Eq. (4.85)):
#
#   Y⁽ⁱ⁾_n'm'(r + R) = Σ_{j ≤ i} Σ_{n,m} (-1)^{n+m} α⁽ⁱ⁾⁽ʲ⁾_n'nm'm(R) y⁽ʲ⁾_nm(r) ,
#
#   α⁽¹⁾⁽¹⁾ = α⁽²⁾⁽²⁾ = α⁽³⁾⁽³⁾ = Y_{n'+n}^{m'-m} ,
#   α⁽²⁾⁽¹⁾ = i (m/n + m'/n') Y_{n'+n-1}^{m'-m} ,           α⁽³⁾⁽²⁾ = -α⁽²⁾⁽¹⁾ ,
#   α⁽³⁾⁽¹⁾ = [1 + m m'/(n n') - (n'²-m'²)/(n'(2n'-1)) - (n²-m²)/(n(2n-1))] Y_{n'+n-2}^{m'-m} .
#
# With x = r + R, U⁽³⁾ = [4(1-ν) B - grad(r·B)] - grad(R·B), B = Y⁽³⁾_n'm'(r + R);
# B expanded as above, the bracket is term by term a regular Lamé solution:
# u⁽³⁾_nm for y⁽³⁾_nm, 4(1-ν) u⁽²⁾_nm for y⁽²⁾_nm (r·y⁽²⁾ = 0) and
# (4(1-ν) - n) u⁽¹⁾_nm for y⁽¹⁾_nm = grad y_n^m (r·grad y_n^m = n y_n^m). The
# last term is the gradient of the harmonic R·B, which translates as
#
#   R·Y⁽³⁾_n'm'(r + R) = Σ_{n,m} (-1)^{n+m} H_n'nm'm(R) y_n^m(r) ,
#   H_n'nm'm = [(R·𝐞₊)(n'+m'-1)(n'+m') Y_{n'-1+n}^{m'-1-m} - (R·𝐞₋)(n'-m'-1)(n'-m') Y_{n'-1+n}^{m'+1-m}
#              + R₃ (n'-m')(n'+m') Y_{n'-1+n}^{m'-m}] / [n'(2n'-1)] .
#
# Hence η = α for i < 3, and η⁽³⁾⁽³⁾ = α⁽³⁾⁽³⁾, η⁽³⁾⁽²⁾ = 4(1-ν) α⁽³⁾⁽²⁾,
# η⁽³⁾⁽¹⁾ = (4(1-ν) - n) α⁽³⁾⁽¹⁾ - H: each coefficient is affine in 4(1-ν).
#
# Over the images of a periodic array, Y_N^M(R) is replaced by the lattice sums
# of `lattice_solid_harmonics` and the products R_k Y_N^M(R) of H by their
# periodic counterparts, from
#
#   (x₁ + i x₂)/2 · Y_N^M = [ r² Y_{N+1}^{M+1} - (N-M)(N-M-1) Y_{N-1}^{M+1}] / [2(2N+1)] ,
#   (x₁ - i x₂)/2 · Y_N^M = [-r² Y_{N+1}^{M-1} + (N+M)(N+M-1) Y_{N-1}^{M-1}] / [2(2N+1)] ,
#   x₃ · Y_N^M           = [ r² Y_{N+1}^M + (N-M)(N+M) Y_{N-1}^M] / (2N+1) ,
#
# in which r² Y_K^μ becomes -(2K-1) B_K^μ, B the table of `lattice_solid_biharmonics`.

"""
    lame_solution_translation(i, j, n′, n, m′, m, Y::SolidHarmonicTable, R) -> (a, b)
    lame_solution_translation(i, j, n′, n, m′, m, Y::SolidHarmonicTable, B::SolidHarmonicTable) -> (a, b)

Coefficient ``\\eta^{(i)(j)}_{n'n\\,m'm} = a + 4(1-\\nu)\\,b`` of the regular Lamé
solution ``\\underline{u}^{(j)}_{nm}`` in the expansion of the decaying one
``\\underline{U}^{(i)}_{n'm'}`` about a center at `R` from its own,

```math
\\underline{U}^{(i)}_{n'm'}(\\underline{r} + \\underline{R}) = \\sum_{j \\le i}\\sum_{n,m}
  (-1)^{n+m}\\,\\eta^{(i)(j)}_{n'n\\,m'm}(\\underline{R})\\;\\underline{u}^{(j)}_{nm}(\\underline{r}),
\\qquad \\|\\underline{r}\\| < \\|\\underline{R}\\| ,
```

from the table `Y` of [`irregular_solid_harmonics`](@ref) at `R`, which reaches
the degree ``n' + n`` [kushch2013; Eq. (4.85)](@cite). The two parts are
returned apart, the Poisson ratio of the solid entering through ``4(1-\\nu)``
alone.

With the tables `Y` and `B` of [`lattice_solid_harmonics`](@ref) and
[`lattice_solid_biharmonics`](@ref) in place of `Y` and `R`, the coefficient is
that of the sum of ``\\underline{U}^{(i)}_{n'm'}`` over the images of a periodic
array, apart from the uniform terms of its background at the lowest degrees.
"""
function lame_solution_translation(
        i::Integer, j::Integer, n′::Integer, n::Integer, m′::Integer, m::Integer,
        Y::SolidHarmonicTable{C}, W::Union{AbstractVector, SolidHarmonicTable}
    ) where {C}
    z = zero(C)
    (j > i || (n == 0 && j < 3)) && return (z, z)
    i == j && return (_entry(Y, n′ + n, m′ - m), z)
    if i - j == 1
        α = im * (m / n + m′ / n′) * _entry(Y, n′ + n - 1, m′ - m)      # α⁽²⁾⁽¹⁾ = -α⁽³⁾⁽²⁾
        return i == 2 ? (α, z) : (z, -α)
    end
    α = (1 + m * m′ / (n * n′) - (n′^2 - m′^2) / (n′ * (2n′ - 1)) - (n^2 - m^2) / (n * (2n - 1))) *
        _entry(Y, n′ + n - 2, m′ - m)
    N = n′ - 1 + n
    H = (
        (n′ + m′ - 1) * (n′ + m′) * _weighted_plus(W, Y, N, m′ - 1 - m) -
            (n′ - m′ - 1) * (n′ - m′) * _weighted_minus(W, Y, N, m′ + 1 - m) +
            (n′ - m′) * (n′ + m′) * _weighted_3(W, Y, N, m′ - m)
    ) / (n′ * (2n′ - 1))
    return (-n * α - H, α)
end

# (R₁ + i R₂)/2 · Y_N^M, (R₁ - i R₂)/2 · Y_N^M and R₃ Y_N^M at the point R of a
# single source, and their sums over the images of a periodic array.
_weighted_plus(R::AbstractVector, Y, N, M) = (R[1] + im * R[2]) / 2 * _entry(Y, N, M)
_weighted_minus(R::AbstractVector, Y, N, M) = (R[1] - im * R[2]) / 2 * _entry(Y, N, M)
_weighted_3(R::AbstractVector, Y, N, M) = R[3] * _entry(Y, N, M)
_weighted_plus(B::SolidHarmonicTable, Y, N, M) =
    -(_entry(B, N + 1, M + 1) + (N - M) * (N - M - 1) * _entry(Y, N - 1, M + 1) / (2N + 1)) / 2
_weighted_minus(B::SolidHarmonicTable, Y, N, M) =
    (_entry(B, N + 1, M - 1) + (N + M) * (N + M - 1) * _entry(Y, N - 1, M - 1) / (2N + 1)) / 2
_weighted_3(B::SolidHarmonicTable, Y, N, M) =
    -_entry(B, N + 1, M) + (N - M) * (N + M) * _entry(Y, N - 1, M) / (2N + 1)

export SolidHarmonicTable, solid_harmonic_index, regular_solid_harmonics, irregular_solid_harmonics, solid_harmonic_norm
export regular_vector_harmonic, irregular_vector_harmonic
export regular_lame_solution, irregular_lame_solution, lame_solution_translation
