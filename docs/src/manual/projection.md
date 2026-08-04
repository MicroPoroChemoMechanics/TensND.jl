# [Projection and symmetry detection](@id man-projection)

Finding the closest tensor of a given material symmetry, and identifying which
symmetry a tensor has. The algorithms are on
[Projection onto a symmetry class](@ref th-projection).

## `proj_tens`

```julia
B, d, drel = proj_tens(sym, A)              # orientation optimized — needs NLopt
B, d, drel = proj_tens(sym, A, n_or_frame)  # orientation given
```

with `sym ∈ (:ISO, :TI, :ORTHO)`, `A` an order-2 or order-4 array or tensor, and

| Returned | Meaning |
| :--- | :--- |
| `B` | the projection, in the storage type of the class |
| `d` | absolute Frobenius distance ``\|B-A\|`` |
| `drel` | relative distance ``d/\|A\|`` |

**Always look at `drel`, not `d`.** The relative distance is dimensionless, so a
tolerance on it is independent of the units of the moduli.

```@example proj
using TensND, LinearAlgebra

n = [0.0, 0.0, 1.0]
C = get_array(tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, n))

B, d, drel = proj_tens(:TI, C, n)
(typeof(B), drel)
```

Projecting onto a class the tensor already belongs to is exact. Projecting about
the wrong axis is not:

```@example proj
[round(proj_tens(:TI, C, [sind(α), 0.0, cosd(α)])[3], digits = 5) for α in (0, 15, 30, 45, 90)]
```

That function of the orientation is what the optimizer minimizes.

## Orientation given, or found

Without `NLopt`, the axis or frame must be supplied. With it, the two-argument
form searches:

```@example proj
using NLopt

θ, ϕ = 0.6, 1.1
n_tilt = [sin(θ)cos(ϕ), sin(θ)sin(ϕ), cos(θ)]
C_tilt = get_array(tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, n_tilt))

B_opt, _, drel_opt = proj_tens(:TI, C_tilt)
(round.(collect(axis(B_opt)), digits = 8), drel_opt)
```

The axis is recovered up to a sign, which is immaterial: ``\underline{n}`` and
``-\underline{n}`` define the same symmetry. Similarly an orthotropic frame is
determined only up to the reorderings and sign flips that map the three
orthogonal planes onto themselves.

The search is a **deterministic** multi-start, so repeated calls return
bit-identical results:

```@example proj
length(unique([proj_tens(:TI, C_tilt)[3] for _ in 1:5]))
```

## `best_sym_tens`

Tries the classes from the most restrictive to the least and returns the first
whose relative error falls below `ε`:

```julia
B, d, drel, sym = best_sym_tens(t; proj = (:ISO, :TI, :ORTHO), ε = 1e-6,
                                optimize_angles = false)
```

`sym` is one of `:ISO`, `:TI`, `:ORTHO`, `:ANISO`. The argument must be an
`AbstractTens`, not a bare array — wrap a raw array with `Tens` first.

```@example proj
best_sym_tens(Tens(C_tilt))[4]
```

| `optimize_angles` | axis / frame | needs NLopt |
| :--- | :--- | :--- |
| `false` (default) | from the tensor if it is a structured container, otherwise from the Kelvin–Mandel eigenstructure | no |
| `true` | found by multi-start optimization | yes |

!!! note "The cheap path is usually enough"
    The eigenstructure candidate is **exact** whenever the tensor genuinely has
    the symmetry sought, so the default path identifies tilted TI and rotated
    orthotropic tensors correctly and roughly an order of magnitude faster.
    Reach for `optimize_angles = true` when the tensor is only *approximately*
    of the class and the best orientation is itself the question.

## Predicates

[`is_ISO`](@ref), [`is_TI`](@ref) and [`is_ORTHO`](@ref) are the same
computation with a boolean answer, and they respect the hierarchy
ISO ⊂ TI ⊂ ORTHO:

```@example proj
𝕀, 𝕁, 𝕂 = ISO(Val(3), Val(Float64))
(is_ISO(𝕀), is_TI(𝕀), is_ORTHO(𝕀))
```

```@example proj
(is_ISO(C), is_TI(C), is_ORTHO(C))
```

Each accepts an optional axis or frame, and the `optimize_angles` keyword.

## Two things worth stating in a report

**Which tensor was projected.** The Euclidean distance is not invariant under
inversion, so projecting a stiffness and projecting the corresponding compliance
give different materials:

```@example proj
ort = TensOrtho(10.0, 8.0, 9.0, 3.0, 2.0, 4.0, 2.5, 3.0, 1.5, CanonicalBasis{3, Float64}())
Biso, _, _ = proj_tens(:ISO, get_array(ort))
Siso, _, _ = proj_tens(:ISO, get_array(inv(ort)))
round(norm(get_array(inv(Biso)) - get_array(Siso)) / norm(get_array(Siso)), digits = 4)
```

Alternatives invariant under inversion are discussed in
[Isotropic tensors](@ref th-isotropic).

**Whether you wanted a projection at all.** Projecting onto the five-parameter
major-symmetric TI subspace is a *best fit*: it forces ``\ell_3=\ell_4`` and
discards ``\ell_7,\ell_8``. If the object is not major-symmetric — a
strain-concentration tensor, typically — an exact **average** over the rotation
group preserves that content where a best fit destroys it. The distinction is
developed on [The extended Walpole algebra](@ref th-walpole-extended) and
[Projection onto a symmetry class](@ref th-projection).
