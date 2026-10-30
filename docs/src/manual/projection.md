# [Projection and symmetry detection](@id man-projection)

<!-- crosslinks: written by .github/scripts/crosslinks.py from docs/crosslinks.toml; do not edit -->
!!! tip "See also"
    **Theory**: [Orthotropy](@ref th-orthotropy) · [Cubic symmetry](@ref th-cubic) · [Projection onto a symmetry class](@ref th-projection)

    **Manual**: [Structured tensors](@ref man-structured)

    **Tutorials**: [Projection and symmetry detection](@ref tut-projection)
<!-- /crosslinks -->

Finding the closest tensor of a given material symmetry, and identifying which
symmetry a tensor has. The algorithms are on
[Projection onto a symmetry class](@ref th-projection).

## `project`

```julia
B, d, drel = project(A, sym)              # orientation searched
B, d, drel = project(A, sym, n_or_frame)  # orientation given
```

with `sym ∈ (:ISO, :TI, :ORTHO, :CUBIC)`, `A` an order-2 or order-4 array or
tensor, and

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
C = get_array(ti_from_constants(10.0, 3.0, 2.5, 12.0, 2.0, n))

B, d, drel = project(C, :TI, n)
(typeof(B), drel)
```

Projecting onto a class the tensor already belongs to is exact. Projecting about
the wrong axis is not:

```@example proj
[round(project(C, :TI, [sind(α), 0.0, cosd(α)])[3], digits = 5) for α in (0, 15, 30, 45, 90)]
```

That function of the orientation is what the search minimizes.

## Orientation given, or found

The three-argument form takes the axis or frame; the two-argument form searches
it:

```@example proj
θ, ϕ = 0.6, 1.1
n_tilt = [sin(θ)cos(ϕ), sin(θ)sin(ϕ), cos(θ)]
C_tilt = get_array(ti_from_constants(10.0, 3.0, 2.5, 12.0, 2.0, n_tilt))

B_opt, _, drel_opt = project(C_tilt, :TI)
(round.(collect(axis(B_opt)), digits = 8), drel_opt)
```

The axis is recovered up to a sign, which is immaterial: ``\underline{n}`` and
``-\underline{n}`` define the same symmetry. Similarly an orthotropic frame is
determined only up to the reorderings and sign flips that map the three
orthogonal planes onto themselves.

The search is **deterministic** — a screening of fixed orientations, then
Newton's method on the angles from the best ones, with no random step — so
repeated calls return bit-identical results:

```@example proj
length(unique([project(C_tilt, :TI)[3] for _ in 1:5]))
```

## `best_sym_tens`

Tries the classes from the most restrictive to the least and returns the first
whose relative error falls below `ε`:

```julia
B, d, drel, sym = best_sym_tens(t; proj = (:ISO, :CUBIC, :TI, :ORTHO), ε = 1e-6,
                                optimize_angles = false)
```

`sym` is one of `:ISO`, `:CUBIC`, `:TI`, `:ORTHO`, `:ANISO`. The default order
is by number of constants — 2, 3, 5, 9 — so the tightest class that fits within
`ε` is the one reported. `:CUBIC` and `:TI` are **incomparable**, neither
containing the other, so a tensor satisfying both is reported cubic, having the
fewer constants. The argument must be an
`AbstractTens`, not a bare array — wrap a raw array with `Tens` first.

```@example proj
best_sym_tens(Tens(C_tilt))[4]
```

| `optimize_angles` | axis / frame |
| :--- | :--- |
| `false` (default) | from the tensor if it is a structured container, otherwise from the Kelvin–Mandel eigenstructure |
| `true` | searched, as by `project(A, sym)` |

!!! note "The cheap path is usually enough"
    The eigenstructure candidate is **exact** whenever the tensor genuinely has
    the symmetry sought, so the default path identifies tilted TI and rotated
    orthotropic tensors correctly and roughly an order of magnitude faster.
    Reach for `optimize_angles = true` when the tensor is only *approximately*
    of the class and the best orientation is itself the question.

## Cubic symmetry

`:CUBIC` takes a cube frame, like `:ORTHO`, and returns a
[`TensCubic`](@ref) — or a [`TensISO`](@ref) at order 2, the two classes
coinciding there.

```@example proj
ℬ = CanonicalBasis{3, Float64}()
cub = cubic_from_constants(10.0, 4.0, 2.0, ℬ)
B, d, drel = project(get_array(cub), :CUBIC, ℬ)
drel
```

Omitting the frame searches the cube orientation, exactly as for a TI axis or
an orthotropic frame. The octahedral group
is discrete but the *orientation* of the cube is not — it is an ordinary
rotation — so the objective is smooth in the Euler angles; the group shows up
only as a 24-fold degeneracy of the minimum, and all 24 frames describe the same
tensor with the same coefficients.

The residual is the number worth reporting. When the morphology and the medium
leave the octahedral group invariant, the answer belongs to the class **by group
theory**, so `drel` is bounded by the discretization error and by nothing else —
an error estimate with no reference solution in it. What lives *inside* the
class, and is therefore not measured by `drel`, is
[`cubic_anisotropy`](@ref): the two together separate a real material anisotropy
from a numerical artifact, since an artifact breaks the symmetry while a real
anisotropy does not.

## Predicates

[`isisotropic`](@ref), [`istransverselyisotropic`](@ref), [`isorthotropic`](@ref) and [`iscubic`](@ref) are the same
computation with a boolean answer, and they respect the hierarchy
ISO ⊂ TI ⊂ ORTHO:

```@example proj
𝕀, 𝕁, 𝕂 = iso_projectors(Val(3), Val(Float64))
(isisotropic(𝕀), istransverselyisotropic(𝕀), isorthotropic(𝕀))
```

```@example proj
(isisotropic(C), istransverselyisotropic(C), isorthotropic(C))
```

Each accepts an optional axis or frame, and the `optimize_angles` keyword.

## Two things worth stating in a report

**Which tensor was projected.** The Euclidean distance is not invariant under
inversion, so projecting a stiffness and projecting the corresponding compliance
give different materials:

```@example proj
ort = TensOrtho(10.0, 8.0, 9.0, 3.0, 2.0, 4.0, 2.5, 3.0, 1.5, CanonicalBasis{3, Float64}())
Biso, _, _ = project(get_array(ort), :ISO)
Siso, _, _ = project(get_array(inv(ort)), :ISO)
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
