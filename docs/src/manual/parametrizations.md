# [Parametrizations and conversions](@id man-parametrizations)

Which five numbers describe a transversely isotropic material, how to move
between the conventions, and how to go in and out of the Kelvin–Mandel matrix.
The matrices themselves are on
[TI parametrizations](@ref th-ti-parametrizations).

## Isotropic moduli

A 3-D isotropic stiffness is ``\mathbb{C} = 3k\,\mathbb{J} + 2\mu\,\mathbb{K}``:
[`iso_stiffness`](@ref) builds it from the bulk and shear moduli, and
[`iso_stiffness_young_poisson`](@ref) from Young's modulus and Poisson's ratio;
[`bulk_shear`](@ref) and [`young_poisson`](@ref) read them back. A compliance is the inverse,
exact in the class, so its moduli are those of `inv(𝕊)`.

```@example par
using TensND, LinearAlgebra

ℂ = iso_stiffness_young_poisson(70.0, 0.33)   # an aluminum, GPa
bulk_shear(ℂ), young_poisson(ℂ)
```

## The three TI conventions

| Constructor | Extractor | Parameters | Builds |
| :--- | :--- | :--- | :--- |
| [`tens_ti`](@ref) | [`ti_args`](@ref) | ``C_{1111},C_{1122},C_{1133},C_{3333},C_{2323}`` | stiffness *or* compliance |
| [`tens_ti_eng`](@ref) | [`ti_args_eng`](@ref) | ``E_1,E_3,\nu_{12},\nu_{31},G_{31}`` | **compliance** |
| [`tens_ti_hoenig`](@ref) | [`ti_args_hoenig`](@ref) | ``E,\nu_1,\nu_2,H,\Gamma`` | **compliance** |

All three take the symmetry axis as their last argument and return a
`TensTI{4,T,5}`.

```@example par
n = [0.0, 0.0, 1.0]
𝕊 = tens_ti_eng(9.0, 140.0, 0.40, 0.30, 4.6, n)   # a carbon/epoxy ply, GPa
ℂ = inv(𝕊)
ti_args_eng(𝕊)
```

## Converting between them

Conversion is a **round trip through the tensor**: build with one constructor,
read with another extractor. Because all three produce the same object, the
conversion is exact rather than fitted.

```@example par
ti_args(𝕊)          # the same compliance, in component form
```

```@example par
ti_args_hoenig(𝕊)   # and in dimensionless Hoenig ratios
```

``H`` is the axial-to-transverse modulus ratio and ``\Gamma`` the shear
anisotropy; an isotropic material sits at ``H=\Gamma=1`` with ``\nu_1=\nu_2``.

!!! warning "The engineering and Hoenig forms build a compliance"
    `tens_ti_eng` and `tens_ti_hoenig` return ``\mathbb{S}``, not
    ``\mathbb{C}``. To get the stiffness, invert — which is exact and stays in
    the class:

    ```julia
    ℂ = inv(tens_ti_eng(E₁, E₃, ν₁₂, ν₃₁, G₃₁, n))
    ```

    Applying `ti_args_eng` to a *stiffness* returns numbers that are not the
    engineering constants of that material.

## Walpole coefficients

[`get_ℓ`](@ref) reads the six classical coefficients, [`get_ℓ8`](@ref) all
eight. Stiffness and compliance are inverse in the synthetic algebra, so their
``2\times2`` blocks are inverse matrices and their ``\ell_5,\ell_6`` reciprocal:

```@example par
ℓC, ℓS = get_ℓ(ℂ), get_ℓ(𝕊)
LC = [ℓC[1] ℓC[3]; ℓC[4] ℓC[2]]
LS = [ℓS[1] ℓS[3]; ℓS[4] ℓS[2]]
round.(LC * LS, digits = 12)
```

```@example par
round(ℓC[5] * ℓS[5], digits = 12), round(ℓC[6] * ℓS[6], digits = 12)
```

## Kelvin–Mandel round trips

| Function | Direction |
| :--- | :--- |
| [`mandel`](@ref) | tensor → matrix, canonical frame |
| [`mandel_material`](@ref) | `TensOrtho` → matrix, **material** frame |
| [`from_mandel`](@ref) | matrix → tensor |

```@example par
mandel(ℂ)
```

```@example par
norm(get_array(from_mandel(mandel(ℂ))) - get_array(ℂ))
```

## Parameters from a matrix, and back

Four exported aliases give direct access to the projection kernels, for when a
``6\times6`` matrix is the natural input — coming from a file, an experiment or
another library:

| Function | Does |
| :--- | :--- |
| [`ti_params_from_mandel`](@ref)`(C)` | ``6\times6`` → ``(\ell_1,\ell_2,\ell_3,\ell_5,\ell_6)`` |
| [`mandel_from_ti_params`](@ref)`(ℓ…)` | the five coefficients → ``6\times6`` |
| [`ortho_params_from_mandel`](@ref)`(C)` | ``6\times6`` → the nine ``C_{ij}`` |
| [`mandel_from_ortho_params`](@ref)`(C…)` | the nine constants → ``6\times6`` |

They assume the matrix is expressed in the frame where the symmetry axis is
``\underline{e}_3`` (respectively the material frame), and they **project**:
given a matrix that is not exactly of the class, they return the closest one.

```@example par
ti_params_from_mandel(mandel(ℂ))
```

```@example par
norm(mandel_from_ti_params(ti_params_from_mandel(mandel(ℂ))...) - mandel(ℂ))
```

For a matrix in an arbitrary frame, use [`proj_tens`](@ref) instead — it handles
the rotation. See [Projection](@ref man-projection).
