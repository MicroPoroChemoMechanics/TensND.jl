# [Parametrizations and conversions](@id man-parametrizations)

!!! tip "See also"
    **Theory**: [The Walpole basis](@ref th-walpole) · [Transversely isotropic parametrizations](@ref th-ti-parametrizations) · [The extended Walpole algebra](@ref th-walpole-extended) · [Isotropic tensors](@ref th-isotropic)

    **Manual**: [Structured tensors](@ref man-structured) · [Symbolic tensors of a material symmetry](@ref man-symbolic-materials) · [Elastic moduli](@ref man-elastic-moduli)

    **Tutorials**: [The Walpole basis end to end](@ref tut-walpole) · [Acoustic and Hill tensors](@ref tut-hill-tensors) · [The eight-dimensional axially invariant space](@ref tut-walpole-extended)

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
Every other pair of the six isotropic moduli, the two dimensions, and plane
strain and plane stress are on [Elastic moduli](@ref man-elastic-moduli).

```@example par
using TensND, LinearAlgebra

ℂ = iso_stiffness_young_poisson(70.0, 0.33)   # an aluminum, GPa
bulk_shear(ℂ), young_poisson(ℂ)
```

## In two dimensions

Every builder of the isotropic moduli takes the dimension as a first argument
`Val(2)`, and every reader the keyword `convention`: `:plane_strain`, the
default, reads ``(E, \nu)`` of the material of space in plane strain,
`:intrinsic` the moduli ``(E_2, \nu_2)`` of the medium of the plane
([Elastic moduli](@ref man-moduli-2d)).

```@example par
ℂ₂ = iso_stiffness_young_poisson(Val(2), 70.0, 0.33)
bulk_shear(ℂ₂), young_poisson(ℂ₂), young_poisson(ℂ₂; convention = :intrinsic)
```

An orthotropic tensor of the plane is given by its four constants
``(C_{11}, C_{22}, C_{12}, C_{66})`` in its frame, which `get_data` reads back,
and its Kelvin–Mandel matrix in that frame by `mandel_material`; a square one
by ``(C_{11}, C_{12}, C_{44})``, read back by `cubic_constants`:

```@example par
𝕆₂ = TensOrtho(10.0, 6.0, 3.0, 2.0, Basis(0.3))
get_data(𝕆₂), cubic_constants(cubic_from_constants(10.0, 4.0, 2.0, Basis(0.3)))
```

The conventions of a transversely isotropic material that follow are those of
dimension 3. In a plane that contains its axis, the material is orthotropic,
of constants ``(n, k + m, l, p)`` in Hill's moduli
([Hill's moduli](@ref man-hill-moduli)), the first axis of the plane being the
axis of the material:

```@example par
get_data(plane_strain_stiffness(ti_from_constants(19.0, 5.0, 11.0, 12.0, 4.0, [0.0, 0.0, 1.0]; form = :hill), Val(:axial)))
```

See [Two dimensions, and any dimension](@ref man-dimensions) for the whole of
what is built in the plane.

## The five TI conventions

| `form` | Parameters | Builds |
| :--- | :--- | :--- |
| `:walpole` | ``\ell_1,\ell_2,\ell_3,\ell_5,\ell_6`` | stiffness *or* compliance |
| `:components` (default) | ``C_{1111},C_{1122},C_{1133},C_{3333},C_{2323}`` | stiffness *or* compliance |
| `:hill` | ``k,l,m,n,p`` | stiffness *or* compliance |
| `:engineering` | ``E_1,E_3,\nu_{12},\nu_{31},G_{31}`` | **compliance** |
| `:hoenig` | ``E,\nu_1,\nu_2,H,\Gamma`` | **compliance** |

[`ti_from_constants`](@ref)`(c₁, …, c₅, n; form)` builds the tensor, with the
symmetry axis as its last argument, and returns a `TensTI{4,T,5}`;
[`ti_constants`](@ref)`(t; form)` reads the five constants back.

```@example par
n = [0.0, 0.0, 1.0]
𝕊 = ti_from_constants(9.0, 140.0, 0.40, 0.30, 4.6, n; form = :engineering)   # a carbon/epoxy ply, GPa
ℂ = inv(𝕊)
ti_constants(𝕊; form = :engineering)
```

## Converting between them

Conversion is a **round trip through the tensor**: build with one constructor,
read with another extractor. Because all five produce the same object, the
conversion is exact rather than fitted.

```@example par
ti_constants(𝕊)          # the same compliance, in component form
```

```@example par
ti_constants(𝕊; form = :hoenig)   # and in dimensionless Hoenig ratios
```

``H`` is the axial-to-transverse modulus ratio and ``\Gamma`` the shear
anisotropy; an isotropic material sits at ``H=\Gamma=1`` with ``\nu_1=\nu_2``.

## [Hill's moduli](@id man-hill-moduli)

Hill's moduli ``(k, l, m, n, p)`` are those of fiber composites: for a
stiffness of axis ``\underline{e}_3``, the plane-strain bulk modulus
``k = \tfrac12(C_{1111}+C_{1122})``, the coupling ``l = C_{1133}``, the shear
modulus in the plane of isotropy ``m = \tfrac12(C_{1111}-C_{1122})``, the axial
modulus ``n = C_{3333}`` and the shear modulus along the axis ``p = C_{2323}``
([Hill form](@ref th-ti-hill)). `form = :hill` takes them in that order, and
[`hill_moduli`](@ref) reads them back by name:

```@example par
ℍ = ti_from_constants(19.0, 5.0, 11.0, 12.0, 4.0, n; form = :hill)    # k, l, m, n, p
hill_moduli(ℍ)
```

```@example par
h = hill_moduli(ℂ)          # of the carbon/epoxy ply above
h.k, h.p                    # its plane-strain bulk modulus and its axial shear modulus
```

The other forms of the same tensor are one call away, in both directions:

```@example par
ti_constants(ℍ; form = :components)        # (k + m, k − m, l, n, p)
```

```@example par
ti_constants(ℍ; form = :engineering, role = :stiffness)   # E₁, E₃, ν₁₂, ν₃₁, G₃₁ of the stiffness ℍ
```

```@example par
E₁, E₃, ν₁₂, ν₃₁, G₃₁ = 9.0, 140.0, 0.40, 0.30, 4.6
hill_moduli(ti_from_constants(E₁, E₃, ν₁₂, ν₃₁, G₃₁, n; form = :engineering, role = :stiffness))
```

Like the components, Hill's moduli are those of the tensor given: `role` has no
effect on them, and the moduli of a compliance are those of its own components.
An isotropic tensor has Hill's moduli about any axis, and symbols go through:

```@example par
using SymPyPythonCall
κ, μ = symbols("kappa mu", positive = true)
hill_moduli(iso_stiffness(κ, μ))
```

```@example par
k, l, m, nₕ, p = symbols("k l m n p", positive = true)
ti_constants(ti_from_constants(k, l, m, nₕ, p, (0, 0, 1); form = :hill); form = :engineering, role = :stiffness)[2:4]   # E₃, ν₁₂, ν₃₁
```

A field transversely isotropic about the radius of a sphere or of a cylinder
takes the same forms: see [A tensor that turns with the point](@ref man-chart-fields).

### Six constants: without the major symmetry

A transversely isotropic tensor without the major symmetry, a concentration
tensor or a product of two stiffnesses, has a sixth constant, its two couplings
``C_{1133}`` and ``C_{3311}`` being different. Three forms take six constants:

| `form` | Parameters |
| :--- | :--- |
| `:hill` | ``k, l, m, n, p, l_t`` with ``l = C_{1133}`` and ``l_t = C_{3311}`` |
| `:components` | ``C_{1111}, C_{1122}, C_{1133}, C_{3311}, C_{3333}, C_{2323}`` |
| `:walpole` | ``\ell_1, \ldots, \ell_6``, with ``\ell_3 = \sqrt2\,C_{3311}`` and ``\ell_4 = \sqrt2\,C_{1133}`` |

and `hill_moduli(t, Val(6))` reads the six Hill moduli back:

```@example par
ℙ = ℍ ⊡ ti_from_constants(8.0, 2.0, 3.0, 15.0, 6.0, n; form = :hill)     # a product of two stiffnesses
h₆ = hill_moduli(ℙ, Val(6))
```

```@example par
ti_from_constants(h₆..., n; form = :hill) ≈ ℙ
```

!!! warning "Which coupling is which"
    A tensor has the major symmetry when ``C_{3311} = C_{1133}``, that is
    ``l_t = l``: a stiffness or a compliance has it, a product of two of them
    or a concentration tensor in general not. `hill_moduli(t)` returns the five
    moduli of a major-symmetric tensor and refuses another one;
    `hill_moduli(t, Val(6))` returns ``l = C_{1133}``, the transverse stress
    due to an axial strain, and ``l_t = C_{3311}``, the axial stress due to a
    transverse strain, in this order. Exchanging the two builds the transpose
    of the tensor, without any error.

## Stiffness or compliance: `role`

The engineering and Hoenig constants are defined on the **compliance**:
`ti_from_constants(…; form = :engineering)` and `form = :hoenig` build
``\mathbb{S}``, and `ti_constants(t; form)` reads `t` as one. The keyword
`role = :stiffness` makes them build, or read, the **stiffness** of the same
material instead. It is never computed through an inverse: each constant has
its closed form on the stiffness components, with
``D = C_{1111}C_{3333} - C_{1133}^2``,

```math
E_1 = \frac{(C_{1111} - C_{1122})(C_{1111}C_{3333} + C_{1122}C_{3333} - 2C_{1133}^2)}{D},
\quad
E_3 = C_{3333} - \frac{2C_{1133}^2}{C_{1111} + C_{1122}},
\quad
\nu_{12} = \frac{C_{1122}C_{3333} - C_{1133}^2}{D},
\quad
\nu_{31} = \frac{C_{1133}}{C_{1111} + C_{1122}},
\quad
G_{31} = C_{2323},
```

and the stiffness has its closed form in the constants, with
``\Delta = E_3(1 - \nu_{12}) - 2E_1\nu_{31}^2``:

```math
C_{1111} = \frac{E_1(E_3 - E_1\nu_{31}^2)}{(1 + \nu_{12})\Delta},
\quad
C_{1122} = \frac{E_1(E_1\nu_{31}^2 + E_3\nu_{12})}{(1 + \nu_{12})\Delta},
\quad
C_{1133} = \frac{E_1E_3\nu_{31}}{\Delta},
\quad
C_{3333} = \frac{E_3^2(1 - \nu_{12})}{\Delta},
\quad
C_{2323} = G_{31}.
```

Hoenig's parameters are those of ``E_1 = E``, ``E_3 = EH``, ``\nu_{12} = \nu_1``,
``\nu_{31} = H\nu_2`` and ``G_{31} = E\Gamma/(2(1 + \nu_1))``.

```@example par
ℂ = ti_from_constants(9.0, 140.0, 0.40, 0.30, 4.6, n; form = :engineering, role = :stiffness)
ℂ ≈ inv(𝕊), ti_constants(ℂ; form = :engineering, role = :stiffness)
```

`role` is `:compliance` by default, as these forms have always been read, and
has no effect on `form = :walpole`, `:components` and `:hill`, whose five numbers
are those of the tensor given, whatever it is.

!!! warning "A role is the caller's to give"
    Nothing in a `TensTI` tells a stiffness from a compliance. Read with the
    wrong `role`, a tensor gives the constants of another material, and a
    stiffness entered by its engineering constants without `role = :stiffness`
    is built as a compliance, without any error.
    [Symbolic tensors of a material symmetry](@ref man-symbolic-materials)
    shows each case on symbols, with the matrix that is built.

## Walpole coefficients

A major-symmetric tensor has five independent Walpole coefficients,
``(\ell_1,\ell_2,\ell_3,\ell_5,\ell_6)``: `form = :walpole` builds it from them and
reads them back, as the other forms do with their constants.

```@example par
ti_constants(ℂ; form = :walpole)
```

```@example par
ti_from_constants(ti_constants(ℂ; form = :walpole)..., n; form = :walpole) == ℂ
```

[`walpole_coeffs`](@ref) reads the same coefficients in three lengths, whatever
the number on which the tensor is stored:

| Call | Returns | For |
| :--- | :--- | :--- |
| `walpole_coeffs(t, Val(5))` | ``(\ell_1,\ell_2,\ell_3,\ell_5,\ell_6)`` | a major-symmetric tensor (refused otherwise): `form = :walpole` |
| `walpole_coeffs(t)` | ``(\ell_1,\ldots,\ell_6)`` | any tensor of the class, ``\ell_3 = \ell_4`` when it is major-symmetric |
| `walpole_coeffs(t, Val(8))` | ``(\ell_1,\ldots,\ell_8)`` | the full axially invariant basis |

!!! warning "Five, six or eight: ask for the length you mean"
    - `walpole_coeffs(t)` returns **six** coefficients even for a
      major-symmetric tensor, ``\ell_3`` being repeated as ``\ell_4``. The form
      `:walpole` has **five**: its fourth number is ``\ell_5``, not ``\ell_4``.
      To feed `ti_from_constants(…; form = :walpole)` or `radial_ti`, read
      with `walpole_coeffs(t, Val(5))` or `ti_constants(t; form = :walpole)`,
      never with `walpole_coeffs(t)`.
    - The five are refused for a tensor that is not major-symmetric, such as
      the product `ℂ ⊡ 𝕊` of two tensors that do not commute: its six
      coefficients have ``\ell_3 \neq \ell_4``.
    - The coefficients are not moduli: ``\ell_1 = C_{3333}`` comes first, and
      ``\ell_2 = 2k``, ``\ell_3 = \sqrt2\,C_{1133}``, ``\ell_5 = 2m``,
      ``\ell_6 = 2p`` carry the factors of the normalized basis.

```@example par
walpole_coeffs(ℂ), walpole_coeffs(ℂ, Val(5))
```

The six are the coefficients on which the products of the class are written: a
product of two major-symmetric tensors has ``\ell_3 \neq \ell_4`` in general, and
the ``2\times2`` block ``\begin{pmatrix}\ell_1 & \ell_3\\ \ell_4 & \ell_2\end{pmatrix}``
multiplies as a matrix. Stiffness and compliance are inverse in the synthetic algebra, so their
``2\times2`` blocks are inverse matrices and their ``\ell_5,\ell_6`` reciprocal:

```@example par
ℓC, ℓS = walpole_coeffs(ℂ), walpole_coeffs(𝕊)
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

For a matrix in an arbitrary frame, use [`project`](@ref) instead — it handles
the rotation. See [Projection](@ref man-projection).
