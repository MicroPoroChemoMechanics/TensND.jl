# [Symbolic tensors of a material symmetry](@id man-symbolic-materials)

!!! tip "See also"
    **Theory**: [The Walpole basis](@ref th-walpole) · [Transversely isotropic parametrizations](@ref th-ti-parametrizations) · [Isotropic tensors](@ref th-isotropic)

    **Manual**: [Structured tensors](@ref man-structured) · [Parametrizations and conversions](@ref man-parametrizations) · [Elastic moduli](@ref man-elastic-moduli)

    **Tutorials**: [The Walpole basis end to end](@ref tut-walpole) · [Acoustic and Hill tensors](@ref tut-hill-tensors)

To test a formula, one needs a stiffness or a compliance of a given symmetry
whose constants are symbols. This page builds one for each class, in each
notation, and shows it as its Kelvin–Mandel matrix, so that what a constructor
builds is seen and not assumed: which constants it takes, in which order, and
whether the result is a stiffness or a compliance. Each tensor built is then
read back by the function that extracts its constants, which must return the
symbols it was built with: it is the check to make before trusting a
construction.

```@example symmat
using TensND, SymPyPythonCall, LinearAlgebra
show_matrix(t) = simplify.(mandel(t))          # the 6×6 Kelvin–Mandel matrix, in the canonical basis
nothing # hide
```

In that matrix the three shear rows and columns carry the factors ``\sqrt2``
and ``2`` of the Kelvin–Mandel convention: its entry ``(4,4)`` is
``2\,C_{2323}``, and the matrix of the inverse tensor is the inverse matrix.

## Isotropic

A stiffness ``\mathbb{C} = 3\kappa\,\mathbb{J} + 2\mu\,\mathbb{K}`` from its bulk
and shear moduli, or from any other pair
([Elastic moduli](@ref man-elastic-moduli)):

```@example symmat
κ, μ = symbols("kappa mu", positive = true)
ℂ = iso_stiffness(κ, μ)
show_matrix(ℂ)
```

Its compliance is its inverse, exact in the class. `iso_compliance` writes it
directly from the same moduli:

```@example symmat
𝕊 = iso_compliance(κ, μ)
show_matrix(𝕊)
```

```@example symmat
inv(ℂ) == 𝕊, simplify.(mandel(ℂ) * mandel(𝕊))
```

Read back, each tensor gives the moduli it was built with, provided it is read
for what it is: `role = :stiffness`, the default of these readers, or
`role = :compliance`.

```@example symmat
bulk_shear(ℂ), simplify.(bulk_shear(𝕊; role = :compliance))
```

```@example symmat
simplify.(young_poisson(ℂ)), simplify.(young_poisson(𝕊; role = :compliance))
```

From Young's modulus and Poisson's ratio, the stiffness and the compliance of
the same material:

```@example symmat
E, ν = symbols("E nu", positive = true)
𝕊ᴱ = iso_compliance_young_poisson(E, ν)
show_matrix(𝕊ᴱ)
```

```@example symmat
simplify.(young_poisson(𝕊ᴱ; role = :compliance)), simplify.(young_poisson(iso_stiffness_young_poisson(E, ν)))
```

## Cubic

Three constants in the frame of the cube, here the canonical one:

```@example symmat
C₁₁, C₁₂, C₄₄ = symbols("C_11 C_12 C_44", positive = true)
cube = cubic_from_constants(C₁₁, C₁₂, C₄₄, CanonicalBasis{3, Sym}())
show_matrix(cube)
```

```@example symmat
cubic_constants(cube)
```

The compliance is the inverse, a tensor of the same class, whose three
constants are read back:

```@example symmat
simplify.(cubic_constants(inv(cube)))
```

## Transversely isotropic

The axis is given last, here ``\underline{e}_3``. Five notations build the
same kind of tensor, a `TensTI{4}` stored on five coefficients; the keyword
`form` says which one the five numbers are in.

### The constants of the tensor itself: components, Hill's moduli, Walpole coefficients

These three notations are linear in the components. They describe the tensor
that is built, whatever it stands for: a stiffness if the numbers are those of
a stiffness, a compliance if they are those of a compliance. The keyword
`role` has no effect on them.

```@example symmat
e₃ = (0, 0, 1)
c₁₁, c₁₂, c₁₃, c₃₃, c₄₄ = symbols("C_1111 C_1122 C_1133 C_3333 C_2323", positive = true)
show_matrix(ti_from_constants(c₁₁, c₁₂, c₁₃, c₃₃, c₄₄, e₃))                 # form = :components, the default
```

```@example symmat
k, l, m, n, p = symbols("k l m n p", positive = true)
ℍ = ti_from_constants(k, l, m, n, p, e₃; form = :hill)                      # Hill's moduli, in this order
show_matrix(ℍ)
```

```@example symmat
ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆ = symbols("l_1 l_2 l_3 l_5 l_6", positive = true)
𝕎 = ti_from_constants(ℓ₁, ℓ₂, ℓ₃, ℓ₅, ℓ₆, e₃; form = :walpole)             # five Walpole coefficients
show_matrix(𝕎)
```

The three matrices show where each constant goes. The Walpole coefficients
start with the axial one, ``\ell_1`` in the entry ``(3,3)``, and carry
factors: ``\ell_2 = C_{1111} + C_{1122}``, ``\ell_3 = \sqrt2\,C_{1133}``,
``\ell_5 = C_{1111} - C_{1122}``, ``\ell_6 = 2\,C_{2323}``.

Each one is read back in its own notation, which returns the five symbols
given:

```@example symmat
ti_constants(ti_from_constants(c₁₁, c₁₂, c₁₃, c₃₃, c₄₄, e₃))
```

```@example symmat
ti_constants(ℍ; form = :hill), hill_moduli(ℍ)
```

```@example symmat
ti_constants(𝕎; form = :walpole)
```

and in any other one, which converts a notation into another:

```@example symmat
ti_constants(ℍ; form = :walpole), ti_constants(𝕎; form = :hill)
```

!!! warning "Five Walpole coefficients, six, or eight"
    `ti_constants(t; form = :walpole)` and `walpole_coeffs(t, Val(5))` return
    the **five** coefficients ``(\ell_1, \ell_2, \ell_3, \ell_5, \ell_6)``.
    `walpole_coeffs(t)` returns **six**, ``\ell_4 = \ell_3`` in fourth place:

    ```julia
    walpole_coeffs(𝕎)            # (l_1, l_2, l_3, l_3, l_5, l_6)
    walpole_coeffs(𝕎, Val(5))    # (l_1, l_2, l_3, l_5, l_6)
    ```

    The six are those of a tensor that need not be major-symmetric, such as a
    product; the five exist for a major-symmetric tensor only.

### The constants of the material: engineering constants, Hoenig's parameters

Young's moduli and Poisson's ratios are not components of a tensor: they are
constants of a **material**, which has a stiffness and a compliance. The
keyword `role` says which of the two is built from them, or read:

| Call | Builds |
| :--- | :--- |
| `ti_from_constants(E₁, E₃, ν₁₂, ν₃₁, G₃₁, n; form = :engineering)` | the **compliance** of the material (`role = :compliance`, the default) |
| `ti_from_constants(E₁, E₃, ν₁₂, ν₃₁, G₃₁, n; form = :engineering, role = :stiffness)` | its **stiffness** |
| `ti_constants(t; form = :engineering)` | reads `t` as a compliance |
| `ti_constants(t; form = :engineering, role = :stiffness)` | reads `t` as a stiffness |

and the same with `form = :hoenig`. The compliance shows the constants
directly:

```@example symmat
E₁, E₃, ν₁₂, ν₃₁, G₃₁ = symbols("E_1 E_3 nu_12 nu_31 G_31", positive = true)
𝕊ₑ = ti_from_constants(E₁, E₃, ν₁₂, ν₃₁, G₃₁, e₃; form = :engineering)       # a compliance
show_matrix(𝕊ₑ)
```

The stiffness of the same material is written in closed form, without an
inverse; the two tensors are inverse of each other:

```@example symmat
ℂₑ = ti_from_constants(E₁, E₃, ν₁₂, ν₃₁, G₃₁, e₃; form = :engineering, role = :stiffness)
show_matrix(ℂₑ)
```

```@example symmat
simplify.(mandel(ℂₑ) * mandel(𝕊ₑ))
```

Reading follows the same rule: each of the two tensors gives back the five
constants when it is read for what it is.

```@example symmat
simplify.(ti_constants(𝕊ₑ; form = :engineering))                           # the compliance, read as one
```

```@example symmat
simplify.(ti_constants(ℂₑ; form = :engineering, role = :stiffness))        # the stiffness, read as one
```

Hoenig's parameters go the same way, here on the stiffness:

```@example symmat
Eₕ, ν₁, ν₂, H, Γ = symbols("E nu_1 nu_2 H Gamma", positive = true)
ℂₕ = ti_from_constants(Eₕ, ν₁, ν₂, H, Γ, e₃; form = :hoenig, role = :stiffness)
simplify.(ti_constants(ℂₕ; form = :hoenig, role = :stiffness))
```

The engineering constants of the stiffness of
Hill's moduli ``(k, l, m, n, p)`` are asked with `role = :stiffness`:

```@example symmat
simplify.(ti_constants(ℍ; form = :engineering, role = :stiffness))        # E₁, E₃, ν₁₂, ν₃₁, G₃₁
```

!!! warning "Nothing in a tensor tells a stiffness from a compliance"
    - The default is `role = :compliance`. To enter a **stiffness** by its
      engineering constants or by Hoenig's parameters, `role = :stiffness`
      must be written; without it the tensor built is the compliance, a valid
      tensor of the class, and no error is raised.
    - Read with the wrong `role`, a tensor gives the constants of another
      material: below, the stiffness `ℂₑ` read as if it were a compliance.
    - `role` does nothing with `form = :components`, `:hill` or `:walpole`:
      those five numbers are those of the tensor given. The stiffness and the
      compliance of a material have different Hill's moduli, and one goes from
      the first to the second by `inv`.

```@example symmat
E₃_right = simplify(ti_constants(ℂₑ; form = :engineering, role = :stiffness)[2])
E₃_wrong = simplify(ti_constants(ℂₑ; form = :engineering)[2])               # ℂₑ taken for a compliance
E₃_right, E₃_wrong
```

The second value is ``1/C_{3333}``: the entry of the tensor that would be
``1/E_3`` if it were a compliance.

```@example symmat
simplify.(values(hill_moduli(inv(ℍ))))          # (k, l, m, n, p) of the compliance of ℍ
```

## Orthotropic

Nine constants in the frame of the material, in the order
``(C_{11}, C_{22}, C_{33}, C_{12}, C_{13}, C_{23}, C_{44}, C_{55}, C_{66})`` of
its Voigt indices:

```@example symmat
o = symbols("C_11 C_22 C_33 C_12 C_13 C_23 C_44 C_55 C_66", positive = true)
ortho = TensOrtho(o..., CanonicalBasis{3, Sym}())
show_matrix(ortho)
```

Its nine constants are those it stores:

```@example symmat
get_data(ortho)
```

## Another axis, another frame

The matrices above are simple because the axis is ``\underline{e}_3`` and the
frame the canonical one. A test should not rely on it: the same constants
about another axis give a full matrix, and an identity that holds for the first
must hold for the second.

```@example symmat
θ = symbols("theta", real = true)
tilted = ti_from_constants(k, l, m, n, p, (0, sin(θ), cos(θ)); form = :hill)
show_matrix(tilted)[1:3, 1:3]
```

The constants are read about the axis the tensor carries, so that they do not
depend on it:

```@example symmat
simplify.(values(hill_moduli(tilted)))
```

## From symbols to numbers, and to derivatives

A symbolic tensor is evaluated by substitution, and the same constructors take
numbers and dual numbers: an identity found on symbols is then a test on
numbers, and a derivative with respect to a constant is exact.

```@example symmat
numbers = (k => 19, l => 5, m => 11, n => 12, p => 4)
Float64.(N.(subs.(mandel(ℍ), numbers...))) ≈ mandel(ti_from_constants(19.0, 5.0, 11.0, 12.0, 4.0, (0.0, 0.0, 1.0); form = :hill))
```

```@example symmat
using ForwardDiff
axial_young(k) = ti_constants(ti_from_constants(k, 5.0, 11.0, 12.0, 4.0, (0.0, 0.0, 1.0); form = :hill); form = :engineering, role = :stiffness)[2]
ForwardDiff.derivative(axial_young, 19.0), 5.0^2 / 19.0^2        # E₃ = n − l²/k
```

## Summary

| Class | Symbolic stiffness | Its compliance | Constants read by |
| :--- | :--- | :--- | :--- |
| isotropic | `iso_stiffness(κ, μ)`, `iso_stiffness_young_poisson(E, ν)`, … | `iso_compliance(κ, μ)`, …, or `inv` | `bulk_shear`, `young_poisson`, … with `role` |
| cubic | `cubic_from_constants(C₁₁, C₁₂, C₄₄, frame)` | `inv` | `cubic_constants` |
| transversely isotropic | `ti_from_constants(c…, n; form)`, with `role = :stiffness` for `:engineering` and `:hoenig` | the same with `role = :compliance`, or `inv` | `ti_constants(t; form, role)`, `hill_moduli`, `walpole_coeffs` |
| orthotropic | `TensOrtho(C₁₁, …, C₆₆, frame)` | `inv` | `get_data` |

[Structured tensors](@ref man-structured) says what the products and the
inverses of these tensors return, and
[Parametrizations and conversions](@ref man-parametrizations) gives the
conversions on a numerical example. A field of such tensors whose axis turns
with the point is built by `radial_ti`, with the same `form` and `role`
([A tensor that turns with the point](@ref man-chart-fields)).
