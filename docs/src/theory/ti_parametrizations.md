# [Transversely isotropic parametrizations](@id th-ti-parametrizations)

!!! tip "See also"
    **Theory**: [The Walpole basis](@ref th-walpole)

    **Manual**: [Structured tensors](@ref man-structured) · [Parametrizations and conversions](@ref man-parametrizations) · [Symbolic tensors of a material symmetry](@ref man-symbolic-materials)

    **Tutorials**: [The Walpole basis end to end](@ref tut-walpole) · [Acoustic and Hill tensors](@ref tut-hill-tensors)

A transversely isotropic tensor with major symmetry has five degrees of freedom.
Which five numbers one uses is a matter of discipline, not of mathematics, and
`TensND` supports five conventions: the Walpole coefficients, the components,
Hill's moduli, the engineering constants and Hoenig's parameters. All of them build a
`TensTI{4,T,5}` and all are interconvertible through the Walpole coefficients
``\ell_i`` of [The Walpole basis](@ref th-walpole).

Every matrix below is the Kelvin–Mandel matrix in the frame where
``\underline{n}=\underline{e}_3``, with the ordering and ``\sqrt2`` factors of
[Kelvin–Mandel representation](@ref th-kelvin-mandel).

## Summary

[`ti_from_constants`](@ref) builds the tensor from five constants and
[`ti_constants`](@ref) reads them back, both according to the keyword `form`:

| `form` | Parameters | They are those of |
| :----- | :--------- | :---------------- |
| `:walpole` | ``\ell_1,\ell_2,\ell_3,\ell_5,\ell_6`` | the tensor given, stiffness or compliance |
| `:components` | ``C_{1111},C_{1122},C_{1133},C_{3333},C_{2323}`` | the tensor given, stiffness or compliance |
| `:hill` | ``k,l,m,n,p`` | the tensor given, stiffness or compliance |
| `:engineering` | ``E_1,E_3,\nu_{12},\nu_{31},G_{31}`` | the material: its compliance, or its stiffness with `role = :stiffness` |
| `:hoenig` | ``E,\nu_1,\nu_2,H,\Gamma`` | the material: its compliance, or its stiffness with `role = :stiffness` |

The first three forms are linear in the components: they describe a tensor,
whatever it stands for. The last two are constants of a material, which its
compliance shows directly and its stiffness through closed forms, never through
an inverse. [`hill_moduli`](@ref) reads Hill's moduli by name, and
[`walpole_coeffs`](@ref) the Walpole coefficients: five with `Val(5)`, those of
`form = :walpole`, or six by default, which also describe a tensor without
major symmetry ([The Walpole basis](@ref th-walpole)).

!!! warning "What to check before giving or reading five constants"
    - **The order of each form is its own.** Hill's moduli are
      ``(k, l, m, n, p)``, the components
      ``(C_{1111}, C_{1122}, C_{1133}, C_{3333}, C_{2323})`` and the Walpole
      coefficients ``(\ell_1, \ell_2, \ell_3, \ell_5, \ell_6)``, which start
      with the axial one: ``\ell_1 = n = C_{3333}``. Five numbers given in the
      order of another form build a valid tensor, of another material, without
      any error.
    - **The Walpole coefficients carry factors**:
      ``\ell_2 = 2k``, ``\ell_3 = \sqrt2\,l``, ``\ell_5 = 2m``, ``\ell_6 = 2p``.
    - **Five Walpole coefficients, not six.** `walpole_coeffs(t)` returns six,
      ``\ell_4 = \ell_3`` included; `form = :walpole` takes and returns five,
      without ``\ell_4``.
    - **`role` matters for two forms only.** The engineering constants and
      Hoenig's parameters build a *compliance* unless `role = :stiffness` is
      given; the Walpole coefficients, the components and Hill's moduli are
      those of the tensor given, and ignore `role`.
    - **A tensor without major symmetry has no five constants**, in any form:
      it is refused.

## Component form

The five independent components of a minor- and major-symmetric TI tensor:

```math
\mathrm{Mat}(\mathbb{C})=
\begin{pmatrix}
C_{1111} & C_{1122} & C_{1133} & 0 & 0 & 0\\
C_{1122} & C_{1111} & C_{1133} & 0 & 0 & 0\\
C_{1133} & C_{1133} & C_{3333} & 0 & 0 & 0\\
0&0&0& 2\,C_{2323} & 0 & 0\\
0&0&0& 0 & 2\,C_{2323} & 0\\
0&0&0& 0 & 0 & C_{1111}-C_{1122}
\end{pmatrix}
```

The ``(6,6)`` entry is *not* independent: in-plane isotropy forces
``2\,C_{1212}=C_{1111}-C_{1122}``, which is what reduces nine orthotropic
constants to five.

### Relation to the Walpole coefficients

```math
\ell_1=C_{3333},\qquad
\ell_2=C_{1111}+C_{1122},\qquad
\ell_3=\ell_4=\sqrt2\,C_{1133},
```

```math
\ell_5=C_{1111}-C_{1122},\qquad
\ell_6=2\,C_{2323}.
```

Inverted:

```math
C_{1111}=\tfrac{\ell_2+\ell_5}{2},\quad
C_{1122}=\tfrac{\ell_2-\ell_5}{2},\quad
C_{1133}=\tfrac{\ell_3}{\sqrt2},\quad
C_{3333}=\ell_1,\quad
C_{2323}=\tfrac{\ell_6}{2}.
```

These are exactly `_build_TI_KM` and its inverse, exported as
[`mandel_from_ti_params`](@ref) and [`ti_params_from_mandel`](@ref).

## [Hill form](@id th-ti-hill)

[hill1964](@citet) writes the law of a fiber-reinforced material on the mean
and the difference of the transverse components. With ``\underline{e}_3`` the
axis, a transversely isotropic stiffness relates the stress to the strain by

```math
\tfrac12(\sigma_{11}+\sigma_{22}) = k\,(\varepsilon_{11}+\varepsilon_{22}) + l\,\varepsilon_{33},
\qquad
\sigma_{33} = l\,(\varepsilon_{11}+\varepsilon_{22}) + n\,\varepsilon_{33},
```

```math
\sigma_{11}-\sigma_{22} = 2m\,(\varepsilon_{11}-\varepsilon_{22}),\qquad
\sigma_{12} = 2m\,\varepsilon_{12},\qquad
\sigma_{13} = 2p\,\varepsilon_{13},\qquad
\sigma_{23} = 2p\,\varepsilon_{23} .
```

Each modulus is the stiffness of one elementary loading:

| Modulus | Components | Loading that isolates it |
| :------ | :--------- | :----------------------- |
| ``k`` | ``\tfrac12(C_{1111}+C_{1122})`` | equal strains in the plane of isotropy, none along the axis: the plane-strain bulk modulus |
| ``l`` | ``C_{1133}`` | the coupling between the plane and the axis |
| ``m`` | ``\tfrac12(C_{1111}-C_{1122}) = C_{1212}`` | shear in the plane of isotropy |
| ``n`` | ``C_{3333}`` | strain along the axis, none across it |
| ``p`` | ``C_{2323} = C_{1313}`` | shear along the axis |

The Kelvin–Mandel matrix of the component form is then

```math
\mathrm{Mat}(\mathbb{C})=
\begin{pmatrix}
k+m & k-m & l & 0 & 0 & 0\\
k-m & k+m & l & 0 & 0 & 0\\
l & l & n & 0 & 0 & 0\\
0&0&0& 2p & 0 & 0\\
0&0&0& 0 & 2p & 0\\
0&0&0& 0 & 0 & 2m
\end{pmatrix} .
```

### Relation to the Walpole coefficients

Hill's moduli are the Walpole coefficients up to a factor, since both
decompose the tensor on the plane, the axis and their coupling:

```math
\ell_1 = n,\qquad \ell_2 = 2k,\qquad \ell_3=\ell_4 = \sqrt2\,l,\qquad
\ell_5 = 2m,\qquad \ell_6 = 2p .
```

The algebra of the Walpole basis thus acts on ``(k, l, m, n, p)`` directly:
``m`` and ``p`` are inverted alone, and ``(k, l, n)`` through the matrix
``\begin{pmatrix} n & \sqrt2\,l\\ \sqrt2\,l & 2k\end{pmatrix}``, of determinant
``2(kn - l^2)``. The tensor is positive definite when

```math
m > 0,\qquad p > 0,\qquad k > 0,\qquad k\,n > l^2 ,
```

and its inverse has, on its own components, the "moduli"

```math
k' = \frac{n}{4(kn-l^2)},\qquad l' = \frac{-l}{2(kn-l^2)},\qquad
n' = \frac{k}{kn-l^2},\qquad m' = \frac{1}{4m},\qquad p' = \frac{1}{4p} .
```

### Relation to the engineering constants

The compliance above is that of the engineering form of the next section, read
term by term: ``1/E_3 = n'``, ``-\nu_{31}/E_3 = l'``,
``1/E_1 = k'+m'`` and ``-\nu_{12}/E_1 = k'-m'``. Hence, for a **stiffness** of
Hill's moduli ``(k, l, m, n, p)``,

```math
E_3 = n - \frac{l^2}{k},\qquad \nu_{31} = \frac{l}{2k},\qquad G_{31} = p,\qquad
G_{12} = m ,
```

```math
\frac{4}{E_1} = \frac1m + \frac1k + \frac{4\nu_{31}^2}{E_3},\qquad
\nu_{12} = \frac{(k-m)\,n - l^2}{(k+m)\,n - l^2} ,
```

and conversely

```math
m = \frac{E_1}{2(1+\nu_{12})},\qquad p = G_{31},\qquad
k = \frac{E_1E_3}{2\bigl(E_3(1-\nu_{12}) - 2E_1\nu_{31}^2\bigr)},\qquad
l = 2k\,\nu_{31},\qquad n = E_3 + 4k\,\nu_{31}^2 .
```

Hoenig's parameters follow by ``E_1 = E``, ``E_3 = EH``,
``\nu_{12} = \nu_1``, ``\nu_{31} = H\nu_2`` and
``G_{31} = E\Gamma/(2(1+\nu_1))``; in particular ``\Gamma = p/m``.

These relations are obtained by building the tensor in one form and reading it
in another, on symbols:

```@example hill
using TensND, SymPyPythonCall

k, l, m, n, p = symbols("k l m n p", positive = true)
ℂ = ti_from_constants(k, l, m, n, p, (0, 0, 1); form = :hill)
E₁, E₃, ν₁₂, ν₃₁, G₃₁ = simplify.(ti_constants(ℂ; form = :engineering, role = :stiffness))
(E₃, ν₃₁, G₃₁)
```

```@example hill
simplify(4 / E₁ - (1 / m + 1 / k + 4 * ν₃₁^2 / E₃)), simplify(ν₁₂ - ((k - m) * n - l^2) / ((k + m) * n - l^2))
```

```@example hill
simplify.(values(hill_moduli(inv(ℂ))))        # (k′, l′, m′, n′, p′) of the compliance
```

```@example hill
E, F, ν, ν′, G = symbols("E_1 E_3 nu_12 nu_31 G_31", positive = true)
h = hill_moduli(ti_from_constants(E, F, ν, ν′, G, (0, 0, 1); form = :engineering, role = :stiffness))
simplify(h.l - 2 * h.k * ν′), simplify(h.n - (F + 4 * h.k * ν′^2)), simplify(h.m - E / (2 * (1 + ν)))
```

### The isotropic case

An isotropic stiffness ``3\kappa\,\mathbb{J} + 2\mu\,\mathbb{K}`` is transversely
isotropic about any axis, with

```math
k = \kappa + \frac{\mu}{3},\qquad l = \kappa - \frac{2\mu}{3},\qquad
m = p = \mu,\qquad n = \kappa + \frac{4\mu}{3} ,
```

so that isotropy is the plane ``m = p``, ``k - l = m``, ``n - l = 2m`` of the
space of Hill's moduli.

```@example hill
κ, μ = symbols("kappa mu", positive = true)
hill_moduli(iso_stiffness(κ, μ))
```

## Engineering form

The convention of composite mechanics, for the **compliance**:

- ``E_1`` — transverse Young's modulus, in the isotropy plane;
- ``E_3`` — axial Young's modulus, along ``\underline{n}``;
- ``\nu_{12}`` — in-plane Poisson's ratio;
- ``\nu_{31}`` — axial–transverse Poisson's ratio, with the reciprocity
  ``\nu_{31}/E_3=\nu_{13}/E_1``;
- ``G_{31}`` — axial shear modulus.

```math
\mathrm{Mat}(\mathbb{S})=
\begin{pmatrix}
\frac{1}{E_1} & \frac{-\nu_{12}}{E_1} & \frac{-\nu_{31}}{E_3} & 0&0&0\\
\frac{-\nu_{12}}{E_1} & \frac{1}{E_1} & \frac{-\nu_{31}}{E_3} & 0&0&0\\
\frac{-\nu_{31}}{E_3} & \frac{-\nu_{31}}{E_3} & \frac{1}{E_3} & 0&0&0\\
0&0&0& \frac{1}{2G_{31}} & 0 & 0\\
0&0&0& 0 & \frac{1}{2G_{31}} & 0\\
0&0&0& 0&0& \frac{1+\nu_{12}}{E_1}
\end{pmatrix}
```

The ``(6,6)`` entry encodes the in-plane shear modulus
``G_{12}=E_1/\bigl(2(1+\nu_{12})\bigr)``, which is *determined* by ``E_1`` and
``\nu_{12}`` — the isotropy of the transverse plane again.

!!! note "Kelvin–Mandel, not Voigt"
    The shear entries read ``1/(2G_{31})`` and ``(1+\nu_{12})/E_1``, where an
    engineering (Voigt) compliance matrix would show ``1/G_{31}`` and
    ``2(1+\nu_{12})/E_1``. The factor two is the Kelvin–Mandel convention, and
    it is what keeps ``\mathrm{Mat}(\mathbb{S})=\mathrm{Mat}(\mathbb{C})^{-1}``
    an ordinary matrix inverse.

## Hoenig form

A dimensionless parametrization introduced for crack problems in an anisotropic
medium [hoenig1978](@cite), convenient when anisotropy *ratios* matter
independently of the overall stiffness scale:

| Parameter | Definition |
| :-------- | :--------- |
| ``E`` | transverse Young's modulus, ``=1/S_{1111}`` |
| ``\nu_1`` | in-plane Poisson's ratio, ``=-E\,S_{1122}`` |
| ``\nu_2`` | axial–transverse Poisson's ratio, ``=-E\,S_{1133}`` |
| ``H`` | axial-to-transverse modulus ratio, ``=1/(E\,S_{3333})`` |
| ``\Gamma`` | shear anisotropy parameter, ``=(1+\nu_1)/(2E\,S_{2323})`` |

so that the compliance components are

```math
S_{1111}=\frac{1}{E},\quad
S_{1122}=\frac{-\nu_1}{E},\quad
S_{1133}=\frac{-\nu_2}{E},\quad
S_{3333}=\frac{1}{E\,H},\quad
S_{2323}=\frac{1+\nu_1}{2E\,\Gamma}.
```

Isotropy is the point ``H=\Gamma=1`` with ``\nu_1=\nu_2``: ``H`` measures the
departure of the axial stiffness from the transverse one, and ``\Gamma`` that of
the axial shear modulus from the in-plane one. The companion paper
[hoenig1979](@citet) applies the parametrization to the effective moduli of a
non-randomly cracked body.

## Choosing between them

| Use | `form` |
| :-- | :----- |
| products, inverses and contractions written on the coefficients of the class | `:walpole` |
| a stiffness known by its components, or symbolic work | `:components` |
| fiber composites, layered cylinders and spheres, plane-strain problems | `:hill` |
| material data sheets | `:engineering` |
| anisotropy ratios, crack and inclusion problems | `:hoenig` |

Conversion between any two is a round trip through the tensor itself: build with
one form, read with another. Since all five produce the same
`TensTI{4,T,5}`, the conversion is exact rather than fitted — unlike the
*projection* of a general tensor onto the class, which is the subject of
[Projection onto a symmetry class](@ref th-projection).
