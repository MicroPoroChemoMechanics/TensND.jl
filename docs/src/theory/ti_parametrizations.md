# [Transversely isotropic parametrizations](@id th-ti-parametrizations)

A transversely isotropic tensor with major symmetry has five degrees of freedom.
Which five numbers one uses is a matter of discipline, not of mathematics, and
`TensND` supports three conventions. All of them build a `TensTI{4,T,5}` and all
are interconvertible through the Walpole coefficients ``\ell_i`` of
[The Walpole basis](@ref th-walpole).

Every matrix below is the Kelvin–Mandel matrix in the frame where
``\underline{n}=\underline{e}_3``, with the ordering and ``\sqrt2`` factors of
[Kelvin–Mandel representation](@ref th-kelvin-mandel).

## Summary

| Constructor | Extractor | Parameters | Builds a |
| :---------- | :-------- | :--------- | :------- |
| [`tens_TI`](@ref) | [`arg_TI`](@ref) | ``C_{1111},C_{1122},C_{1133},C_{3333},C_{2323}`` | stiffness *or* compliance |
| [`tens_TI_eng`](@ref) | [`arg_TI_eng`](@ref) | ``E_1,E_3,\nu_{12},\nu_{31},G_{31}`` | **compliance** |
| [`tens_TI_Hoenig`](@ref) | [`arg_TI_Hoenig`](@ref) | ``E,\nu_1,\nu_2,H,\Gamma`` | **compliance** |

To obtain a stiffness from either compliance form, invert: `inv(tens_TI_eng(…))`.
Inversion is exact and stays in the class, by the synthetic rule of
[The Walpole basis](@ref th-walpole).

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
[`KM_from_ti_params`](@ref) and [`ti_params_from_KM`](@ref).

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
the axial shear modulus from the in-plane one. The same author's companion paper
[hoenig1979](@cite) applies the parametrization to the effective moduli of a
non-randomly cracked body.

## Choosing between them

| Use | Form |
| :-- | :--- |
| a stiffness known by its components, or symbolic work | `tens_TI` |
| material data sheets, composites | `tens_TI_eng` |
| anisotropy ratios, crack and inclusion problems | `tens_TI_Hoenig` |

Conversion between any two is a round trip through the tensor itself: build with
one constructor, extract with another extractor. Since all three produce the same
`TensTI{4,T,5}`, the conversion is exact rather than fitted — unlike the
*projection* of a general tensor onto the class, which is the subject of
[Projection onto a symmetry class](@ref th-projection).
