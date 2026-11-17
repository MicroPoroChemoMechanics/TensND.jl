# [Solid harmonics, vector harmonics and the Lamé solutions](@id th-solid-harmonics)

This page defines the complex solid harmonics that expand a harmonic field
about a point, the vector harmonics and the displacements of an isotropic
elastic solid built on them, the formulas that carry each family from one
center to another, and their sums over a periodic array. The meridian parts of
the same harmonics in the spherical chart are those of
[Spheroidal harmonics and their re-expansion](@ref th-spheroidal-harmonics);
here they are written with their azimuthal factor and in the normalization in
which the translation formulas have no other factor than the harmonics
themselves [kushch2013](@cite).

## 1. Solid harmonics

Let ``(r, \theta, \varphi)`` be the spherical coordinates of ``\underline{x}``
about ``\underline{e}_3`` and ``P_n^m`` the Ferrers function of
[`legendre_P`](@ref), for which ``P_1^1(p) = -\sqrt{1-p^2}``. The surface
harmonic of degree ``n \ge 0`` and order ``0 \le m \le n`` is

```math
\chi_n^m(\theta, \varphi) = (-1)^m\,P_n^m(\cos\theta)\,e^{\mathrm{i}m\varphi} ,
```

so that ``\chi_1^1 = \sin\theta\,e^{\mathrm{i}\varphi}``, and the regular and
decaying solid harmonics are

```math
y_n^m(\underline{x}) = \frac{r^n}{(n+m)!}\,\chi_n^m , \qquad
Y_n^m(\underline{x}) = \frac{(n-m)!}{r^{n+1}}\,\chi_n^m ,
```

with ``y_n^{-m} = (-1)^m\,\overline{y_n^m}`` and the same for ``Y``. The first
ones are ``y_0^0 = 1``, ``y_1^0 = x_3``, ``y_1^1 = (x_1 + \mathrm{i}x_2)/2``,
``Y_0^0 = 1/r``, ``Y_1^0 = x_3/r^3``, ``Y_1^1 = (x_1 + \mathrm{i}x_2)/r^3``. Up to
the factors ``(-1)^m e^{\mathrm{i}m\varphi}/(n+m)!`` and
``(-1)^m (n-m)!\,e^{\mathrm{i}m\varphi}``, they are the harmonics
[`interior_harmonics`](@ref) and [`exterior_harmonics`](@ref) of the
[`SphericalChart`](@ref).

[`regular_solid_harmonics`](@ref) and [`irregular_solid_harmonics`](@ref)
compute all the harmonics of degree ``n \le N`` at a point by recurrences on
its Cartesian components, with ``w = x_1 + \mathrm{i}x_2``,

```math
y_m^m = \frac{w}{2m}\,y_{m-1}^{m-1} , \qquad
y_{n+1}^m = \frac{(2n+1)\,x_3\,y_n^m - r^2\,y_{n-1}^m}{(n-m+1)(n+m+1)} ,
```

```math
Y_m^m = \frac{(2m-1)\,w}{r^2}\,Y_{m-1}^{m-1} , \qquad
Y_{n+1}^m = \frac{(2n+1)\,x_3\,Y_n^m - (n-m)(n+m)\,Y_{n-1}^m}{r^2} ,
```

which involve neither an angle nor a factorial and hold at the origin for
``y``. The two families are related by
``Y_n^m = (n-m)!\,(n+m)!\,y_n^m/r^{2n+1}``; the factor
``\sqrt{(n-m)!\,(n+m)!}`` of [`solid_harmonic_norm`](@ref) brings the
coefficients of an expansion to the same order at every degree.

## 2. From one center to another

Away from its center, a decaying harmonic is a regular field. About a second
center, at ``\underline{R}`` from the first, it expands on the regular
harmonics [kushch2013](@cite):

```math
Y_{n'}^{m'}(\underline{r} + \underline{R}) =
  \sum_{n \ge 0}\,\sum_{|m| \le n} (-1)^{n+m}\,Y_{n'+n}^{m'-m}(\underline{R})\;y_n^m(\underline{r}) ,
\qquad \|\underline{r}\| < \|\underline{R}\| ,
```

a harmonic whose order exceeds its degree being zero. The series converges in
the ball about the second center that reaches the first one.

## 3. Vector harmonics

Three families of harmonic vectors are built on each family of solid
harmonics [kushch2013](@cite). The decaying ones are

```math
\underline{Y}^{(1)}_{nm} = \underline{\operatorname{grad}}\,Y_n^m , \qquad
\underline{Y}^{(2)}_{nm} = \frac{1}{n}\,\underline{\operatorname{grad}}\,Y_n^m\times\underline{x} ,
```

a gradient and a rotation field, both of zero divergence, and a third family
whose Cartesian components are harmonics of degree ``n-1``, such that

```math
\underline{\operatorname{curl}}\,\underline{Y}^{(3)}_{nm} = \underline{Y}^{(2)}_{nm} , \qquad
\operatorname{div}\,\underline{Y}^{(3)}_{nm} = -Y_n^m , \qquad
\underline{x}\cdot\underline{Y}^{(3)}_{nm} = \frac{r^2\,Y_n^m}{2n-1} .
```

The regular families are built in the same way on ``y_n^m``:
``\underline{y}^{(1)}_{nm} = \underline{\operatorname{grad}}\,y_n^m``,
``\underline{y}^{(2)}_{nm} = -\underline{\operatorname{grad}}\,y_n^m\times\underline{x}/(n+1)``,
and ``\underline{y}^{(3)}_{nm}`` with the same curl and divergence and
``\underline{x}\cdot\underline{y}^{(3)}_{nm} = -r^2\,y_n^m/(2n+3)``. In the complex
basis ``\underline{e}_\pm = (\underline{e}_1 \pm \mathrm{i}\,\underline{e}_2)/2``,
``\underline{e}_3``, the components of each are scalar solid harmonics of a
neighboring degree and order (Eqs. (4.14)–(4.15) of [kushch2013](@citet)), which
is how [`irregular_vector_harmonic`](@ref) and [`regular_vector_harmonic`](@ref)
read them from a table.

## 4. Solutions of the Lamé equation

The displacement ``\underline{u}`` of an isotropic solid of Poisson ratio ``\nu``
without body force solves

```math
\frac{2(1-\nu)}{1-2\nu}\,\underline{\operatorname{grad}}\,\operatorname{div}\,\underline{u}
  - \underline{\operatorname{curl}}\,\underline{\operatorname{curl}}\,\underline{u} = \underline{0} .
```

The first two families of vector harmonics are of zero divergence and their
curl is a gradient or zero: they solve it as they are, whatever ``\nu``. For
any harmonic vector ``\underline{B}``, the field of Papkovich and Neuber
``4(1-\nu)\,\underline{B} - \underline{\operatorname{grad}}\,(\underline{x}\cdot\underline{B})``
solves it too, hence the three decaying solutions

```math
\underline{U}^{(1)}_{nm} = \underline{Y}^{(1)}_{nm} , \qquad
\underline{U}^{(2)}_{nm} = \underline{Y}^{(2)}_{nm} , \qquad
\underline{U}^{(3)}_{nm} = 4(1-\nu)\,\underline{Y}^{(3)}_{nm}
  - \underline{\operatorname{grad}}\,\big(\underline{x}\cdot\underline{Y}^{(3)}_{nm}\big) ,
```

and the three regular ones ``\underline{u}^{(i)}_{nm}`` built in the same way on
``\underline{y}^{(i)}_{nm}``. The Poisson ratio enters through the factor
``4(1-\nu)`` alone. With the last relation of section 3,
``\underline{\operatorname{grad}}\,(\underline{x}\cdot\underline{Y}^{(3)}_{nm})
= (2\,\underline{x}\,Y_n^m + r^2\,\underline{Y}^{(1)}_{nm})/(2n-1)``, which
[`irregular_lame_solution`](@ref) and [`regular_lame_solution`](@ref) use.

The lowest degrees are the uniform strains and the rigid motions:
``\underline{u}^{(3)}_{00}`` is the dilatation, ``\underline{u}^{(1)}_{1m}`` the
translations, ``\underline{u}^{(2)}_{1m}`` the rotations and
``\underline{u}^{(1)}_{2m}`` the five uniform deviatoric strains. Among the
decaying ones, ``\underline{U}^{(1)}_{00}`` is the center of dilatation,
``\underline{U}^{(3)}_{1m}`` the point forces and ``\underline{U}^{(2)}_{1m}`` the
point torques.

## 5. Translation of the Lamé solutions

The decaying vector harmonics expand about a second center on the regular ones
of the same and of the lower families [kushch2013](@cite):

```math
\underline{Y}^{(i)}_{n'm'}(\underline{r} + \underline{R}) =
  \sum_{j \le i}\,\sum_{n,m} (-1)^{n+m}\,\alpha^{(i)(j)}_{n'n\,m'm}\;\underline{y}^{(j)}_{nm}(\underline{r}) ,
```

```math
\begin{aligned}
\alpha^{(1)(1)}_{n'n\,m'm} &= \alpha^{(2)(2)}_{n'n\,m'm} = \alpha^{(3)(3)}_{n'n\,m'm} = Y_{n'+n}^{m'-m} , \qquad
\alpha^{(2)(1)}_{n'n\,m'm} = -\alpha^{(3)(2)}_{n'n\,m'm} = \mathrm{i}\,\Big(\frac{m}{n} + \frac{m'}{n'}\Big)\,Y_{n'+n-1}^{m'-m} ,\\
\alpha^{(3)(1)}_{n'n\,m'm} &= \Big[1 + \frac{m\,m'}{n\,n'} - \frac{n'^2-m'^2}{n'(2n'-1)} - \frac{n^2-m^2}{n(2n-1)}\Big]\,Y_{n'+n-2}^{m'-m} ,
\end{aligned}
```

the harmonics being taken at ``\underline{R}``. The first two Lamé solutions
translate as these vectors do. For the third, the position from the first
center is ``\underline{x} = \underline{r} + \underline{R}``, and with
``\underline{B} = \underline{Y}^{(3)}_{n'm'}(\underline{r} + \underline{R})``,

```math
\underline{U}^{(3)}_{n'm'} = \big[4(1-\nu)\,\underline{B} - \underline{\operatorname{grad}}\,(\underline{r}\cdot\underline{B})\big]
  - \underline{\operatorname{grad}}\,(\underline{R}\cdot\underline{B}) .
```

``\underline{B}`` being expanded as above, the bracket is term by term a regular
Lamé solution about the second center: ``\underline{u}^{(3)}_{nm}`` for
``\underline{y}^{(3)}_{nm}``; ``4(1-\nu)\,\underline{u}^{(2)}_{nm}`` for
``\underline{y}^{(2)}_{nm}``, which is tangent to the spheres; and
``(4(1-\nu) - n)\,\underline{u}^{(1)}_{nm}`` for
``\underline{y}^{(1)}_{nm} = \underline{\operatorname{grad}}\,y_n^m``, since
``\underline{r}\cdot\underline{\operatorname{grad}}\,y_n^m = n\,y_n^m``. The last
term is the gradient of the harmonic function ``\underline{R}\cdot\underline{B}``,
whose components translate by section 2: with ``H_{n'n\,m'm}`` its coefficient
on ``(-1)^{n+m}\,y_n^m``,

```math
H_{n'n\,m'm} = \frac{(R_1 + \mathrm{i}R_2)\,(n'+m'-1)(n'+m')\,Y_{n'-1+n}^{m'-1-m}
  - (R_1 - \mathrm{i}R_2)\,(n'-m'-1)(n'-m')\,Y_{n'-1+n}^{m'+1-m}
  + 2R_3\,(n'-m')(n'+m')\,Y_{n'-1+n}^{m'-m}}{2\,n'\,(2n'-1)} ,
```

it contributes ``-H_{n'n\,m'm}\,\underline{u}^{(1)}_{nm}``. Hence

```math
\underline{U}^{(i)}_{n'm'}(\underline{r} + \underline{R}) =
  \sum_{j \le i}\,\sum_{n,m} (-1)^{n+m}\,\eta^{(i)(j)}_{n'n\,m'm}\;\underline{u}^{(j)}_{nm}(\underline{r}) ,
\qquad \|\underline{r}\| < \|\underline{R}\| ,
```

with ``\eta^{(i)(j)} = \alpha^{(i)(j)}`` for ``i = 1, 2`` and

```math
\eta^{(3)(3)}_{n'n\,m'm} = \alpha^{(3)(3)}_{n'n\,m'm} , \qquad
\eta^{(3)(2)}_{n'n\,m'm} = 4(1-\nu)\,\alpha^{(3)(2)}_{n'n\,m'm} , \qquad
\eta^{(3)(1)}_{n'n\,m'm} = \big(4(1-\nu) - n\big)\,\alpha^{(3)(1)}_{n'n\,m'm} - H_{n'n\,m'm} .
```

At ``n = 0`` the only regular solution is ``\underline{u}^{(3)}_{00}``. Every
coefficient is affine in ``4(1-\nu)``; [`lame_solution_translation`](@ref)
returns its two parts.

## 6. Periodic arrays

On a periodic array of sources the decaying harmonics are replaced by their
sums over the images. The sum of the harmonics of degree 2, which decrease as
``r^{-3}``, converges only conditionally, a difficulty met by
[rayleigh1892](@citet) and settled by [mcphedran1978](@citet) for the simple
cubic array. It is avoided by starting from a function periodic by
construction: the potential ``S`` of a unit source repeated on the lattice of
periods ``(L_1, L_2, L_3)``, at the points ``\underline{x}_i``, neutralized by
a uniform background [hasimoto1959](@cite),

```math
\Delta S = -4\pi\Big[\sum_{i}\delta(\underline{x} - \underline{x}_i) - \frac{1}{V}\Big] ,
\qquad \langle S\rangle = 0 , \qquad V = L_1 L_2 L_3 .
```

A decaying harmonic is a derivative of ``1/r``,
``Y_n^m = \mathcal{D}_n^m(1/r)`` with
``\mathcal{D}_n^m = (-1)^n (n-m)!\,(n+m)!\,y_n^m(\underline{\operatorname{grad}})/(2n-1)!!``
a harmonic polynomial of the partial derivatives, and its periodic counterpart
is ``Y_n^{*m} = \mathcal{D}_n^m S``, computed by
[`lattice_solid_harmonics`](@ref). ``S`` is split as [ewald1921](@citet) did,
with a parameter ``\alpha``,

```math
S(\underline{x}) = \sum_i \frac{\operatorname{erfc}(\alpha r_i)}{r_i}
  + \frac{4\pi}{V}\sum_{\underline{k}\ne\underline{0}} \frac{e^{-k^2/4\alpha^2}\,e^{\mathrm{i}\underline{k}\cdot\underline{x}}}{k^2}
  - \frac{\pi}{\alpha^2 V} ,
```

``\underline{k}`` running over the reciprocal lattice, and the derivative goes
through each sum term by term. On a radial function a harmonic polynomial of
the derivatives acts as
``\mathcal{Y}(\underline{\operatorname{grad}})\,f(r) = [(1/r)\,\mathrm{d}/\mathrm{d}r]^n f\;\mathcal{Y}(\underline{x})``,
which gives ``\mathcal{D}_n^m[\operatorname{erfc}(\alpha r)/r] = w_n(r)\,Y_n^m``,
with ``w_0 = \operatorname{erfc}(\alpha r)`` and
``w_n = w_{n-1} + (2\alpha^2 r^2)^n e^{-\alpha^2 r^2}/[(2n-1)!!\,\sqrt\pi\,\alpha r]``;
on a plane wave it multiplies by ``\mathcal{D}_n^m(\mathrm{i}\underline{k})``.
Both sums converge as a Gaussian at every degree, and their truncation follows
the precision of the number type. At a lattice point the source of that point
is left out, ``S - 1/r`` being regular there.

The background has one consequence. Near a lattice point ``S - 1/r`` is
regular but not harmonic, its Laplacian being ``4\pi/V``: its expansion holds
the term ``2\pi r^2/(3V)`` beside the harmonic ones, and the regular field the
images of a periodic dipole send to a lattice point holds the uniform term
``-4\pi\underline{d}/(3V)`` beside the translation formula with ``Y^*`` in place
of ``Y``.

For the Lamé solutions the translation also uses the products ``R_k\,Y_N^M``,
which follow from ``r^2 Y_K^\mu``; their periodic counterpart is
``B_n^{*m} = \mathcal{D}_n^m\rho``, ``\rho`` the periodic function of zero mean
with ``\Delta\rho = 2S``, the counterpart of ``r`` (``\Delta r = 2/r``,
``\mathcal{D}_n^m r = -r^2 Y_n^m/(2n-1)``), computed with ``Y^*`` by
[`lattice_solid_biharmonics`](@ref).
