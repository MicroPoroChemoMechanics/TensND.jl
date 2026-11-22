# [Distances between elasticity tensors](@id th-tensor-distances)

!!! tip "See also"
    **Manual**: [Distances and closest tensors](@ref man-tensor-distances)

    **Tutorials**: [The isotropic tensor closest to a measured stiffness](@ref tut-closest-isotropic)

How far is a stiffness from another, and which isotropic stiffness is the
closest to a measured one? The answer depends on the distance. This page
defines the distances of the package between two symmetric tensors, says which
ones treat a stiffness and its compliance alike, and derives the tensor of a
symmetry class closest to a given one in the sense of each.

## 1. A function of a symmetric tensor

A stiffness ``\mathbb{C}`` has the minor and the major symmetries. In the
[Kelvin-Mandel notation](@ref th-kelvin-mandel) it is a symmetric matrix, of
size 6 in three dimensions, whose eigenvalues ``\lambda_i`` and eigenvectors
give the spectral decomposition

```math
\mathbb{C}=\sum_i\lambda_i\,\mathbb{P}_i ,\qquad
\mathbb{P}_i:\mathbb{P}_j=\delta_{ij}\,\mathbb{P}_i,\qquad\sum_i\mathbb{P}_i=\mathbb{I} ,
```

with ``\mathbb{P}_i`` the projectors on its eigenspaces. A function ``f`` of
the real variable extends to such a tensor by its eigenvalues,

```math
f(\mathbb{C})=\sum_i f(\lambda_i)\,\mathbb{P}_i .
```

The inverse is the case ``f(x)=1/x``; the logarithm and the square root exist
when ``\mathbb{C}`` is positive-definite, and ``\exp(\ln\mathbb{C})=\mathbb{C}``.
For an isotropic tensor the projectors are ``\mathbb{J}`` and ``\mathbb{K}``:

```math
f(\alpha\,\mathbb{J}+\beta\,\mathbb{K})=f(\alpha)\,\mathbb{J}+f(\beta)\,\mathbb{K} .
```

[`tensor_function`](@ref) computes ``f(\mathbb{C})``.

## 2. The Euclidean distance, and what it lacks

The Euclidean norm of a tensor is ``\lVert\mathbb{A}\rVert^2=A_{ijkl}A_{ijkl}``,
the sum of the squares of the entries of its Kelvin-Mandel matrix, and the
Euclidean distance between two stiffnesses is

```math
d_{\mathrm{E}}(\mathbb{C}_1,\mathbb{C}_2)=\lVert\mathbb{C}_1-\mathbb{C}_2\rVert .
```

The two materials are described as well by their compliances
``\mathbb{S}_i=\mathbb{C}_i^{-1}``, and ``\lVert\mathbb{S}_1-\mathbb{S}_2\rVert`` is
another number, in another unit. For two moduli ``E_1`` and ``E_2`` the first
distance is ``|E_1-E_2|`` and the second ``|1/E_1-1/E_2|``: a pair of stiff
materials is far apart for the first and close for the second. A distance
between two *materials* should not depend on that choice. It is said to be
invariant under inversion when

```math
d(\mathbb{C}_1,\mathbb{C}_2)=d(\mathbb{C}_1^{-1},\mathbb{C}_2^{-1}) .
```

## 3. Distances built on a function

Let ``f`` be strictly monotone on the positive reals. Then

```math
d_f(\mathbb{C}_1,\mathbb{C}_2)=\lVert f(\mathbb{C}_1)-f(\mathbb{C}_2)\rVert
```

is a distance: it is the Euclidean distance between the images of the two
tensors by a one-to-one map. If moreover ``f(x)+f(1/x)`` is a constant ``c``,
then ``f(\mathbb{C}^{-1})=c\,\mathbb{I}-f(\mathbb{C})`` by section 1, the
constant cancels in the difference, and ``d_f`` is invariant under inversion
[morin2020](@cite).

**The logarithm.** ``\ln x+\ln(1/x)=0``, which gives the log-Euclidean
distance [arsigny2006](@cite)

```math
d_{\ln}(\mathbb{C}_1,\mathbb{C}_2)=\lVert\ln\mathbb{C}_1-\ln\mathbb{C}_2\rVert .
```

Since ``\ln(a\,\mathbb{C})=\ln a\;\mathbb{I}+\ln\mathbb{C}``, it does not change
when both tensors are multiplied by a number: it does not depend on the unit
of the moduli.

**A gauge.** Any other function needs dimensionless arguments. The tensors are
divided by a modulus ``E_0``, the gauge, and the distance is

```math
d_f^{E_0}(\mathbb{C}_1,\mathbb{C}_2)=\Bigl\lVert f\Bigl(\frac{\mathbb{C}_1}{E_0}\Bigr)-f\Bigl(\frac{\mathbb{C}_2}{E_0}\Bigr)\Bigr\rVert ,
```

which depends on ``E_0``. The compliances are multiplied by ``E_0``: the
invariance under inversion reads
``d_f^{E_0}(\mathbb{C}_1,\mathbb{C}_2)=d_f^{1/E_0}(\mathbb{S}_1,\mathbb{S}_2)``,
the gauge of the compliances being the compliance ``1/E_0``. Two functions
proposed in [morin2020](@cite) satisfy the condition:

- ``f=\arctan``, since ``\arctan x+\arctan(1/x)=\pi/2`` for ``x>0``. For two
  moduli the distance is the angle between the two straight lines of slopes
  ``E_1/E_0`` and ``E_2/E_0`` in the plane of the strain and of the stress
  divided by ``E_0``.
- ``f(x)=x^n-x^{-n}``, ``n>0``, for which ``f(x)+f(1/x)=0``.

## 4. The Riemannian distance

The positive-definite tensors form an open cone, on which the metric
``\mathrm{d}s^2=\lVert\mathbb{C}^{-1/2}:\mathrm{d}\mathbb{C}:\mathbb{C}^{-1/2}\rVert^2``
defines the geodesic distance [moakher2006spd](@cite)

```math
d_{\mathrm{R}}(\mathbb{C}_1,\mathbb{C}_2)=\bigl\lVert\ln\bigl(\mathbb{C}_1^{-1/2}:\mathbb{C}_2:\mathbb{C}_1^{-1/2}\bigr)\bigr\rVert
=\Bigl(\sum_i\ln^2\mu_i\Bigr)^{1/2} ,
```

where the ``\mu_i`` are the eigenvalues of ``\mathbb{C}_1^{-1}:\mathbb{C}_2``.
The eigenvalues of ``\mathbb{C}_1:\mathbb{C}_2^{-1}`` are the ``1/\mu_i``, with
the same squared logarithms: the distance is invariant under inversion. It is
also unchanged when both tensors are multiplied by a number. When the two
tensors have the same eigenprojectors, as two isotropic tensors do,
``\ln(\mathbb{C}_1^{-1}:\mathbb{C}_2)=\ln\mathbb{C}_2-\ln\mathbb{C}_1`` and the
Riemannian and the log-Euclidean distances coincide.

| distance | invariant under inversion | depends on the unit |
| :--- | :---: | :---: |
| Euclidean | no | yes |
| log-Euclidean | yes | no |
| arctan- and power-Euclidean | yes, with the gauge inverted | through the gauge |
| Riemannian | yes | no |

## 5. The closest tensor of a symmetry class

Let ``\mathcal{S}`` be a symmetry class with its axes given: the isotropic
tensors, the transversely isotropic ones about an axis, the cubic or the
orthotropic ones in a frame. It is a linear space, and the tensor of
``\mathcal{S}`` closest to ``\mathbb{C}`` for the Euclidean distance is its
orthogonal projection ``\Pi(\mathbb{C})``, which [`project`](@ref) computes
[gazis1963](@cite). On the isotropic tensors,

```math
\Pi(\mathbb{C})=(\mathbb{J}::\mathbb{C})\,\mathbb{J}+\frac{\mathbb{K}::\mathbb{C}}{5}\,\mathbb{K} .
```

The projection of the compliance is not the inverse of that of the stiffness.

**For a distance built on a function.** The tensors of ``\mathcal{S}`` are
those left unchanged by the rotations of the class, and
``f(\mathbb{Q}\star\mathbb{C})=\mathbb{Q}\star f(\mathbb{C})`` for a rotation
``\mathbb{Q}``: ``f`` maps ``\mathcal{S}`` into itself, and one-to-one onto the
part of ``\mathcal{S}`` whose eigenvalues lie in the range of ``f``. Minimizing
``\lVert f(\mathbb{C}/E_0)-f(\mathbb{B}/E_0)\rVert`` over ``\mathbb{B}`` in
``\mathcal{S}`` is then minimizing the Euclidean distance to
``f(\mathbb{C}/E_0)`` over ``\mathbb{X}=f(\mathbb{B}/E_0)`` in ``\mathcal{S}``,
whose solution is the projection. That projection is an average over the
rotations of the class, so that its eigenvalues stay between the extreme ones
of ``f(\mathbb{C}/E_0)``, in the range of ``f``. Hence

```math
\mathbb{B}=E_0\,f^{-1}\Bigl(\Pi\bigl(f(\mathbb{C}/E_0)\bigr)\Bigr) .
```

For the closest isotropic tensor ``3k\,\mathbb{J}+2\mu\,\mathbb{K}`` this is
[morin2020](@cite)

```math
f\Bigl(\frac{3k}{E_0}\Bigr)=\mathbb{J}::f\Bigl(\frac{\mathbb{C}}{E_0}\Bigr),
\qquad
f\Bigl(\frac{2\mu}{E_0}\Bigr)=\frac15\,\mathbb{K}::f\Bigl(\frac{\mathbb{C}}{E_0}\Bigr) ,
```

and, for the logarithm, ``\ln\mathbb{B}=\Pi(\ln\mathbb{C})``
[moakher2006](@cite). When ``f(x)+f(1/x)`` is a constant the tensor closest to
the compliance, with the gauge ``1/E_0``, is the inverse of ``\mathbb{B}``.

**For the Riemannian distance** no such form exists in general. The unknown is
written ``\mathbb{B}=\exp\mathbb{X}`` with ``\mathbb{X}`` in ``\mathcal{S}``,
which keeps it positive-definite, and the squared distance is the sum of the
squares of the entries of
``\ln\bigl(\mathbb{C}^{-1/2}:\exp\mathbb{X}:\mathbb{C}^{-1/2}\bigr)``: a
least-squares problem in the coefficients of ``\mathbb{X}``, solved by the
method of Gauss and Newton from the log-Euclidean solution. When
``\mathbb{C}`` commutes with the tensors of the class, as a cubic tensor does
with the isotropic ones, the two solutions coincide (section 4).

[`closest_tensor`](@ref) returns the closest tensor and its distance for each
of these distances.

## 6. Which distance?

The closest isotropic tensor depends on the distance, and for the arctan- and
power-Euclidean ones on the gauge. As the gauge tends to zero or to infinity,
the closest isotropic tensors of the arctan-Euclidean distance and of the
power-Euclidean one of exponent 1 tend to the two Euclidean ones, that of the
stiffness and the inverse of that of the compliance [morin2020](@cite); with
another exponent the limits are others. The log-Euclidean
and the Riemannian distances involve no gauge; the first has a closed form on
every class and with symbols, the second is invariant under a wider group. The
tutorial
[The isotropic tensor closest to a measured stiffness](@ref tut-closest-isotropic)
compares them on a measured tensor.
