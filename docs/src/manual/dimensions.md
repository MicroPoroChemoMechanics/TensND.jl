# [Two dimensions, and any dimension](@id man-dimensions)

!!! tip "See also"
    **Theory**: [Tensors of the plane, of space and of any dimension](@ref th-dimensions)

    **Tutorials**: [The symmetry classes of a plane section](@ref tut-plane-classes) · [Elastic Green's functions](@ref tut-green-function)

The examples of the Manual are written in dimension 3, where every type and
every function exists. This page gives the calls of **dimension 2**, then those
of **any dimension** ``d``, and lists what exists in dimension 3 only. What
changes with the dimension, and why, is derived in
[Tensors of the plane, of space and of any dimension](@ref th-dimensions).

## What exists in which dimension

The table is the summary; each line is detailed below with its call.

| | Dimension 2 | Dimension 3 | Dimension 4 and above |
| :--- | :--- | :--- | :--- |
| tensors of any order, their products and contractions | yes | yes | yes |
| bases from the matrix of their vectors | yes | yes | yes |
| rotated bases | one angle, `Basis(θ)` | three angles, `Basis(θ, ϕ, ψ)` | an orthogonal matrix, `Basis(Q)` |
| random rotations and bases | yes | yes | yes |
| Kelvin–Mandel matrix | ``3\times3`` | ``6\times6`` | ``n\times n``, ``n = d(d+1)/2`` |
| isotropic tensors `TensISO` and their algebra | yes | yes | yes |
| elastic moduli of an isotropic tensor | two conventions | yes | those of the medium of dimension ``d`` |
| transversely isotropic class | no ([why](@ref th-dimensions)) | yes | no |
| cubic class | the square class, `TensCubic{2}` | yes | no |
| orthotropic class | yes, 4 or 5 constants | yes, 9 or 12 constants | no |
| projection onto a class | frame in closed form | frame searched | isotropic class only |
| coordinate systems | Cartesian, polar | all | Cartesian |

## Dimension 2

### Bases

A basis of the plane is the canonical one, `CanonicalBasis{2, T}()`, a rotation
of the canonical one by an angle, `Basis(θ)`, or the matrix of its two vectors:

```@example dims
using TensND, LinearAlgebra, SymPyPythonCall, Random
ℬ = Basis(π / 6)
```

```@example dims
angles(ℬ), isorthogonal(Basis([1.0 0.5; 0.0 1.0]))
```

### Tensors

The dimension of a tensor is that of its array of components: a ``2\times2``
matrix is a tensor of order 2 of the plane, on the canonical basis or on the
basis given.

```@example dims
ε = Tens([1.0 0.2; 0.2 -0.5])
get_order(ε), get_dim(ε), mandel(ε)
```

```@example dims
t = Tens([1.0 0.0; 0.0 2.0], ℬ)          # components in the rotated basis
canonical_components(t)
```

The products, the contractions, `tr`, `det`, `norm` and the changes of basis
([Tensors](@ref man-tensors)) are those of dimension 3.

### Isotropic tensors

`TensISO{2}` is the isotropic class of the plane: of order 2 with one
coefficient, of order 4 with two, on the projectors ``\mathbb{J}`` and
``\mathbb{K}`` of the plane. The identities and the projectors take the
dimension as `Val(2)`:

```@example dims
𝟏, 𝕀, 𝕁, 𝕂 = tens_Id2(Val(2), Val(Float64)), tens_Id4(Val(2), Val(Float64)), tens_J4(Val(2), Val(Float64)), tens_K4(Val(2), Val(Float64))
mandel(𝕁), mandel(𝕂)
```

### Elastic moduli and the two conventions

Every builder of an isotropic stiffness or compliance takes the dimension as a
first argument, `Val(2)`, and every reader takes the keyword `convention`:

```@example dims
ℂ = iso_stiffness_young_poisson(Val(2), 70.0, 0.33)
young_poisson(ℂ), young_poisson(ℂ; convention = :intrinsic)
```

!!! warning "Two pairs of moduli for one tensor of the plane"
    With `convention = :plane_strain`, the default in the plane, a
    `TensISO{4,2}` is the stiffness **in plane strain** of a material of
    space: ``E``, ``\nu``, ``\lambda``, ``\mu`` and ``M`` are those of that
    material, and the bulk modulus is the in-plane one, ``k = \lambda + \mu``,
    so that ``\mathbb{C} = 2k\,\mathbb{J} + 2\mu\,\mathbb{K}``. With
    `convention = :intrinsic` the same tensor is read as a material of the
    plane, of moduli ``E_2 = 4k\mu/(k + \mu)`` and
    ``\nu_2 = (k - \mu)/(k + \mu)``, related to the first pair by
    ``E_2 = E/(1 - \nu^2)`` and ``\nu_2 = \nu/(1 - \nu)``. The moduli
    ``k``, ``\lambda``, ``\mu`` and ``M`` do not depend on the convention. In
    plane stress the tensor is that of another material, given by
    `plane_stress_stiffness` ([Elastic moduli](@ref man-moduli-2d)).

With symbols, the two pairs read on ``\mathbb{C} = 2k\,\mathbb{J} + 2\mu\,\mathbb{K}``:

```@example dims
k, μ = symbols("k mu", positive = true)
ℂˢ = iso_stiffness(Val(2), k, μ)
simplify.(young_poisson(ℂˢ)), simplify.(young_poisson(ℂˢ; convention = :intrinsic))
```

### The orthotropic and the square classes

An orthotropic tensor of the plane has four constants
``(C_{11}, C_{22}, C_{12}, C_{66})``, with ``C_{66} = C_{1212}``, in a frame of
the plane; its Kelvin–Mandel matrix in that frame is

```math
\begin{pmatrix} C_{11} & C_{12} & 0\\ C_{12} & C_{22} & 0\\ 0 & 0 & 2C_{66}\end{pmatrix} :
```

```@example dims
𝕆 = TensOrtho(10.0, 6.0, 3.0, 2.0, ℬ)
mandel_material(𝕆)
```

The canonical matrix is that matrix turned by the Kelvin–Mandel rotation of the
frame:

```@example dims
Q = mandel_rotation(π / 6)
mandel(𝕆) ≈ Q * mandel_material(𝕆) * transpose(Q)
```

The inverse keeps the class and the frame; the product of two tensors on the
same frame is orthotropic, with five constants, ``C_{21}`` last, since it has
no major symmetry in general:

```@example dims
𝕆₂ = TensOrtho(4.0, 8.0, 1.0, 1.5, ℬ)
mandel_material(inv(𝕆)) ≈ inv(mandel_material(𝕆)), get_data(𝕆 ⊡ 𝕆₂)
```

The constants are read back on the stored data, and those of the square class,
``C_{11} = C_{22}``, by `cubic_constants`:

```@example dims
ℚ = cubic_from_constants(10.0, 4.0, 2.0, ℬ)       # C₁₁, C₁₂, C₄₄
get_data(𝕆), cubic_constants(ℚ)
```

With symbols, in a frame at the angle ``\theta``, the material matrix is
recovered from the canonical one by the rotation back:

```@example dims
C₁₁, C₂₂, C₁₂, C₆₆ = symbols("C_11 C_22 C_12 C_66", positive = true)
θ = symbols("theta", real = true)
𝕆ˢ = TensOrtho(C₁₁, C₂₂, C₁₂, C₆₆, Basis(θ))
Qˢ = mandel_rotation(θ)
simplify.(transpose(Qˢ) * mandel(𝕆ˢ) * Qˢ)
```

### From space to the plane and back

A plane of symmetry of a material of space carries two plane problems.
`plane_strain_stiffness` and `plane_stress_stiffness` give them for the
orthotropic class (in the plane normal to the axis `k` of the frame, `k = 3` by
default), the cubic class and the transversely isotropic class:

```@example dims
ℂ₃ = TensOrtho(10.0, 8.0, 9.0, 3.0, 2.0, 4.0, 2.5, 3.0, 1.5, CanonicalBasis{3, Float64}())
get_data(plane_strain_stiffness(ℂ₃)), get_data(plane_stress_stiffness(ℂ₃))
```

For any tensor of space, `plane_block(t, k)` is the restriction to the plane
normal to ``\underline{e}_k`` (the plane-strain stiffness, whatever the
anisotropy) and `plane_stress_block(ℂ, k)` the plane-stress stiffness. Their
Kelvin–Mandel matrices are a block of the matrix of space and the inverse of a
block of its inverse:

```@example dims
𝔾 = rand(Xoshiro(1), Tens{4, 3}; definite = true)
p = [1, 2, 6]                      # the slots 11, 22, √2·12 of the plane normal to e₃
mandel(plane_block(𝔾, 3)) ≈ mandel(𝔾)[p, p], mandel(plane_stress_block(𝔾, 3)) ≈ inv(inv(mandel(𝔾))[p, p])
```

`space_embedding` goes the other way: the tensor of space whose restriction is a
given tensor of the plane, uncoupled from the normal axis, with the entries
`axial` and `shear` along it (Kelvin–Mandel entries):

```@example dims
𝔼 = space_embedding(iso_stiffness(Val(2), 3.0, 1.0); axial = 5.0, shear = 2.0)
mandel(𝔼)
```

### Projections and symmetry detection

The plane has three classes, isotropic, square and orthotropic. Without a
frame, `project` finds the best one in closed form
([Projections in the plane](@ref th-projection-plane)):

```@example dims
𝕋 = Tens(get_array(cubic_from_constants(10.0, 4.0, 2.0, Basis(0.4))))   # square, its class forgotten
𝔹, d, drel = project(𝕋, :CUBIC)
cubic_constants(𝔹), drel
```

With a frame given, the projection is made in that frame, and `best_sym_tens`
returns the most constrained class within its tolerance:

```@example dims
project(𝕋, :ORTHO, Basis(0.1))[3], best_sym_tens(𝕋)[4], best_sym_tens(Tens(get_array(𝕆)))[4]
```

`:TI` is refused in the plane, which has no transversely isotropic class.
`transverse_isotropify(t)` is the average of a tensor of the plane over the
rotations of the plane ([Projection and symmetry detection](@ref man-projection)).

### Coordinates

The Cartesian system of the plane, `coorsys_cartesian(Val(2))`, and the polar
one give the differential operators in dimension 2
([Coordinate systems](@ref man-coorsystems)):

```@example dims
Plane = coorsys_cartesian(Val(2))
x, y = get_coords(Plane)
laplacian(x^3 * y, Plane)
```

```@example dims
Polar = coorsys_polar()
r, ϑ = get_coords(Polar)
simplify(laplacian(r^2 * cos(2ϑ), Polar))
```

The pointwise systems take a field as a Julia function and return a function of
the point ([Coordinate systems (numerical)](@ref man-coorsystems-num)):

```@example dims
Planeₚ = coorsys_cartesian(PointwiseGeometry, Val(2))
laplacian(X -> X[1]^3 * X[2], Planeₚ)([1.0, 2.0])
```

## Any dimension

### Tensors and bases

A tensor takes the dimension of its array, whatever it is. A basis of
dimension ``d`` is given by the matrix of its vectors; an orthogonal matrix
gives a rotated basis, without angles beyond dimension 3:

```@example dims
A = Tens(rand(Xoshiro(2), 4, 4))
R = Matrix(random_rotation(Xoshiro(3), Val(4)))      # Haar measure of SO(4)
ℬ₄ = Basis(R)
get_dim(A), typeof(ℬ₄), angles(ℬ₄)
```

```@example dims
a = rand(Xoshiro(4), 4, 4)
canonical_components(Tens(a, ℬ₄)) ≈ R * a * transpose(R)
```

### Kelvin–Mandel notation in dimension 4

A symmetric tensor of order 2 has ``n = d(d+1)/2`` Kelvin–Mandel components,
listed with the diagonal first:

```@example dims
mandel_length(4), mandel_pairs(Val(4))
```

The matrix of a tensor of order 4 with the minor symmetries is ``n\times n``,
and a rotation of the basis acts on it by the ``n\times n`` matrix
`mandel_rotation(R)`:

```@example dims
𝔾₄ = rand(Xoshiro(5), Tens{4, 4})
Q₄ = mandel_rotation(R, Val(4))
size(mandel(𝔾₄)), transpose(Q₄) * mandel(𝔾₄) * Q₄ ≈ mandel(express(𝔾₄, ℬ₄))
```

### Isotropic tensors and their moduli

`TensISO{d}` and the builders of the moduli with `Val(d)` give
``\mathbb{C} = d\,k\,\mathbb{J} + 2\mu\,\mathbb{K}``; Young's modulus and
Poisson's ratio are those of the medium of dimension ``d``:

```@example dims
ℂ₄ = iso_stiffness(Val(4), 5.0, 3.0)
young_poisson(ℂ₄), get_data(inv(ℂ₄) ⊡ ℂ₄)
```

```@example dims
simplify.(young_poisson(iso_stiffness(Val(4), k, μ)))
```

### The Cartesian system in dimension 4

`coorsys_cartesian(Val(d))` names its coordinates ``x_1, \dots, x_d``:

```@example dims
Cart₄ = coorsys_cartesian(Val(4))
x₄ = get_coords(Cart₄)
laplacian(x₄[1]^2 * x₄[4] + x₄[2] * x₄[3]^3, Cart₄)
```

## What exists in dimension 3 only

These objects need the three dimensions of space
([What is three-dimensional by nature](@ref th-dimensions)):

- the three Euler angles of a frame, `Basis(θ, ϕ, ψ)`, where the plane has one
  angle and dimension 4 an orthogonal matrix;
- the transversely isotropic class `TensTI` and its Walpole basis, since an
  axis of the plane leaves no rotation free;
- the cubic group of the cube, `TensCubic{3}`, whose plane counterpart is the
  square class;
- the curvilinear systems other than the Cartesian and the polar ones
  (cylindrical, spherical, spheroidal, ellipsoidal, toroidal) and their
  harmonics, built on the coordinates of space.

## Where to go next

- [Bases](@ref man-bases) and [Tensors](@ref man-tensors): the constructions
  shared by every dimension.
- [Structured tensors](@ref man-structured): the plane classes next to those of
  space.
- [Elastic moduli](@ref man-moduli-2d): every pair of moduli in two
  dimensions, the two conventions, plane strain and plane stress.
- [Projection and symmetry detection](@ref man-projection): the classes of the
  plane.
- [Symbolic tensors of a material symmetry](@ref man-symbolic-materials): the
  plane tensors with symbols, read back.
