# [Two dimensions, and any dimension](@id man-dimensions)

The examples of the Manual are mostly in dimension 3. This page says what is
built in **dimension 2**, with the call that builds it, then what holds in
**any dimension**, and what does not exist outside dimension 3. Each of the
other pages of the Manual has a section *In two dimensions* that refers here.

| | Dimension 2 | Dimension 3 | Dimension 4 and above |
| :--- | :--- | :--- | :--- |
| tensors of any order on the canonical basis | yes | yes | yes |
| general bases, from the matrix of their vectors | yes | yes | yes |
| rotated bases, from angles | yes, one angle | yes, three angles | no |
| isotropic tensors `TensISO`, their products and inverses | yes | yes | yes |
| builders and readers of the elastic moduli | yes, in plane strain | yes | no |
| Kelvin–Mandel matrix | yes, ``3\times3`` | yes, ``6\times6`` | no |
| transversely isotropic, cubic, orthotropic classes | no | yes | no |
| closest isotropic tensor | yes | yes | yes |
| coordinate systems and differential operators | Cartesian, polar | all | no |

## Dimension 2

### Bases

The canonical basis of the plane is `CanonicalBasis{2, T}()`, and a basis
rotated by one angle is `Basis(θ)`:

```@example dims
using TensND, LinearAlgebra
ℬ = Basis(π / 6)
```

### Tensors

The dimension of a tensor is that of its array of components: a ``2\times2``
matrix is a tensor of order 2 of the plane, on the canonical basis or on the
basis given.

```@example dims
ε = Tens([1.0 0.2; 0.2 -0.5])
get_order(ε), get_dim(ε)
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

The Kelvin–Mandel matrix of a tensor of order 4 of the plane is ``3\times3``,
on the components ``(11, 22, \sqrt2\,12)``.

### Elastic moduli

Every builder of an isotropic stiffness or compliance takes the dimension as a
first argument, `Val(2)`, and every reader accepts the result:

```@example dims
ℂ = iso_stiffness_young_poisson(Val(2), 70.0, 0.33)
young_poisson(ℂ), bulk_shear(ℂ), lame_shear(ℂ)
```

```@example dims
𝕊 = iso_compliance_young_poisson(Val(2), 70.0, 0.33)
inv(ℂ) ≈ 𝕊, young_poisson(𝕊; role = :compliance)
```

!!! warning "What the moduli of a plane tensor are"
    A `TensISO{4,2}` built by these functions is the stiffness **in plane
    strain** of a 3-D material: ``E``, ``\nu``, ``\lambda``, ``\mu`` and ``M``
    are those of that 3-D material. Its bulk modulus is the **in-plane** one,
    ``k = \lambda + \mu``, so that
    ``\mathbb{C} = 2k\,\mathbb{J} + 2\mu\,\mathbb{K}``, and not the 3-D bulk
    modulus ``\lambda + \tfrac23\mu``. In plane stress the tensor is that of
    another 3-D material, given by `plane_stress_stiffness` and
    `plane_stress_to_plane_strain`
    ([Elastic moduli](@ref man-moduli-2d)).

With symbols, the matrix shows where the moduli go:

```@example dims
using SymPyPythonCall
k, μ = symbols("k mu", positive = true)
ℂˢ = iso_stiffness(Val(2), k, μ)
mandel(ℂˢ)
```

```@example dims
bulk_shear(ℂˢ), simplify.(young_poisson(ℂˢ))
```

The plane-strain and the plane-stress stiffnesses of a 3-D material, isotropic
or transversely isotropic about the normal to the plane, are plane tensors:

```@example dims
E, ν = symbols("E nu", positive = true)
ℂ₃ = iso_stiffness_young_poisson(E, ν)
simplify.(mandel(plane_strain_stiffness(ℂ₃)))
```

```@example dims
simplify.(mandel(plane_stress_stiffness(ℂ₃)))
```

### Anisotropic tensors of the plane

The classes `TensTI`, `TensCubic` and `TensOrtho` are those of dimension 3:
the plane has no structured class but the isotropic one. An anisotropic
stiffness of the plane is a general tensor, entered by its components or by
its ``3\times3`` Kelvin–Mandel matrix; its products and its inverse are
computed on the full components.

```@example dims
M = [10.0 3.0 0.0; 3.0 6.0 0.0; 0.0 0.0 4.0]        # orthotropic in the plane, in its axes
ℂₒ = from_mandel(M)
ℂₒ ⊡ ε
```

```@example dims
mandel(inv(ℂₒ)) ≈ inv(M)
```

The closest isotropic tensor of the plane is found as in dimension 3:

```@example dims
best_fit_iso(ℂₒ)
```

### Coordinates

The Cartesian system of the plane and the polar one give the differential
operators in dimension 2 ([Coordinate systems](@ref man-coorsystems)), and
`PolarChart()` the conversions of coordinates ([Charts](@ref man-charts)):

```@example dims
Polar = coorsys_polar()
r, θ = get_coords(Polar)
simplify(laplacian(r^2 * cos(2θ), Polar))
```

## Any dimension

A tensor takes the dimension of its array, whatever it is, on the canonical
basis; the isotropic tensors, their products and their inverses are written
for any dimension as well.

```@example dims
A = Tens(rand(4, 4))                       # order 2, dimension 4
get_dim(A), tr(A) ≈ sum(A[i, i] for i in 1:4)
```

```@example dims
𝕃 = TensISO{4}(3.0, 2.0)                   # 3 𝕁 + 2 𝕂 in dimension 4
get_dim(𝕃), inv(𝕃)
```

A basis is given by the matrix of its vectors, and the closest isotropic
tensor is found, in any dimension:

```@example dims
ℬ₄ = Basis([1.0 0.3 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 2.0 0.0; 0.0 0.0 0.0 1.0])
get_dim(Tens(rand(4, 4), ℬ₄)), best_fit_iso(Tens(rand(4, 4, 4, 4))) isa TensISO{4, 4}
```

Outside dimensions 2 and 3, the following are not available: the bases
rotated by angles, the Kelvin–Mandel matrices, the builders and the readers of
the elastic moduli, the anisotropic classes and the coordinate systems.

## Where to go next

- [Bases](@ref man-bases) and [Tensors](@ref man-tensors): the constructions
  shared by every dimension.
- [Elastic moduli](@ref man-moduli-2d): every pair of moduli in two
  dimensions, plane strain and plane stress.
- [Symbolic tensors of a material symmetry](@ref man-symbolic-materials): the
  plane isotropic tensors with symbols, read back.
- [Coordinate systems (symbolic)](@ref man-coorsystems): the polar system.
