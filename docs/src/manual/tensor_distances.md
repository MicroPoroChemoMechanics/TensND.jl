# [Distances and closest tensors](@id man-tensor-distances)

!!! tip "See also"
    **Theory**: [Distances between elasticity tensors](@ref th-tensor-distances)

    **Tutorials**: [The isotropic tensor closest to a measured stiffness](@ref tut-closest-isotropic)

Three functions measure how far two symmetric tensors are from each other and
find the tensor of a symmetry class closest to a given one, for a distance
that is chosen.

| to … | write |
| :--- | :--- |
| take a function of a tensor | `tensor_function(log, ℂ)`, `tensor_function(sqrt, ℂ)`, … |
| measure a distance | `tensor_distance(ℂ₁, ℂ₂, d)` |
| find the closest tensor of a class | `closest_tensor(ℂ, :ISO; distance = d)`, `closest_tensor(ℂ, :TI, axis; distance = d)`, … |
| choose the distance `d` | `EuclideanDistance()` (the default), `LogEuclideanDistance()`, `RiemannianDistance()`, `ArctanEuclideanDistance(E₀)`, `PowerEuclideanDistance(n, E₀)` |

```@setup distances
using TensND, LinearAlgebra
```

## A function of a tensor

[`tensor_function`](@ref) applies a function to the eigenvalues of a symmetric
tensor and keeps its eigenprojectors. On a structured tensor it acts on the
coefficients, and the class is kept:

```@example distances
ℂ = iso_stiffness_young_poisson(70.0, 0.33)
tensor_function(log, ℂ)
```

On a general tensor it goes through the eigenvalues of the Kelvin-Mandel
matrix. The square root of a stiffness, contracted with itself, gives the
stiffness back:

```@example distances
ℂₐ = Tens(get_array(ℂ)) + Tens(get_array(ti_from_constants(30.0, 10.0, 10.0, 90.0, 20.0, [0.0, 0.0, 1.0])))   # a stiffness with an axis
ℝ = tensor_function(sqrt, ℂₐ)
norm(get_array(ℝ ⊡ ℝ) - get_array(ℂₐ))
```

## Distances

[`tensor_distance`](@ref) takes the two tensors and the distance. The
Euclidean distance, the default, is not the same on two stiffnesses and on
their compliances; the three others below are:

```@example distances
ℂ₁, ℂ₂ = iso_stiffness(10.0, 4.0), iso_stiffness(14.0, 3.0)
[(tensor_distance(ℂ₁, ℂ₂, d), tensor_distance(inv(ℂ₁), inv(ℂ₂), d)) for d in (EuclideanDistance(), LogEuclideanDistance(), RiemannianDistance())]
```

A distance built on another function than the logarithm needs a gauge, a
modulus that makes the tensors dimensionless. For the compliances the gauge is
its inverse:

```@example distances
E₀ = 20.0
tensor_distance(ℂ₁, ℂ₂, ArctanEuclideanDistance(E₀)), tensor_distance(inv(ℂ₁), inv(ℂ₂), ArctanEuclideanDistance(1 / E₀))
```

Any strictly monotone function may be given with its inverse,
`FunctionEuclideanDistance(f, f⁻¹, E₀)`; the distance is invariant under
inversion when `f(x) + f(1/x)` is a constant.

With symbolic moduli the distance between two isotropic tensors is a closed
form:

```@example distances
using SymPyPythonCall
k₁, μ₁, k₂, μ₂ = symbols("k_1 μ_1 k_2 μ_2", positive = true)
tensor_distance(iso_stiffness(k₁, μ₁), iso_stiffness(k₂, μ₂), LogEuclideanDistance())
```

## The closest tensor of a symmetry class

[`closest_tensor`](@ref) returns the tensor of a class closest to a tensor,
and the distance between them. The class is given as for
[`project`](@ref): `:ISO`, `:TI` with its axis, `:CUBIC` or `:ORTHO` with their
frame. With the default distance the result is that of `project`:

```@example distances
ℂc = TensCubic(30.0, 8.0, 12.0, CanonicalBasis{3, Float64}())    # coefficients on 𝕁 and on the two cubic projectors
[bulk_shear(first(closest_tensor(ℂc, :ISO; distance = d))) for d in (EuclideanDistance(), LogEuclideanDistance(), RiemannianDistance())]
```

The Euclidean shear modulus is the arithmetic mean of the two shear moduli of
the cube with the weights 2/5 and 3/5, the two others their geometric mean
with the same weights. With a distance invariant under inversion the closest
tensor does not depend on whether the stiffness or the compliance is given:

```@example distances
B, _ = closest_tensor(ℂc, :ISO; distance = LogEuclideanDistance())
S, _ = closest_tensor(inv(ℂc), :ISO; distance = LogEuclideanDistance())
bulk_shear(B), bulk_shear(inv(S))
```

The same with symbols:

```@example distances
α, β, γ = symbols("α β γ", positive = true)
simplify.(get_data(first(closest_tensor(TensCubic(α, β, γ, CanonicalBasis{3, Float64}()), :ISO; distance = LogEuclideanDistance()))))
```

## What each distance accepts

| distance | numbers | duals | symbols |
| :--- | :---: | :---: | :---: |
| Euclidean | yes | yes | yes |
| log-, arctan-, power-Euclidean | yes | yes | structured tensors |
| Riemannian, the distance | yes | yes | isotropic tensors |
| Riemannian, the closest tensor | yes | no | no |

The functions of a general tensor are taken on its eigenvalues, which symbols
do not give; a structured tensor needs none. The closest tensor in the
Riemannian sense is the result of an iteration.
