# [Random tensors](@id man-random-tensors)

[`rand`](@ref Base.rand(::Type{<:TensND._RandomClass})) draws a tensor of any class: isotropic, transversely
isotropic, orthotropic, cubic or generic, with random coefficients and, for a
class with a direction, a random orientation. The orientations are uniform:
no direction of space is favored. A random tensor is what a test needs to
reach the cases a hand-picked one misses: an axis that lies along no
coordinate axis, a frame that shares nothing with the canonical one, a
coupling of an unexpected sign.

## Drawing a tensor

The class is the type of the tensor, with its default shape (the one with the
major symmetry) and `Float64` coefficients, or written in full for another shape
or element type. Without a generator, Julia's default one is used:

```@example rnd
using TensND, LinearAlgebra, Random

rand(TensCubic; definite = true)
```

A generator of one's own, passed first, makes the draws reproducible, as a
test wants them:

```@example rnd
rng = Xoshiro(2026)
ℂ = rand(rng, TensCubic)              # a cubic tensor on a random cube
𝕋 = rand(rng, TensTI{4})              # about a random axis, 5 coefficients
𝔸 = rand(rng, TensTI{4, Float64, 8})  # the shape without major symmetry
𝕆 = rand(rng, TensOrtho{4})           # 9 constants, a random frame
𝔾 = rand(rng, Tens{4, 3})             # generic, minor and major symmetries
typeof(𝕋), axis(𝕋)
```

As for numbers, a size after the class draws an array, every entry from the
default law of the class; with keywords, the array is a comprehension:

```@example rnd
crystals = rand(rng, TensCubic, 3)                                   # three random cubes
stiff = [rand(rng, TensCubic; definite = true) for _ in 1:3]         # three positive definite ones
typeof(crystals), all(t -> minimum(eigvals(Symmetric(Matrix(mandel(t))))) > 0, stiff)
```

| class | shapes | direction |
|:--|:--|:--|
| `TensISO{order, dim}` | order 2 or 4 | none |
| `TensTI{4}`, `TensTI{2}` | `N = 5, 6, 8` and `N = 2, 3` | `axis` or `frame` |
| `TensOrtho{4}`, `TensOrtho{2}` | `N = 9, 12` and `N = 3` | `frame` |
| `TensCubic` | | `frame` |
| `Tens{order, dim}` | `symmetry = Val(:major)`, `Val(:minor)`, `Val(:none)` at order 4; `Val(:symmetric)`, `Val(:none)` at order 2 | none |

The rotations, axes and frames are drawn by [`random_rotation`](@ref),
[`random_axis`](@ref) and [`random_basis`](@ref), which serve on their own as
well.

## A fixed direction

Without `axis` or `frame`, the orientation is random. Passing one fixes it, and
only the coefficients are drawn; a transversely isotropic tensor takes either,
its axis being the third vector of the frame.

```@example rnd
𝕋₃ = rand(rng, TensTI{4}; axis = (0.0, 0.0, 1.0))         # about e₃
ℂ₀ = rand(rng, TensCubic; frame = CanonicalBasis{3, Float64}())   # on the canonical cube
axis(𝕋₃), frame(ℂ₀)
```

## Positive definite tensors

By default the coefficients take either sign. With `definite = true` they make
the tensor positive definite, as a stiffness or a conductivity is, whatever its
orientation; only the shapes with the major symmetry accept it.

```@example rnd
𝕊 = rand(rng, TensTI{4}; definite = true)
minimum(eigvals(Symmetric(Matrix(mandel(𝕊))))) > 0
```

A generic positive definite tensor of order 4 is the random stiffness of a
triclinic material:

```@example rnd
ℂ₁ = rand(rng, Tens{4, 3}; definite = true)
eigvals(Symmetric(Matrix(mandel(ℂ₁))))
```

## The laws

**Orientations.** A rotation follows the Haar measure of ``SO(3)``, the law
that no rotation of space changes, and an axis is uniform on the unit sphere.
Drawing the three Euler angles uniformly does not give these laws: the axis
would fall near the poles as often as near the equator. Why, and how the
uniform law is drawn, is on [Uniform random rotations](@ref th-random-rotations).

**Coefficients.** The Kelvin–Mandel matrix of a tensor in its own frame has
the blocks of its class: the ``2\times 2`` block coupling the axial and the
transverse normal strains of a transversely isotropic tensor, the
``3\times 3`` normal block of an orthotropic one, the eigenvalues of a cubic one,
the whole ``6\times 6`` matrix of a generic one. A symmetric block is drawn from
the Gaussian orthogonal ensemble, normal entries of variance 1 on the diagonal
and ``1/2`` off it; positive definite, it has Haar eigenvectors and log-normal
eigenvalues; a block without symmetry has independent standard normal
entries. Each law is unchanged by an orthogonal change of basis of its block.
Since a rotation acts on the Kelvin–Mandel matrix as an orthogonal change of
basis, the law of a generic tensor is invariant under rotations, and needs no
random frame:

```@example rnd
ℛ = random_basis(rng)
km(t) = Matrix(mandel(t))
m(f) = sum(_ -> km(f()) .^ 2, 1:4000) / 4000      # mean squares of the entries
a = m(() -> rand(rng, Tens{4, 3}))                     # read in the canonical basis
b = m(() -> express(rand(rng, Tens{4, 3}), ℛ))     # read in a rotated one
round.(a; digits = 1), maximum(abs, a - b) < 0.2
```

## In tests

A random tensor makes a test independent of the axes it was written in, and
the same seed makes it reproducible. The algebra of the structured types is
checked this way, each pair of classes on independent random frames against
the generic route:

```@example rnd
A = rand(rng, TensCubic)
B = rand(rng, TensTI{4})
Array(canonical_components(A + B)) ≈ Array(canonical_components(A)) + Array(canonical_components(B))
```
