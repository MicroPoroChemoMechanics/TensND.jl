```@meta
EditURL = "../../../../scripts/02_tensor_algebra.jl"
```

# Tensor products, contractions and their identities

The products of [Tensor algebra](@ref th-tensor-algebra), checked rather than
asserted: the index formulas, the eight algebraic identities, the two order-4
identities, and the isotropic projection they make possible.

The operators, once:

| symbol | function | index formula |
|:--|:--|:--|
| `⊗` | `otimes` | ``(\mathcal{T}\otimes\mathcal{T}')_{i\ldots j\ldots}=\mathcal{T}_{i\ldots}\mathcal{T}'_{j\ldots}`` |
| `⊗ˢ` | `sotimes` | symmetrized over the last index of the left and the first of the right |
| `⊠` | `otimesu` | ``(\boldsymbol{a}\boxtimes\boldsymbol{b})_{ijkl}=a_{ik}b_{jl}`` |
| `⊠ˢ` | `otimesul` | ``(a_{ik}b_{jl}+a_{il}b_{jk})/2`` |
| `⋅` | `dot` | one contracted index |
| `⊡` | `dcontract` | two contracted indices (pair-wise convention) |
| `⊙` | `qcontract` | four contracted indices |

````@example tensor_algebra
using TensND
using LinearAlgebra
using SymPy
using Tensors
using Random

Random.seed!(20260804)
````

## Symmetrized tensor product

``\underline{u}\stackrel{s}{\otimes}\underline{v}
 =(\underline{u}\otimes\underline{v}+\underline{v}\otimes\underline{u})/2``:

````@example tensor_algebra
a = Tens(Vec{3}((i,) -> symbols("a$i", real = true)))
b = Tens(Vec{3}((i,) -> symbols("b$i", real = true)))
a ⊗ b
````

````@example tensor_algebra
a ⊗ˢ b
````

Their difference is the antisymmetric part, which ``\otimes^s`` removes:

````@example tensor_algebra
tsimplify(get_array(a ⊗ b) - get_array(a ⊗ˢ b))
````

## The two order-4 identities

``\mathbb{1}=\boldsymbol{1}\boxtimes\boldsymbol{1}`` is the identity of *all*
order-2 tensors; ``\mathbb{I}=\boldsymbol{1}\stackrel{s}{\boxtimes}\boldsymbol{1}``
only of the **symmetric** ones.

````@example tensor_algebra
𝟏 = Matrix(1.0I, 3, 3)
𝕀 = get_array(tens_Id4(Val(3), Val(Float64)))

norm(𝕀 - 𝟏 ⊠ˢ 𝟏)
````

On a *non*-symmetric argument they differ: ``\mathbb{I}`` returns the
symmetric part, ``\mathbb{1}`` the tensor itself.

````@example tensor_algebra
m = rand(3, 3)
norm((𝟏 ⊠ 𝟏) ⊡ m - m), norm((𝟏 ⊠ˢ 𝟏) ⊡ m - m), norm((𝟏 ⊠ˢ 𝟏) ⊡ m - (m + m') / 2)
````

## The identities of the box algebra

All five hold to machine precision on random arguments.

````@example tensor_algebra
A, B, C, D = rand(3, 3), rand(3, 3), rand(3, 3), rand(3, 3)

ids = [
    "(a⊠b):(c⊠d) = (a·c)⊠(b·d)" => norm((A ⊠ B) ⊡ (C ⊠ D) - (A * C) ⊠ (B * D)),
    "(a⊠b):c     = a·c·ᵗb" => norm((A ⊠ B) ⊡ C - A * C * B'),
    "(a⊠ˢb):c    = (a·c·ᵗb + a·ᵗc·ᵗb)/2" => norm((A ⊠ˢ B) ⊡ C - (A * C * B' + A * C' * B') / 2),
    "(a⊗b):(c⊗d) = (b:c) a⊗d" => norm((A ⊗ B) ⊡ (C ⊗ D) - sum(B .* C) * (A ⊗ D)),
    "(a⊠b):(a⁻¹⊠b⁻¹) = 1⊠1" => norm((A ⊠ B) ⊡ (inv(A) ⊠ inv(B)) - 𝟏 ⊠ 𝟏),
]
for (name, residual) in ids
    println(rpad(name, 38), " residual = ", residual)
end
````

## The one that fails

There is **no** termwise inverse for ``\stackrel{s}{\boxtimes}``. Equality
requires ``\boldsymbol{a}`` and ``\boldsymbol{b}`` to be *proportional* —
commuting is not enough, and taking ``\boldsymbol{b}=\boldsymbol{1}`` does not
help either.

````@example tensor_algebra
chk(x, y) = norm((x ⊠ˢ y) ⊡ (inv(x) ⊠ˢ inv(y)) - 𝟏 ⊠ˢ 𝟏)
D1, D2 = diagm(rand(3)), diagm(rand(3))     # diagonal ⟹ they commute

for (name, r) in [
        "b = a          " => chk(A, A),
        "b = 3a         " => chk(A, 3A),
        "a, b commuting " => chk(D1, D2),
        "b = 1          " => chk(A, 𝟏),
    ]
    println(name, " residual = ", round(r, sigdigits = 4))
end
````

Only proportionality works. This is why inversion is implemented per symmetry
class rather than by one generic formula.

## Quadruple contraction is the Frobenius product

``\mathbb{J}`` and ``\mathbb{K}`` are complementary orthogonal projectors, and
their norms are the dimensions of the subspaces they project onto.

````@example tensor_algebra
𝕁 = tens_J4(Val(3), Val(Float64))
𝕂 = tens_K4(Val(3), Val(Float64))

(𝕁 ⊙ 𝕁, 𝕂 ⊙ 𝕂, 𝕁 ⊙ 𝕂)
````

``\mathbb{J}::\mathbb{J}=1`` (one spherical direction),
``\mathbb{K}::\mathbb{K}=5`` (five deviatoric ones), and they are orthogonal.
In dimension ``d`` the second is ``d(d+1)/2-1``:

````@example tensor_algebra
[(d, tens_K4(Val(d), Val(Float64)) ⊙ tens_K4(Val(d), Val(Float64)), d * (d + 1) ÷ 2 - 1) for d in 2:3]
````

## Isotropic projection

The closest isotropic tensor for the Frobenius distance, obtained by dividing
each scalar product by the corresponding norm:

```math
\mathrm{ISO}(\mathbb{T})=(\mathbb{T}::\mathbb{J})\,\mathbb{J}
+\frac{\mathbb{T}::\mathbb{K}}{5}\,\mathbb{K}
```

On a genuinely isotropic input it is exact — here a stiffness with bulk
modulus ``k`` and shear modulus ``\mu``:

````@example tensor_algebra
k, μ = symbols("k μ", positive = true)
𝕀s, 𝕁s, 𝕂s = iso_projectors(Val(3), Val(Sym))
ℂ = 3k * 𝕁s + 2μ * 𝕂s

(tsimplify(ℂ ⊙ 𝕁s / 3), tsimplify(ℂ ⊙ 𝕂s / 10))
````

recovering ``k`` and ``\mu`` exactly.

On an anisotropic input it is a genuine approximation. Take an orthotropic
stiffness and isotropize it:

````@example tensor_algebra
t = TensOrtho(10.0, 8.0, 9.0, 3.0, 2.0, 4.0, 2.5, 3.0, 1.5, CanonicalBasis{3, Float64}())
Ct = get_array(t)
kiso = (Ct ⊙ get_array(𝕁)) / 3
μiso = (Ct ⊙ get_array(𝕂)) / 10
(kiso, μiso)
````

````@example tensor_algebra
Biso, d, drel = proj_tens(:ISO, Ct)
println("closest isotropic tensor : ", Biso)
println("relative distance        : ", round(drel, sigdigits = 4))
````

The projection agrees with the closed form:

````@example tensor_algebra
norm(get_array(Biso) - (3kiso * get_array(𝕁) + 2μiso * get_array(𝕂)))
````

## Isotropization does not commute with inversion

Projecting a stiffness and projecting its compliance give **different**
isotropic materials, because the Euclidean distance is not invariant under
inversion.

````@example tensor_algebra
Siso, _, _ = proj_tens(:ISO, get_array(inv(t)))
lhs = get_array(inv(Biso))          # ISO(ℂ) then invert
rhs = get_array(Siso)               # invert then ISO

println("‖inv(ISO(ℂ)) − ISO(inv(ℂ))‖ = ", round(norm(lhs - rhs), sigdigits = 4))
println("relative                    = ", round(norm(lhs - rhs) / norm(rhs), sigdigits = 4))
````

A reported isotropic estimate must therefore always say **which** of the
stiffness or the compliance was projected. Distances that repair this are
discussed in [Isotropic tensors](@ref th-isotropic).

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

