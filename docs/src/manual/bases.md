# [Bases](@id man-bases)

How to build a basis and read what it stores. The mathematics — dual basis,
metric, variance, and why the type hierarchy exists — is on
[Bases and variance](@ref th-bases-variance).

## Constructors

[`Basis`](@ref) always returns the **most specific type** that applies, so the
cost of every later component conversion is decided here.

| Call | Returns | Meaning |
| :--- | :--- | :--- |
| `Basis{dim,T}()` | [`CanonicalBasis`](@ref) | the reference frame |
| `Basis(θ)` | [`RotatedBasis`](@ref) | 2-D rotation by one angle |
| `Basis(θ, ϕ, ψ = 0)` | [`RotatedBasis`](@ref) | 3-D rotation, Z–Y–Z Euler angles |
| `Basis(ℬ, χᵢ)` | [`OrthogonalBasis`](@ref) | `ℬ` scaled by one factor per direction |
| `Basis(v::AbstractMatrix, var)` | most specific | from a matrix of `:cov` or `:cont` vectors |
| `Basis(eᵢ, eⁱ, gᵢⱼ, gⁱʲ)` | [`Basis`](@ref) | all four matrices given explicitly |

The matrix constructors take the basis vectors as **columns**, expressed in the
canonical basis.

```@example bases
using TensND, SymPy, LinearAlgebra

ℬ = Basis(Sym[1 0 0; 0 1 0; 0 1 1])
```

```@example bases
typeof(ℬ), isorthogonal(ℬ), isorthonormal(ℬ)
```

An orthonormal basis is detected automatically and gets the cheaper type:

```@example bases
θ, ϕ, ψ = symbols("θ ϕ ψ", real = true)
typeof(Basis(θ, ϕ, ψ))
```

## Accessors

| Function | Returns |
| :--- | :--- |
| [`vecbasis`](@ref)`(ℬ, :cov)` | matrix of the ``\underline{e}_i`` |
| [`vecbasis`](@ref)`(ℬ, :cont)` | matrix of the dual ``\underline{e}^i`` |
| [`metric`](@ref)`(ℬ, :cov)` | ``[g_{ij}]`` |
| [`metric`](@ref)`(ℬ, :cont)` | ``[g^{ij}]`` |
| [`angles`](@ref)`(ℬ)` | Euler angles of a rotated basis |
| [`get_dim`](@ref)`(ℬ)` | the dimension |
| [`isorthogonal`](@ref), [`isorthonormal`](@ref) | predicates |

```@example bases
metric(ℬ, :cov)
```

```@example bases
tsimplify(metric(ℬ, :cont) * metric(ℬ, :cov))
```

The product of the two metrics is the identity, by construction — the cheapest
check that a basis is consistent.

## Normalization

`LinearAlgebra.normalize``(ℬ)` divides each vector by its norm. It
removes the scaling but **not** the obliquity: the metric acquires a unit
diagonal while the off-diagonal terms, which measure the angles between the
vectors, survive.

```@example bases
metric(normalize(ℬ), :cov)
```

## Predefined coordinates and basis vectors

The `init_*` helpers return a triple `(coordinates, vectors, basis)`:

| Function | Coordinates |
| :--- | :--- |
| [`init_cartesian`](@ref)`(dim)` | ``(x,y)`` or ``(x,y,z)`` |
| [`init_polar`](@ref) | ``(r,\theta)`` |
| [`init_cylindrical`](@ref) | ``(r,\theta,z)`` |
| [`init_spherical`](@ref) | ``(\theta,\varphi,r)`` |
| [`init_rotated`](@ref) | ``(\theta,\varphi,\psi)`` |

```@example bases2
using TensND, SymPy
(r, θ, z), (𝐞ʳ, 𝐞ᶿ, 𝐞ᶻ), ℬᶜ = init_cylindrical()
ℬᶜ
```

!!! note "The spherical ordering is ``(\theta,\varphi,r)``"
    Not ``(r,\theta,\varphi)``. The ordering is chosen so that
    ``\theta=\varphi=0`` reproduces the canonical basis **in the canonical
    order**, which makes the spherical frame a genuine
    [`RotatedBasis`](@ref). The same applies to
    [`coorsys_spherical`](@ref) — see
    [Curvilinear differential calculus](@ref th-curvilinear).

### The `canonical` keyword

Each `init_*` takes `canonical::Bool`. It decides whether the returned vectors
carry their components **in the canonical basis** (`true`) or **in the local
rotated basis** (`false`, the default):

```@example bases2
(θs, ϕs, rs), (𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ), ℬˢ = init_spherical()
components_canon(𝐞ʳ)
```

The default is usually what you want: further calculations then stay in the
rotated frame, where components are the physical ones.
