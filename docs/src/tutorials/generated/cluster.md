```@meta
EditURL = "../../../../scripts/21_cluster.jl"
```

# A two-point Green tensor between offset frames

The Green operator ``\mathbb{\Gamma}(\underline{x},\underline{x}')`` linking two
points that belong to **different** spherical frames, offset along
``\underline{e}_3``. This is the kernel used when scattering by a cluster of
inclusions is treated inclusion by inclusion, each with its own local frame.

It is also a compact demonstration that tensors carrying different bases can
be combined: the position vectors live in two distinct spherical charts and
are subtracted after conversion to a common frame.

Theory: [Bases and variance](@ref th-bases-variance),
[Curvilinear differential calculus](@ref th-curvilinear).

````@example cluster
using TensND
using LinearAlgebra
using SymPy

𝕀, 𝕁, 𝕂 = iso_projectors(Val(3), Val(Sym))
𝟏 = tens_Id2(Val(3), Val(Sym))

E = symbols("E", positive = true)
ν = symbols("ν", real = true)
k = E / (3(1 - 2ν))
μ = E / (2(1 + ν))
````

## Two spherical frames

The "field" frame ``(\theta,\varphi,r)`` and the "source" frame
``(\theta',\varphi',r')``, the latter translated by ``R`` along
``\underline{e}_3``.

````@example cluster
θ, ϕ, r = symbols("θ", real = true), symbols("ϕ", real = true), symbols("r", positive = true)
S = coorsys_spherical((θ, ϕ, r))
𝐱 = getOM(S)

θ′, ϕ′, r′ = symbols("θ′", real = true), symbols("ϕ′", real = true), symbols("r′", positive = true)
S′ = coorsys_spherical((θ′, ϕ′, r′))

Cartesian = coorsys_cartesian(symbols("x y z", real = true))
𝐞₁, 𝐞₂, 𝐞₃ = unitvec(Cartesian)
R = symbols("R", positive = true)

𝐱′ = getOM(S′) + R * 𝐞₃
````

## The separation vector

Both position vectors are expressed in their own chart; subtracting them
requires a common frame, which `TensND` supplies automatically.

````@example cluster
Δ𝐱 = 𝐱 - 𝐱′
ρ = tsimplify(norm(Δ𝐱))
````

The distance between the two points, as a function of both sets of spherical
coordinates and the offset:

````@example cluster
ρ
````

The unit vector along the separation, brought back to the canonical basis:

````@example cluster
𝐍 = Δ𝐱 / ρ
𝐧 = tsimplify(change_tens_canon(𝐍))
````

## The Green operator

The same closed form as in [Elastic Green's functions](@ref), with the
separation direction ``\underline{n}`` in place of the radial vector and
``\rho`` in place of ``r``:

```math
\mathbb{\Gamma}=\frac{1}{16\pi\mu(1-\nu)\rho^{3}}
\Bigl(-3\mathbb{J}+2(1-2\nu)\mathbb{I}
+3(\boldsymbol{1}\otimes\underline{n}\otimes\underline{n}
  +\underline{n}\otimes\underline{n}\otimes\boldsymbol{1})
+12\nu\,\underline{n}\stackrel{s}{\otimes}\boldsymbol{1}\stackrel{s}{\otimes}\underline{n}
-15\,\underline{n}^{\otimes4}\Bigr)
```

````@example cluster
ℾ = 1 / (16PI * μ * (1 - ν) * ρ^3) * (
    -3𝕁 + 2(1 - 2ν) * 𝕀 + 3(𝟏 ⊗ 𝐧 ⊗ 𝐧 + 𝐧 ⊗ 𝐧 ⊗ 𝟏)
        + 12ν * 𝐧 ⊗ˢ 𝟏 ⊗ˢ 𝐧 - 15𝐧 ⊗ 𝐧 ⊗ 𝐧 ⊗ 𝐧
)
````

## Checks

``\underline{n}`` is a unit vector, whatever the two sets of coordinates:

````@example cluster
tsimplify(𝐧 ⋅ 𝐧)
````

On the axis ``\theta=\theta'=0`` the two frames align and the separation
reduces to ``|r-r'-R|``:

````@example cluster
tsimplify(subs(ρ, θ => Sym(0), θ′ => Sym(0)))
````

and the separation direction becomes ``\pm\underline{e}_3``:

````@example cluster
tsimplify(subs(get_array(𝐧)[3], θ => Sym(0), θ′ => Sym(0)))
````

## Coincident frames

Setting ``R=0`` and ``(\theta',\varphi',r')=(\theta,\varphi,r')`` reduces the
two-point kernel to the ordinary radial one — the consistency check that the
construction must pass:

````@example cluster
ρ_same = tsimplify(subs(ρ, R => Sym(0), θ′ => θ, ϕ′ => ϕ))
tsimplify(ρ_same - abs(r - r′))
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

