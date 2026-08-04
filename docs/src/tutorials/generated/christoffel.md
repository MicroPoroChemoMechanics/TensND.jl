```@meta
EditURL = "../../../../scripts/11_christoffel.jl"
```

# [Christoffel symbols](@id tut-christoffel)

The connection coefficients are what distinguish a derivative in a curvilinear
chart from a plain partial derivative. This tutorial computes them from their
definition, checks that definition term by term, tabulates them for every
predefined system, confirms that the symbolic and the automatic-differentiation
routes agree, and builds a non-trivial custom chart.

```math
\frac{\partial\underline{a}_j}{\partial q^i}=\Gamma^k_{ij}\,\underline{a}_k
\qquad\Longleftrightarrow\qquad
\Gamma^k_{ij}=\frac{\partial\underline{a}_j}{\partial q^i}\cdot\underline{a}^k
```

Theory: [Curvilinear differential calculus](@ref th-curvilinear).

!!! note "Storage convention"
    [`Christoffel`](@ref) returns an array indexed
    `Γ[i,j,k]` ``=\Gamma^k_{ij}`` — the contravariant index **last**.

````@example christoffel
using TensND
using LinearAlgebra
using SymPy
using Printf
````

## Checking the definition term by term

Compute ``\partial_i\underline{a}_j\cdot\underline{a}^k`` directly from the
natural basis and compare with the stored array. In spherical coordinates:

````@example christoffel
Spherical = coorsys_spherical()
coords = getcoords(Spherical)
Γ = Christoffel(Spherical)
d = 3

residuals = [
    tsimplify(
        ∂(natvec(Spherical, j, :cov), coords[i], Spherical) ⋅ natvec(Spherical, k, :cont)
            - Γ[i, j, k],
    )
        for i in 1:d, j in 1:d, k in 1:d
]
println("all Γ recovered from the definition ? ", all(iszero, residuals))
````

Symmetry in the lower pair, ``\Gamma^k_{ij}=\Gamma^k_{ji}``, which follows
from ``\partial_i\partial_j\underline{OM}=\partial_j\partial_i\underline{OM}``:

````@example christoffel
tsimplify(maximum(abs.(Γ - permutedims(Γ, (2, 1, 3)))))
````

## Every predefined system

Only the non-vanishing symbols are listed.

````@example christoffel
function christoffel_table(name, CS)
    c = getcoords(CS)
    G = Christoffel(CS)
    dim = length(c)
    println("── ", name, "   coords = ", c, "   Lamé = ", Lame(CS))
    any_nz = false
    for i in 1:dim, j in 1:dim, k in 1:dim
        v = tsimplify(G[i, j, k])
        if !iszero(v)
            @printf "     Γ^%s_%s%s = %s\n" string(c[k]) string(c[i]) string(c[j]) string(v)
            any_nz = true
        end
    end
    any_nz || println("     (all vanish — this is a Cartesian chart)")
    return println()
end

christoffel_table("cartesian", coorsys_cartesian())
christoffel_table("polar", coorsys_polar())
christoffel_table("cylindrical", coorsys_cylindrical())
christoffel_table("spherical", coorsys_spherical())
````

The Cartesian chart is precisely the one whose connection vanishes; everything
else pays a ``\Gamma`` term.

## Why they matter

The divergence of the radial unit vector is ``2/r`` in spherical coordinates,
although the field has constant components ``(0,0,1)`` in the normalized
basis. The whole answer comes from the connection.

````@example christoffel
θ, ϕ, r = getcoords(Spherical)
𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical)
@set_coorsys Spherical

get_array(𝐞ʳ), tsimplify(DIV(𝐞ʳ))
````

## Symbolic and numerical routes agree

[`CoorSystemNum`](@ref) evaluates the same geometry pointwise with
`ForwardDiff` instead of SymPy. The Christoffel closure uses the same
``\Gamma[i,j,k]=\Gamma^k_{ij}`` convention, so the two can be compared
directly.

````@example christoffel
CSnum = coorsys_spherical_num()
x₀ = [0.7, 1.1, 2.3]                       # (θ, ϕ, r)

Γnum = Christoffel(CSnum, x₀)
Γsym = [Float64(tsimplify(Γ[i, j, k]).subs(Dict(θ => x₀[1], ϕ => x₀[2], r => x₀[3]))) for i in 1:3, j in 1:3, k in 1:3]

println("‖Γ_num − Γ_sym‖ at (θ,ϕ,r) = ", x₀, " : ", norm(Γnum - Γsym))
````

The same for the Lamé coefficients ``\chi=(r,\,r\sin\theta,\,1)``:

````@example christoffel
(Lame(CSnum, x₀), (x₀[3], x₀[3] * sin(x₀[1]), 1.0))
````

## A custom chart: prolate spheroidal coordinates

Non-trivial charts produce metric expressions SymPy will not reduce on its
own. [`CoorSystemSym`](@ref) therefore accepts auxiliary variables and rewrite
rules — `tmp_coords` standing for compound expressions, `tmp_var` substituting
them in, and `to_coords` eliminating them again before differentiation.

Prolate spheroidal coordinates ``(\varphi,p,q)`` with focal distance ``c``:

```math
\underline{OM}=c\bigl(\bar p\,\bar q\cos\varphi,\;
                      \bar p\,\bar q\sin\varphi,\;
                      p\,q\bigr),
\qquad
\bar p=\sqrt{1-p^2},\quad \bar q=\sqrt{q^2-1}.
```

````@example christoffel
Spheroidal = coorsys_spheroidal()
getcoords(Spheroidal)
````

````@example christoffel
components_canon(getOM(Spheroidal))
````

The Lamé coefficients, which without the simplification machinery come out as
unreadable nested radicals:

````@example christoffel
Lame(Spheroidal)
````

A harmonic check: the first Cartesian coordinate squared has Laplacian 2,
whatever the chart used to compute it.

````@example christoffel
simplify(LAPLACE(getOM(Spheroidal)[1]^2, Spheroidal))
````

## Spheroidal harmonics

Associated Legendre products ``P_n^m(p)\,P_n^m(q)\cos m\varphi`` are harmonic
in this chart — a genuine test of the connection, since the result is zero
only if every ``\Gamma`` term is right.

````@example christoffel
ϕs, ps, qs = getcoords(Spheroidal)
P = sympy.assoc_legendre

for (n, m) in ((2, 0), (2, 1), (3, 2))
    T = P(n, m, ps) * P(n, m, qs) * cos(m * ϕs)
    lap = simplify(LAPLACE(T, Spheroidal))
    @printf "  ∇²[P_%d^%d(p) P_%d^%d(q) cos(%dφ)] = %s\n" n m n m m string(lap)
end
````

## Building a chart from an `OM` function alone

The numerical route needs no symbolic setup at all: hand
[`CoorSystemNum`](@ref) a differentiable position map and it derives the Lamé
coefficients, the rotation to the normalized frame and the Christoffel symbols
by automatic differentiation.

Elliptic coordinates in the plane,
``\underline{OM}(\mu,\nu)=a(\cosh\mu\cos\nu,\;\sinh\mu\sin\nu)``:

````@example christoffel
a_ell = 2.0
OM_elliptic = x -> [a_ell * cosh(x[1]) * cos(x[2]), a_ell * sinh(x[1]) * sin(x[2])]
CS_ell = CoorSystemNum(OM_elliptic, 2)

x_ell = [1.0, 0.8]
χ_ell = Lame(CS_ell, x_ell)
χ_exact = a_ell * sqrt(sinh(x_ell[1])^2 + sin(x_ell[2])^2)
println("Lamé  : ", round.(χ_ell, digits = 10), "   exact (both equal) : ", round(χ_exact, digits = 10))
````

``\cosh\mu\cos\nu`` is harmonic in elliptic coordinates — it is the real part
of ``\cosh(\mu+i\nu)``:

````@example christoffel
lap_harm = LAPLACE(x -> cosh(x[1]) * cos(x[2]), CS_ell)(x_ell)
@printf "∇²(cosh μ cos ν) = %.3e   (should vanish)\n" lap_harm
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

