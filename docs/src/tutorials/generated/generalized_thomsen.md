```@meta
EditURL = "../../../../scripts/27_generalized_thomsen.jl"
```

# [Weak anisotropy of any symmetry](@id tut-generalized-thomsen)

!!! tip "See also"
    **Theory**: [Plane waves](@ref th-plane-waves) · [Weak anisotropy of any symmetry](@ref th-weak-anisotropy) · [Thomsen's parameters](@ref th-thomsen)

    **Manual**: [Plane waves](@ref man-plane-waves)

    **Tutorials**: [Thomsen's parameters and weak anisotropy](@ref tut-thomsen) · [Stress-induced anisotropy](@ref tut-stress-anisotropy)

Thomsen's three parameters describe a transversely isotropic medium. A rock
measured in every direction may have no symmetry at all: its stiffness then
has 21 constants. To first order in the deviation from isotropy, the
velocities of such a medium depend on combinations of these constants that
generalize Thomsen's [mensch1997](@cite). This tutorial takes two triclinic
rocks, a limestone and a sandstone, and compares their exact velocities with
the first-order ones in every direction.

````@example generalized_thomsen
using TensND
using LinearAlgebra
using JSON
using Printf
using Plots
using SymPyPythonCall: symbols
gr()
nothing #hide
````

## Two triclinic rocks

[mensch1997](@citet) list in their Appendix B the stiffness matrices, in the
two-index notation (pairs 11, 22, 33, 23, 13, 12) and in GPa, of a
water-saturated Lavoux limestone and of a dry Vosges sandstone, with their
densities in kg/m³. They are stored with their source under
`data/literature`; `stiffness_from_voigt` turns a matrix into the tensor
``\mathbb{C}``.

````@example generalized_thomsen
data = JSON.parsefile(joinpath(pkgdir(TensND), "data", "literature", "menschGJI1997.json"))
voigt(key) = Float64[row[c] for row in data["tables"][key], c in ("1", "2", "3", "4", "5", "6")]
density(key) = data["quantities"][key]["value"] / 1000                  # g/cm³: GPa with km/s
rocks = [
    (name = "Lavoux limestone", ℂ = stiffness_from_voigt(voigt("lavoux_stiffness")), ρ = density("lavoux_density")),
    (name = "Vosges sandstone", ℂ = stiffness_from_voigt(voigt("vosges_stiffness")), ρ = density("vosges_density")),
]
voigt("vosges_stiffness")                       ## the matrix of the sandstone, as printed

Main.show_literature("menschGJI1997") # hide
````

## The generalized parameters

Measured from the vertical waves (``M = C_{33}`` for P, ``C_{44}`` and
``C_{55}`` for the two S waves), seventeen dimensionless numbers fix the first-order
velocities: ``\varepsilon_x``, ``\varepsilon_y``, ``\delta_x``,
``\delta_y``, ``\delta_z`` play the roles of Thomsen's ``\varepsilon`` and
``\delta`` in the three coordinate planes, ``\chi_{x,y,z}`` and the six
``\varepsilon_{IJ}`` vanish with the orthotropic symmetry, and
``\gamma_{x,y,z}`` measure the shear anisotropy
([Weak anisotropy of any symmetry](@ref th-weak-anisotropy)).

````@example generalized_thomsen
for r in rocks
    println(r.name, ": ", map(x -> round(x; digits = 4), generalized_thomsen_parameters(r.ℂ)))
end
````

## The parameters and the P wave with symbols

The parameters are ratios of the components ``C_{IJ}``, and they take
symbols. A monoclinic stiffness whose mirror is the plane normal to
``\underline{e}_3`` has thirteen components; of the parameters, ``\chi_z``,
``\varepsilon_{16}`` and ``\varepsilon_{26}`` survive next to those of
orthotropy:

````@example generalized_thomsen
C₁₁, C₁₂, C₁₃, C₁₆, C₂₂, C₂₃, C₂₆, C₃₃, C₃₆, C₄₄, C₄₅, C₅₅, C₆₆ = symbols("C_11 C_12 C_13 C_16 C_22 C_23 C_26 C_33 C_36 C_44 C_45 C_55 C_66", real = true)
o = zero(C₁₁)
monoclinic = stiffness_from_voigt(
    [
        C₁₁ C₁₂ C₁₃ o o C₁₆
        C₁₂ C₂₂ C₂₃ o o C₂₆
        C₁₃ C₂₃ C₃₃ o o C₃₆
        o o o C₄₄ C₄₅ o
        o o o C₄₅ C₅₅ o
        C₁₆ C₂₆ C₃₆ o o C₆₆
    ]
)
gₛ = map(tsimplify, generalized_thomsen_parameters(monoclinic))
````

To first order the modulus of the P wave along a unit ``\underline{n}`` is the
quartic form of [Weak anisotropy of any symmetry](@ref th-weak-anisotropy),
``\rho V_P^2 \simeq C_{33}\,[1 + q(\underline{n})]``, with

```math
\begin{aligned}
q(\underline{n}) = {} & 2\varepsilon_y\,n_1^4 + 2\varepsilon_x\,n_2^4
 + 2\delta_z\,n_1^2n_2^2 + 2\delta_y\,n_1^2n_3^2 + 2\delta_x\,n_2^2n_3^2\\
&+ 4\bigl(\varepsilon_{15}\,n_1^3n_3 + \varepsilon_{16}\,n_1^3n_2 + \varepsilon_{24}\,n_2^3n_3 + \varepsilon_{26}\,n_2^3n_1 + \varepsilon_{34}\,n_3^3n_2 + \varepsilon_{35}\,n_3^3n_1\bigr)\\
&+ 4\bigl(\chi_x\,n_1^2n_2n_3 + \chi_y\,n_1n_2^2n_3 + \chi_z\,n_1n_2n_3^2\bigr) .
\end{aligned}
```

````@example generalized_thomsen
function quartic(g, (n₁, n₂, n₃))
    return 2 * (g.εy * n₁^4 + g.εx * n₂^4 + g.δz * n₁^2 * n₂^2 + g.δy * n₁^2 * n₃^2 + g.δx * n₂^2 * n₃^2) +
        4 * (g.ε15 * n₁^3 * n₃ + g.ε16 * n₁^3 * n₂ + g.ε24 * n₂^3 * n₃ + g.ε26 * n₂^3 * n₁ + g.ε34 * n₃^3 * n₂ + g.ε35 * n₃^3 * n₁) +
        4 * (g.χx * n₁^2 * n₂ * n₃ + g.χy * n₁ * n₂^2 * n₃ + g.χz * n₁ * n₂ * n₃^2)
end
nothing #hide
````

`weak_plane_waves` normalizes its direction: for ``\underline{n}`` of any
length ``r``, ``\rho V_P^2\,r^4 = C_{33}\,[r^4 + q(\underline{n})]``, ``q`` being
homogeneous of degree four. Along a direction of three symbolic components the
difference reduces to zero:

````@example generalized_thomsen
n₁, n₂, n₃ = symbols("n_1 n_2 n_3", real = true)
ρₛ = symbols("ρ", positive = true)
ρV² = ρₛ * weak_plane_waves(monoclinic, ρₛ, [n₁, n₂, n₃]).qP[1]^2
r² = n₁^2 + n₂^2 + n₃^2
tsimplify(ρV² * r²^2 - C₃₃ * (r²^2 + quartic(gₛ, (n₁, n₂, n₃))))
````

## Exact and first-order velocities

`plane_waves` solves the Christoffel problem exactly and sorts the three
waves by decreasing velocity (qP, qS1, qS2); `weak_plane_waves` gives them to
first order, the P modulus being ``\underline{n}\cdot\boldsymbol{\Gamma}(\underline{n})\cdot\underline{n}``
and the two S moduli the eigenvalues of ``\boldsymbol{\Gamma}(\underline{n})``
restricted to the plane normal to ``\underline{n}``. A direction is given by
its angle ``\theta`` to ``\underline{e}_3`` and its azimuth ``\phi``
(`direction_from_angles`).

````@example generalized_thomsen
exact_velocities(r, n) = plane_waves(r.ℂ, r.ρ, n).velocity
weak_velocities(r, n) = collect(map(first, values(weak_plane_waves(r.ℂ, r.ρ, n))))
nothing #hide
````

Along the meridian ``\phi = 0`` (the plane of ``\underline{e}_1`` and
``\underline{e}_3``) and around the vertical axis in the horizontal plane
``\theta = 90^\circ``:

````@example generalized_thomsen
θs = range(0, π; length = 121)
ϕs = range(0, 2π; length = 121)
function velocity_panel(r, angles, direction, xlabel)
    V = [exact_velocities(r, direction(a))[k] for a in angles, k in 1:3]
    W = [weak_velocities(r, direction(a))[k] for a in angles, k in 1:3]
    p = plot(; title = r.name, titlefontsize = 10, xlabel, ylabel = "phase velocity (km/s)", xticks = 0:45:360, legend = :right)
    plot!(p, rad2deg.(angles), V; lw = 2, color = [1 2 3], label = ["qP" "qS1" "qS2"])
    plot!(p, rad2deg.(angles), W; lw = 2, ls = :dash, color = [1 2 3], label = ["first order" "" ""])
    return p
end
panels = [
    velocity_panel(r, a, d, x) for r in rocks
        for (a, d, x) in ((θs, θ -> direction_from_angles(θ, 0.0), "θ (°), ϕ = 0"), (ϕs, ϕ -> direction_from_angles(π / 2, ϕ), "ϕ (°), θ = 90°"))
]
plot(panels...; layout = (2, 2), size = (980, 680), left_margin = 4Plots.mm, bottom_margin = 4Plots.mm)
````

## The error in every direction

A wave travels at the same velocity along ``\underline{n}`` and
``-\underline{n}``, so that the upper hemisphere suffices. The relative error
``(V_{\text{first order}} - V)/V`` of each wave, on a grid of directions:

````@example generalized_thomsen
θ_grid = range(0, π / 2; length = 46)
ϕ_grid = range(0, 2π; length = 91)
relative_error(r, n) = weak_velocities(r, n) ./ exact_velocities(r, n) .- 1
errors = [[relative_error(r, direction_from_angles(θ, ϕ)) for θ in θ_grid, ϕ in ϕ_grid] for r in rocks]
@printf "%-18s %9s %9s %9s\n" "largest |error|" "qP" "qS1" "qS2"
for (r, E) in zip(rocks, errors)
    @printf "%-18s %8.2f%% %8.2f%% %8.2f%%\n" r.name (100 * maximum(e -> abs(e[k]), E) for k in 1:3)...
end
````

The map of the qP error, in percent:

````@example generalized_thomsen
maps = [
    heatmap(
        rad2deg.(ϕ_grid), rad2deg.(θ_grid), 100 .* getindex.(E, 1); title = r.name, titlefontsize = 10,
        xlabel = "ϕ (°)", ylabel = "θ (°)", xticks = 0:90:360, color = :balance, clims = (-1, 1) .* maximum(abs, 100 .* getindex.(E, 1))
    ) for (r, E) in zip(rocks, errors)
]
plot(maps...; layout = (1, 2), size = (980, 340), left_margin = 4Plots.mm, bottom_margin = 4Plots.mm)
````

## Shear-wave splitting

The two shear waves of a direction travel at different velocities: a shear
pulse splits into a fast and a slow one. `shear_splitting` gives the
relative difference ``(V_{S1} - V_{S2})/V_{S1}``, the polarization of the
fast wave and the delay per unit length of path
([Shear-wave splitting](@ref th-splitting)):

````@example generalized_thomsen
shear_splitting(rocks[2].ℂ, rocks[2].ρ, [0.0, 0.0, 1.0])
````

Its map, in percent; it vanishes in the directions where the two shear waves
have the same velocity (acoustic axes):

````@example generalized_thomsen
maps = [
    heatmap(
        rad2deg.(ϕ_grid), rad2deg.(θ_grid), [100 * shear_splitting(r.ℂ, r.ρ, direction_from_angles(θ, ϕ)).splitting for θ in θ_grid, ϕ in ϕ_grid];
        title = r.name, titlefontsize = 10, xlabel = "ϕ (°)", ylabel = "θ (°)", xticks = 0:90:360, color = :viridis
    ) for r in rocks
]
plot(maps...; layout = (1, 2), size = (980, 340), left_margin = 4Plots.mm, bottom_margin = 4Plots.mm)
````

## A transversely isotropic case

For a transversely isotropic medium of axis ``\underline{e}_3`` the
generalized parameters reduce to Thomsen's: ``\varepsilon_x = \varepsilon_y = \varepsilon``,
``\gamma_y = \gamma``, and ``\delta_x = \delta_y`` equals ``\delta`` to first
order ([Weak anisotropy of any symmetry](@ref th-weak-anisotropy)). The
Mesaverde sandstone of the same appendix, matrix (B1), checks it; its
five independent entries give the `TensTI` that `thomsen_parameters` reads:

````@example generalized_thomsen
C = voigt("mesaverde_sandstone_stiffness")
ρ = density("mesaverde_sandstone_density")
ℂ = ti_from_constants(C[1, 1], C[1, 2], C[1, 3], C[3, 3], C[4, 4], [0.0, 0.0, 1.0]; form = :components, role = :stiffness)
g = generalized_thomsen_parameters(stiffness_from_voigt(C))
t = thomsen_parameters(ℂ, ρ)
(εx = g.εx, εy = g.εy, ε = t.ε), (δx = g.δx, δy = g.δy, δ = t.δ), (γy = g.γy, γ = t.γ)
````

The parameters printed for this rock in section 4.1 of the article:

````@example generalized_thomsen
printed(name) = data["quantities"]["mesaverde_sandstone_$name"]["value"]
(ε = printed("epsilon"), δ = printed("delta"), γ = printed("gamma"))
````

``\varepsilon`` and ``\delta`` agree with those computed from the matrix to
the digits printed; ``\gamma = (C_{66} - C_{44})/(2C_{44})`` does not. The
entry ``C_{66}`` of the matrix, next to ``C_{44}(1 + 2\gamma)`` and
``C_{44}(1 + \gamma)`` with the printed ``\gamma``:

````@example generalized_thomsen
C[6, 6], C[4, 4] * (1 + 2 * printed("gamma")), C[4, 4] * (1 + printed("gamma"))
````

The computations of this page use the matrix as printed.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

