```@meta
EditURL = "../../../../scripts/07_plane_classes.jl"
```

# [The symmetry classes of a plane section](@id tut-plane-classes)

!!! tip "See also"
    **Theory**: [Tensors of the plane, of space and of any dimension](@ref th-dimensions)

    **Manual**: [Two dimensions, and any dimension](@ref man-dimensions)

    **Tutorials**: [Elastic Green's functions](@ref tut-green-function)

A stiffness of the plane has six constants, and four symmetry classes:
isotropic, square, orthotropic, or none. This tutorial takes a stiffness
identified by ultrasonic measurements [francois1998](@cite), cuts it by its
three coordinate planes, and asks of each section how far it is from each
class and in which frame. In the plane the answer needs no search: the frame
of the closest tensor of a class is known in closed form
([Projections in the plane](@ref th-projection-plane)).

````@example plane_classes
using TensND
using LinearAlgebra
using JSON
using Plots
gr()
nothing #hide
````

## The measured stiffness and its three sections

The tensor is given by its Kelvin-Mandel matrix, in GPa, stored with its
source under `data/literature`.

````@example plane_classes
data = JSON.parsefile(joinpath(pkgdir(TensND), "data", "literature", "francoisIJSS1998.json"))
components = ("11", "22", "33", "23", "13", "12")
ℂ = from_mandel(Float64[row[c] for row in data["tables"]["kelvin_mandel"], c in components])
nothing #hide

Main.show_literature("francoisIJSS1998") # hide
````

In plane strain normal to the axis ``\underline{e}_k`` the strains that
carry the index ``k`` vanish, and the stiffness of the plane is the block of
``\mathbb{C}`` on the in-plane components: `plane_block(ℂ, k)`, a tensor of
the plane, whose axes are the two other axes in the cyclic order.

````@example plane_classes
sections = [plane_block(ℂ, k) for k in 1:3]
mandel(sections[3])
````

## How far from each class

`project` returns the closest tensor of a class and its relative distance.
For the square and the orthotropic classes the frame is free:

````@example plane_classes
relative_distance(t, class) = round(project(t, class)[3]; digits = 3)
[(k = k, iso = relative_distance(t, :ISO), square = relative_distance(t, :CUBIC), ortho = relative_distance(t, :ORTHO)) for (k, t) in enumerate(sections)]
````

The classes are nested, so that the distances decrease from left to right.
None of the sections is close to isotropy, and each is within two percent of
the square class: the measured crystal is close to cubic, which its sections
show in a frame that is not that of the measurement.

## The frame, in closed form

Under a rotation of the frame by ``\theta`` a stiffness of the plane splits
into two invariants and two pairs that turn by ``2\theta`` and ``4\theta``.
The distance to the square class involves the second pair only, and is least
at ``\theta=\operatorname{atan2}(a,b)/4``; the distance to the orthotropic
class is least at a root of a quartic. Both are exact, which the distance in
every frame confirms for the section normal to ``\underline{e}_3``:

````@example plane_classes
t = sections[3]
angle_of(B) = (E = vecbasis(TensND.frame(B), :cov); mod(atan(E[2, 1], E[1, 1]), π / 2))
θ_ortho, θ_square = angle_of(project(t, :ORTHO)[1]), mod(angle_of(project(t, :CUBIC)[1]), π / 4)
θs = range(0, π / 2; length = 361)
d_ortho = [project(t, :ORTHO, Basis(θ))[3] for θ in θs]
d_square = [project(t, :CUBIC, Basis(θ))[3] for θ in θs]
plot(rad2deg.(θs), [d_ortho d_square]; label = ["orthotropic" "square"], lw = 2, xlabel = "angle of the frame (°)",
    ylabel = "relative distance to the class", xticks = 0:15:90, size = (640, 360), left_margin = 4Plots.mm, bottom_margin = 4Plots.mm)
scatter!([rad2deg(θ_ortho)], [project(t, :ORTHO)[3]]; label = "closed form", color = 1)
scatter!(rad2deg.([θ_square, θ_square + π / 4]), fill(project(t, :CUBIC)[3], 2); label = "", color = 2)
````

The orthotropic class has a period of a quarter turn, and the square class of
an eighth: the diagonals of a square are mirror lines too. The two curves
nearly coincide: the pair that turns by ``2\theta`` is small in this section.

## What the projections keep

The directional modulus ``E(\varphi)=1/(\underline{n}\otimes\underline{n}:\mathbb{S}:\underline{n}\otimes\underline{n})``,
with ``\underline{n}=(\cos\varphi,\sin\varphi)`` and ``\mathbb{S}`` the
compliance, shows what each projection keeps of the section. It is read on
the Kelvin-Mandel matrix of the compliance:

````@example plane_classes
function directional_modulus(t, φ)
    n̂ = [cos(φ)^2, sin(φ)^2, √2 * sin(φ) * cos(φ)]          # Kelvin-Mandel components of n⊗n
    return 1 / (n̂' * operator_matrix(inv(t)) * n̂)
end
φs = range(0, 2π; length = 361)
curves = (("section", t), ("orthotropic", project(t, :ORTHO)[1]), ("square", project(t, :CUBIC)[1]), ("isotropic", project(t, :ISO)[1]))
plt = plot(; proj = :polar, size = (520, 520), legend = :outerbottom, legend_columns = 4)
for (name, s) in curves
    plot!(plt, φs, [directional_modulus(s, φ) for φ in φs]; label = name, lw = name == "section" ? 3 : 2)
end
plt
````

The orthotropic and the square projections follow the section closely, the
second with the symmetry of a quarter turn; the isotropic one is a circle,
far from the section.

## Detecting the class

`best_sym_tens` tries the classes from the most to the least symmetric and
returns the first one within a tolerance, here of five percent:

````@example plane_classes
[last(best_sym_tens(s; ε = 0.05)) for s in sections]
````

A tensor built in a class is recognized whatever the frame it was built in,
and its frame is recovered:

````@example plane_classes
𝕆 = TensOrtho(10.0, 6.0, 3.0, 2.0, Basis(0.7))
𝔹, _, _, class = best_sym_tens(Tens(get_array(𝕆)))
class, angle_of(𝔹), get_data(𝔹)
````

## Plane strain and plane stress

In plane stress the stresses that carry the index ``k`` vanish instead of the
strains, and the stiffness of the plane is the Schur complement of the
out-of-plane block, `plane_stress_block`. The two sections of the same plane
are different tensors of the plane; here the in-plane bulk and shear moduli
of their closest isotropic tensors, in GPa:

````@example plane_classes
strain, stress = plane_block(ℂ, 3), plane_stress_block(ℂ, 3)
moduli(s) = round.(bulk_shear(project(s, :ISO)[1]); digits = 1)
moduli(strain), moduli(stress)
````

## Where to go next

- [Tensors of the plane, of space and of any dimension](@ref th-dimensions):
  the classes of the plane, the harmonics and the closed forms used here.
- [Two dimensions, and any dimension](@ref man-dimensions): every call that
  takes a tensor of the plane.
- [The isotropic tensor closest to a measured stiffness](@ref tut-closest-isotropic):
  the same measured tensor, in space, for several distances.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

