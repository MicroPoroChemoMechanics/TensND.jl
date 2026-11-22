```@meta
EditURL = "../../../../scripts/06_closest_isotropic.jl"
```

# [The isotropic tensor closest to a measured stiffness](@id tut-closest-isotropic)

!!! tip "See also"
    **Theory**: [Distances between elasticity tensors](@ref th-tensor-distances)

    **Manual**: [Distances and closest tensors](@ref man-tensor-distances)

A stiffness measured on a material has 21 constants, and a first computation
often wants two: a bulk and a shear modulus. Which isotropic tensor is the
closest to the measured one depends on what "close" means. This tutorial
takes a stiffness identified by ultrasonic measurements
[francois1998](@cite) and compares the isotropic tensors closest to it for
the distances of the package, as [Morin, Gilormini and Derrien](@cite morin2020)
do for the distances they introduce.

````@example closest_isotropic
using TensND
using LinearAlgebra
using JSON
using Printf
using Plots
gr()
nothing #hide
````

## The measured stiffness

The tensor is given by its Kelvin-Mandel matrix, in GPa. It is stored with
its source under `data/literature`.

````@example closest_isotropic
data = JSON.parsefile(joinpath(pkgdir(TensND), "data", "literature", "francoisIJSS1998.json"))
components = ("11", "22", "33", "23", "13", "12")
rows = data["tables"]["kelvin_mandel"]
ℂ = from_mandel(Float64[row[c] for row in rows, c in components])
round.(mandel(ℂ); digits = 1)

Main.show_literature("francoisIJSS1998") # hide
````

It is positive-definite, and far from isotropic: its relative Euclidean
distance to its isotropic projection is

````@example closest_isotropic
eigvals(Symmetric(Matrix(mandel(ℂ)))), last(project(ℂ, Val(:ISO)))
````

## Two Euclidean answers

The Euclidean projection of the stiffness and that of the compliance are two
isotropic materials. Their bulk moduli are close; their shear moduli are not.

````@example closest_isotropic
𝕊 = inv(ℂ)
primal = first(closest_tensor(ℂ, :ISO))
dual = inv(first(closest_tensor(𝕊, :ISO)))
bulk_shear(primal), bulk_shear(dual)
````

These are the values published for this tensor [moakher2006](@cite):

````@example closest_isotropic
published(name) = data["quantities"][name]["value"]
(published("k_primal"), published("mu_primal")), (published("k_dual"), published("mu_dual"))
````

The dual shear modulus computed here, 55.78 GPa, differs from the published
one in its fourth digit; the three other moduli agree to the digits printed.
The first tensor is the mean of Voigt, the second that of Reuss: an average
of the stiffness over the orientations, and of the compliance.

## Distances that do not choose

With a distance invariant under inversion the closest tensor is the same
whether the stiffness or the compliance is given. The log-Euclidean and the
Riemannian distances have no parameter:

````@example closest_isotropic
for (name, distance) in ("log-Euclidean" => LogEuclideanDistance(), "Riemannian" => RiemannianDistance())
    B, d = closest_tensor(ℂ, :ISO; distance)
    S, _ = closest_tensor(𝕊, :ISO; distance)
    k, μ = bulk_shear(B)
    kₛ, μₛ = bulk_shear(inv(S))
    @printf "%-14s  k = %.2f GPa  μ = %.2f GPa   from the compliance: k = %.2f  μ = %.2f   distance %.4f\n" name k μ kₛ μₛ d
end
````

Both fall between the two Euclidean answers, and they differ from each
other, the tensor not commuting with the isotropic ones.

## The effect of a gauge

The arctan- and the power-Euclidean distances compare dimensionless tensors,
divided by a gauge ``E_0``. Their closest isotropic tensor moves with it:

````@example closest_isotropic
gauges = 10 .^ range(0, 4; length = 61)
families = [
    "arctan" => E₀ -> ArctanEuclideanDistance(E₀), "power, n = 1/2" => E₀ -> PowerEuclideanDistance(1 / 2, E₀),
    "power, n = 1" => E₀ -> PowerEuclideanDistance(1, E₀), "power, n = 2" => E₀ -> PowerEuclideanDistance(2, E₀),
]
moduli = Dict(name => [bulk_shear(first(closest_tensor(ℂ, :ISO; distance = family(E₀)))) for E₀ in gauges] for (name, family) in families)
kₗ, μₗ = bulk_shear(first(closest_tensor(ℂ, :ISO; distance = LogEuclideanDistance())))
panels = map((1, 2), ("k (GPa)", "μ (GPa)"), (kₗ, μₗ), (bulk_shear(primal), bulk_shear(dual))) do i, label, logarithmic, _
    panel = plot(; xscale = :log10, xlabel = "gauge E₀ (GPa)", ylabel = label, legend = i == 2 ? :outerright : false)
    for (name, _) in families
        plot!(panel, gauges, getindex.(moduli[name], i); lw = 2, label = name)
    end
    hline!(panel, [logarithmic]; color = :black, label = "log-Euclidean")
    hline!(panel, [bulk_shear(primal)[i], bulk_shear(dual)[i]]; color = :gray, ls = :dash, label = i == 2 ? "Euclidean (stiffness, compliance)" : "")
    panel
end
plot(panels...; layout = grid(1, 2; widths = [0.4, 0.6]), size = (980, 360), left_margin = 4Plots.mm, bottom_margin = 4Plots.mm)
````

For a small gauge the dimensionless stiffness is large, ``x - 1/x`` is close
to ``x``, and the power-Euclidean distance of exponent 1 is the Euclidean
distance on the stiffnesses; for a large gauge it is that on the
compliances. The arctan-Euclidean distance has the same two limits, in the
other order:

````@example closest_isotropic
limits(family) = (bulk_shear(first(closest_tensor(ℂ, :ISO; distance = family(1.0e-4)))), bulk_shear(first(closest_tensor(ℂ, :ISO; distance = family(1.0e7)))))
limits(E₀ -> PowerEuclideanDistance(1, E₀)), limits(E₀ -> ArctanEuclideanDistance(E₀))
````

The closest isotropic material then depends on a modulus chosen freely. For
this reason [morin2020](@citet) retain, among the distances built on a
function, the log-Euclidean one, which no gauge enters.

## Other symmetry classes

The same function gives the closest tensor of the other classes, with their
axis or their frame. The distances to the classes are nested as the classes
are:

````@example closest_isotropic
frame = CanonicalBasis{3, Float64}()
for (name, distance) in ("Euclidean" => EuclideanDistance(), "log-Euclidean" => LogEuclideanDistance(), "Riemannian" => RiemannianDistance())
    d = (
        last(closest_tensor(ℂ, :ORTHO, frame; distance)), last(closest_tensor(ℂ, :CUBIC, frame; distance)),
        last(closest_tensor(ℂ, :TI, [0.0, 0.0, 1.0]; distance)), last(closest_tensor(ℂ, :ISO; distance)),
    )
    @printf "%-14s  orthotropic %.4f   cubic %.4f   transversely isotropic %.4f   isotropic %.4f\n" name d...
end
````

The distances of one line are in the unit of that line: GPa for the
Euclidean one, none for the two others.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

