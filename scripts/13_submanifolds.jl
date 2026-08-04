# # Surfaces: fundamental forms and curvature
#
# A [`SubManifoldSym`](@ref) is a hypersurface embedded in ``\mathbb{R}^d``: one
# coordinate fewer than the ambient space, with the unit normal completing the
# frame. It carries both geometries — the metric it inherits, and the way it
# bends inside the ambient space.
#
# Surface indices run from ``1`` to ``d-1`` and are written ``\alpha,\beta,\gamma``;
# ambient indices run to ``d`` and are written ``i,j,k``.
#
# ```math
# a_{\alpha\beta}=\underline{a}_\alpha\cdot\underline{a}_\beta
# \qquad\text{(first fundamental form)}
# ```
#
# ```math
# b_{\alpha\beta}=\underline{n}\cdot\partial_\alpha\underline{a}_\beta
# \qquad\text{(second fundamental form)}
# ```
#
# Theory: [Submanifolds](@ref th-submanifolds).

import Pkg                                                          #jl
Pkg.activate(joinpath(@__DIR__, "..", "docs"); io = devnull)        #jl
## The `docs` environment declares every dependency the tutorials use     #jl
## (SymPy, Symbolics, NLopt, BenchmarkTools, ...) and resolves TensND      #jl
## itself through `[sources]`. These four lines are stripped from the      #jl
## generated page and notebook, which already run inside that project.     #jl

using TensND
using LinearAlgebra
using SymPy

# ## Resolving ``|\sin\theta|``
#
# SymPy cannot know that ``\theta`` lies in ``(0,\pi)``: declaring the symbol
# `real`, or even `positive`, still leaves ``|\sin\theta|`` irreducible, and
# every curvature comes back multiplied by ``\sin\theta/|\sin\theta|`` —
# correct, but unreadable. The `rules` argument of [`SubManifoldSym`](@ref) is
# exactly the mechanism for this: a rewrite applied after each simplification.

R = symbols("R", positive = true)
θ = symbols("θ", positive = true)
ϕ = symbols("ϕ", real = true)
z = symbols("z", real = true)

polar_rules = Dict(abs(sin(θ)) => sin(θ))

# ## The sphere
#
# The classical umbilical surface: every direction curves identically.

Sphere = SubManifoldSym(
    Tens(R * [sin(θ)cos(ϕ), sin(θ)sin(ϕ), cos(θ)]), (θ, ϕ), (), (R,);
    rules = polar_rules,
)

# The unit normal is the radial vector:

tsimplify.(components_canon(normal(Sphere)))

# First fundamental form — the familiar ``\mathrm{d}s^2=R^2\mathrm{d}\theta^2
# +R^2\sin^2\theta\,\mathrm{d}\varphi^2``:

𝐚 = submetric(Sphere)
get_array(𝐚)

# Second fundamental form:

𝐛 = curvature(Sphere)
tsimplify.(get_array(𝐛))

# For a sphere with the **outward** normal, ``\boldsymbol{b}=-\boldsymbol{a}/R``:

tsimplify.(get_array(𝐛) + get_array(𝐚) / R)

# ## Curvatures
#
# The shape operator is the second form with one index raised,
# ``b_\alpha{}^{\beta}=a^{\beta\gamma}b_{\gamma\alpha}``; its determinant is the Gaussian curvature and half
# its trace the mean curvature.

function curvatures(SM)
    dim = length(get_array(submetric(SM))[1, :])
    tang = 1:(dim - 1)
    a = tsimplify.(get_array(submetric(SM))[tang, tang])
    b = tsimplify.(get_array(curvature(SM))[tang, tang])
    S = tsimplify.(inv(a) * b)
    return (K = tsimplify(det(S)), H = tsimplify(tr(S) / (dim - 1)), shape = tsimplify.(diag(S)))
end

cs = curvatures(Sphere)
println("sphere : K = ", cs.K, "   H = ", cs.H, "   principal curvatures = ", cs.shape)

# ``K=1/R^2`` and ``H=-1/R``: both principal curvatures equal ``-1/R``, the sign
# following the outward orientation of ``\underline{n}``.

# ## The cylinder
#
# One principal curvature vanishes, so ``K=0``: the cylinder is **intrinsically
# flat** and can be unrolled onto a plane without distortion, although it is
# visibly curved in space. This is Gauss's *Theorema Egregium* in one line.

Cylinder = SubManifoldSym(Tens([R * cos(ϕ), R * sin(ϕ), z]), (ϕ, z), (), (R,))

tsimplify.(get_array(submetric(Cylinder))), tsimplify.(get_array(curvature(Cylinder)))

#-

cc = curvatures(Cylinder)
println("cylinder : K = ", cc.K, "   H = ", cc.H, "   principal curvatures = ", cc.shape)

# ## The plane
#
# The degenerate case: the normal is constant, so ``\boldsymbol{b}=0`` and both
# curvatures vanish.

x, y = symbols("x y", real = true)
Plane = SubManifoldSym(Tens([x, y, zero(x)]), (x, y))

tsimplify.(get_array(curvature(Plane)))

#-

cp = curvatures(Plane)
println("plane : K = ", cp.K, "   H = ", cp.H)

# ## A paraboloid
#
# A surface whose curvature actually varies with position.

Paraboloid = SubManifoldSym(Tens([x, y, (x^2 + y^2) / (2R)]), (x, y), (), (R,))
pb = curvatures(Paraboloid)

# The closed forms are large, so they are more usefully read as values. At the
# apex ``x=y=0`` the surface osculates a sphere of radius ``R``, giving
# ``|K|=1/R^2`` and ``|H|=1/R``; both fall off as one moves outwards.
#
# Note the **sign**: the ``(x,y)`` parametrization orients the normal towards
# the concave side, opposite to the outward normal used for the sphere above, so
# ``H`` comes out positive here where the sphere gave ``-1/R``. The Gaussian
# curvature ``K``, being a determinant, is insensitive to the orientation.

println("  s = √(x²+y²)      K·R²           H·R")
for s in (0.0, 0.5, 1.0, 2.0, 4.0)
    sub = Dict(x => Sym(s), y => Sym(0), R => Sym(1))
    Kv = Float64(simplify(pb.K.subs(sub)))
    Hv = Float64(simplify(pb.H.subs(sub)))
    println("     ", rpad(s, 12), rpad(round(Kv, digits = 8), 15), round(Hv, digits = 8))
end

# The apex values, symbolically:

(simplify(pb.K.subs(Dict(x => Sym(0), y => Sym(0)))), simplify(pb.H.subs(Dict(x => Sym(0), y => Sym(0)))))

# ## The Weingarten relation as a gradient
#
# Because a `SubManifoldSym` is a coordinate system, the differential operators
# apply to fields defined on the surface. The identity tying the two geometries
# together is
#
# ```math
# \mathrm{GRAD}(\underline{n})=-\,\boldsymbol{b}.
# ```

G = GRAD(normal(Sphere), Sphere)

# !!! warning "Compare tensors, not stored arrays"
#     The two sides are stored on different bases and variances, so subtracting
#     their `get_array`s gives a nonzero — and rather convincing — wrong answer.
#     On a sphere of radius ``R`` the raw arrays differ by a factor ``R^2``, so
#     the discrepancy even vanishes for ``R=1``.

tsimplify.(get_array(G)), tsimplify.(get_array(𝐛))

# Brought to a **common** basis and variance the difference is exactly zero, in
# every variance:

ℬnat = natural_basis(Sphere)
for var in ((:cov, :cov), (:cont, :cov), (:cov, :cont), (:cont, :cont))
    resid = tsimplify.(components(G, ℬnat, var) + components(𝐛, ℬnat, var))
    println("  variance ", var, " : all zero ? ", all(iszero, resid))
end

# This is the practical face of [Bases and variance](@ref th-bases-variance): a
# tensor is not its component array, and only components in the same basis may
# be compared.

# ## Intrinsic connection coefficients
#
# [`connection`](@ref) returns the **connection coefficients** ``\Gamma^\gamma_{\alpha\beta}``
# of the induced metric — Christoffel symbols, never a Riemann curvature tensor.
# (The accessor used to be called `Riemann`, which named the wrong object.)
# For the sphere:

Γ = connection(Sphere)
println("Γ^θ_ϕϕ = ", tsimplify(Γ[2, 2, 1]))
println("Γ^ϕ_θϕ = ", tsimplify(Γ[1, 2, 2]), "     Γ^ϕ_ϕθ = ", tsimplify(Γ[2, 1, 2]))

# They are not tensor components: they change inhomogeneously under a change of
# chart and can be made to vanish at any single point, which no curvature can.
# The intrinsic curvature is the ``K`` computed above from ``\boldsymbol{a}``
# and ``\boldsymbol{b}``.

# ## Summary
#
# | Surface | ``K`` | ``H`` | comment |
# |:--|:--|:--|:--|
# | plane | ``0`` | ``0`` | ``\boldsymbol{b}=0`` |
# | sphere, radius ``R`` | ``1/R^2`` | ``-1/R`` | umbilical |
# | cylinder, radius ``R`` | ``0`` | ``-1/(2R)`` | intrinsically flat |
# | paraboloid | varies | varies | ``\to1/R^2`` at the apex |
