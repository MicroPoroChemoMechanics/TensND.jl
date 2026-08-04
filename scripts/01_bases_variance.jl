# # Bases, variance and the metric
#
# The one notion `TensND` has that a Cartesian-only tensor library does not: a
# tensor carries a **basis**, which need not be orthonormal, and a **variance**
# saying how each index transforms. This tutorial builds an oblique basis, moves
# components between variances, and checks the identities of
# [Bases and variance](@ref th-bases-variance) — symbolically and numerically.
#
# The definitions in one place: for a basis ``(\underline{e}_i)`` the dual basis
# ``(\underline{e}^i)`` satisfies ``\underline{e}^i\cdot\underline{e}_j=\delta^i_j``,
# and the two metrics are
#
# ```math
# g_{ij}=\underline{e}_i\cdot\underline{e}_j,
# \qquad
# g^{ij}=\underline{e}^i\cdot\underline{e}^j,
# \qquad
# g^{ik}g_{kj}=\delta^i_j .
# ```

import Pkg                                                          #jl
Pkg.activate(joinpath(@__DIR__, "..", "docs"); io = devnull)        #jl
## The `docs` environment declares every dependency the tutorials use     #jl
## (SymPy, Symbolics, NLopt, BenchmarkTools, ...) and resolves TensND      #jl
## itself through `[sources]`. These four lines are stripped from the      #jl
## generated page and notebook, which already run inside that project.     #jl

using TensND
using LinearAlgebra
using SymPy
using Tensors

# ## An oblique basis
#
# The constructor takes the new basis vectors as **columns** expressed in the
# canonical basis. Nothing here is orthogonal or normalized.

ℬ = Basis(Sym[1 0 0; 0 1 0; 0 1 1])

# The four stored matrices. `vecbasis(ℬ, :cov)` holds the ``\underline{e}_i``,
# `vecbasis(ℬ, :cont)` the dual vectors ``\underline{e}^i``.

vecbasis(ℬ, :cov)

#-

vecbasis(ℬ, :cont)

# The covariant metric ``g_{ij}`` and its inverse ``g^{ij}``:

metric(ℬ, :cov)

#-

metric(ℬ, :cont)

# ## The defining identities
#
# ``g^{ik}g_{kj}=\delta^i_j`` and
# ``\underline{e}^i\cdot\underline{e}_j=\delta^i_j``:

tsimplify(metric(ℬ, :cont) * metric(ℬ, :cov))

#-

tsimplify(transpose(vecbasis(ℬ, :cont)) * vecbasis(ℬ, :cov))

# The basis is neither orthogonal nor orthonormal, which is exactly why variance
# will matter below:

isorthogonal(ℬ), isorthonormal(ℬ)

# ## Covariant and contravariant components of a vector
#
# The same vector, two sets of components. `TensND` stores a tensor together
# with its basis and variance, and [`components`](@ref) converts.

V = Tens(Tensor{1, 3}(i -> symbols("v$i", real = true)))

# Contravariant components ``v^i`` (the coefficients on ``\underline{e}_i``):

components(V, ℬ, (:cont,))

# Covariant components ``v_i`` (the projections on ``\underline{e}_i``):

components(V, ℬ, (:cov,))

# They are related by the metric, ``v_i = g_{ij}v^j``:

tsimplify(metric(ℬ, :cov) * components(V, ℬ, (:cont,)) - components(V, ℬ, (:cov,)))

# ## Order 2: four variance combinations
#
# An order-``p`` tensor has one choice per index, so an order-2 tensor has four.

T = Tens(Tensor{2, 3}((i, j) -> symbols("t$i$j", real = true)))

#-

components(T, ℬ, (:cov, :cov))

#-

tfactor(tsimplify(components(T, ℬ, (:cont, :cov))))

# ## The cheapest sanity check on a basis
#
# Store the covariant metric as a twice-covariant tensor and raise one index.
# Since ``g^{ik}g_{kj}=\delta^i_j``, the answer *must* be the identity — whatever
# the basis.

G = Tens(metric(ℬ, :cov), ℬ, (:cov, :cov))
components(G, (:cont, :cov))

# ## Normalization removes the scaling, not the obliquity
#
# Dividing each vector by its norm gives unit vectors, so the metric acquires a
# unit diagonal — but the off-diagonal terms, which measure the angles between
# the vectors, survive. Variance still matters.

ℬ̄ = normalize(ℬ)
metric(ℬ̄, :cov)

#-

isorthogonal(ℬ̄), isorthonormal(ℬ̄)

# ## Where variance becomes invisible
#
# On an orthonormal basis ``g_{ij}=g^{ij}=\delta_{ij}``, so the two sets of
# components coincide and the distinction collapses. This is the case for
# [`CanonicalBasis`](@ref) and [`RotatedBasis`](@ref), and it is why tensors on
# those bases carry no variance tuple at all.

θ, ϕ, ψ = symbols("θ ϕ ψ", real = true)
ℬʳ = Basis(θ, ϕ, ψ)
typeof(ℬʳ)

#-

isorthonormal(ℬʳ), tsimplify(metric(ℬʳ, :cov))

# The same vector expressed on a rotated basis: covariant and contravariant
# components are equal.

W = Tens(Tensor{1, 3}(i -> symbols("w$i", real = true)))
tsimplify(components(W, ℬʳ, (:cov,)) - components(W, ℬʳ, (:cont,)))

# ## Numerically, on an orthogonal (scaled) basis
#
# An [`OrthogonalBasis`](@ref) is an orthonormal basis with a scaling factor
# ``\chi_i`` per direction — the situation of the *natural* basis of a
# curvilinear coordinate system, where the ``\chi_i`` are the Lamé coefficients
# (see [Curvilinear differential calculus](@ref th-curvilinear)). Its metric is
# diagonal, ``g_{ij}=\mathrm{diag}(\chi_i^2)``.

ℬ₀ = Basis(0.3, 0.7, 0.0)          # orthonormal, rotated
χ = (2.0, 3.0, 0.5)                 # scaling factors
ℬχ = Basis(ℬ₀, χ)
typeof(ℬχ)

#-

round.(metric(ℬχ, :cov), digits = 12)

# Diagonal, with entries ``\chi_i^2``:

round.(diag(metric(ℬχ, :cov)) .- [χ[i]^2 for i in 1:3], digits = 14)

# and the dual vectors are ``\underline{e}^i=\underline{e}_i/\chi_i^2``:

round.(
    vecbasis(ℬχ, :cont) .- vecbasis(ℬχ, :cov) * diagm([1 / χ[i]^2 for i in 1:3]),
    digits = 14,
)

# ## Summary
#
# | Basis type | metric | variance matters? |
# |:--|:--|:--|
# | [`CanonicalBasis`](@ref) | ``\boldsymbol{1}`` | no |
# | [`RotatedBasis`](@ref) | ``\boldsymbol{1}`` | no |
# | [`OrthogonalBasis`](@ref) | ``\mathrm{diag}(\chi_i^2)`` | yes, diagonally |
# | [`Basis`](@ref) | full | yes |
#
# The constructor always returns the most specific type that applies, so the
# cost of a component conversion is decided by the basis, not by the caller.
