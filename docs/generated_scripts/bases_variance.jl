import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "docs"); io = devnull)
# The `docs` environment declares every dependency the tutorials use
# (SymPy, Symbolics, NLopt, BenchmarkTools, ...) and resolves TensND
# itself through `[sources]`. These four lines are stripped from the
# generated page and notebook, which already run inside that project.

using TensND
using LinearAlgebra
using SymPy
using Tensors

ℬ = Basis(Sym[1 0 0; 0 1 0; 0 1 1])

vecbasis(ℬ, :cov)

vecbasis(ℬ, :cont)

metric(ℬ, :cov)

metric(ℬ, :cont)

tsimplify(metric(ℬ, :cont) * metric(ℬ, :cov))

tsimplify(transpose(vecbasis(ℬ, :cont)) * vecbasis(ℬ, :cov))

isorthogonal(ℬ), isorthonormal(ℬ)

V = Tens(Tensor{1, 3}(i -> symbols("v$i", real = true)))

components(V, ℬ, (:cont,))

components(V, ℬ, (:cov,))

tsimplify(metric(ℬ, :cov) * components(V, ℬ, (:cont,)) - components(V, ℬ, (:cov,)))

T = Tens(Tensor{2, 3}((i, j) -> symbols("t$i$j", real = true)))

components(T, ℬ, (:cov, :cov))

tfactor(tsimplify(components(T, ℬ, (:cont, :cov))))

G = Tens(metric(ℬ, :cov), ℬ, (:cov, :cov))
components(G, (:cont, :cov))

ℬ̄ = normalize(ℬ)
metric(ℬ̄, :cov)

isorthogonal(ℬ̄), isorthonormal(ℬ̄)

θ, ϕ, ψ = symbols("θ ϕ ψ", real = true)
ℬʳ = Basis(θ, ϕ, ψ)
typeof(ℬʳ)

isorthonormal(ℬʳ), tsimplify(metric(ℬʳ, :cov))

W = Tens(Tensor{1, 3}(i -> symbols("w$i", real = true)))
tsimplify(components(W, ℬʳ, (:cov,)) - components(W, ℬʳ, (:cont,)))

ℬ₀ = Basis(0.3, 0.7, 0.0)          # orthonormal, rotated
χ = (2.0, 3.0, 0.5)                 # scaling factors
ℬχ = Basis(ℬ₀, χ)
typeof(ℬχ)

round.(metric(ℬχ, :cov), digits = 12)

round.(diag(metric(ℬχ, :cov)) .- [χ[i]^2 for i in 1:3], digits = 14)

round.(
    vecbasis(ℬχ, :cont) .- vecbasis(ℬχ, :cov) * diagm([1 / χ[i]^2 for i in 1:3]),
    digits = 14,
)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
