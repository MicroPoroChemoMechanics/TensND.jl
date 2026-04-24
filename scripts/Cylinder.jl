# ============================================================================
#  Cylindrical coordinate system — symbolic symmetrised gradient
#
#  Builds an axisymmetric (independent of θ) displacement field
#    𝛏(r, z) = ξʳ(r, z) 𝐞ʳ + ξᶻ(r, z) 𝐞ᶻ
#  and computes its symmetric gradient (strain tensor) in the
#  normalised cylindrical basis.
# ============================================================================

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

using TensND, LinearAlgebra, SymPy, Tensors, OMEinsum, Rotations
sympy.init_printing(use_unicode = true)

# Cylindrical coordinate system (r, θ, z) with its natural frame
CS = coorsys_cylindrical()
r, θ, z = getcoords(CS)
𝐞ʳ, 𝐞ᶿ, 𝐞ᶻ = unitvec(CS)
ℬ = normalized_basis(CS)
@set_coorsys CS

# Axisymmetric displacement field 𝛏(r, z)
ξʳ, ξᶻ = SymFunction("ξʳ, ξᶻ", real = true)
𝛏 = ξʳ(r, z) * 𝐞ʳ + ξᶻ(r, z) * 𝐞ᶻ

# Symmetric gradient → 2nd-order strain tensor
𝛜 = SYMGRAD(𝛏)
