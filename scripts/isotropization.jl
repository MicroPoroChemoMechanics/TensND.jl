# ============================================================================
#  Isotropic projections of a generic 4th-order tensor
#
#  Uses the J/K projectors to extract the isotropic (bulk + shear)
#  part of a symbolic 4th-order tensor ℂ:
#     k = ⟨ℂ, 𝕁⟩ / 3       (bulk modulus)
#     μ = ⟨ℂ, 𝕂⟩ / 10      (shear modulus, dim=3)
# ============================================================================

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

using TensND, LinearAlgebra, SymPy, Tensors, OMEinsum, Rotations
sympy.init_printing(use_unicode = true)

# Symmetric-identity, spherical, deviatoric projectors (symbolic, dim=3)
𝕀, 𝕁, 𝕂 = iso_projectors(Val(3), Val(Sym))

# --- 1.  Analytic isotropic tensor  ℂ = 3k𝕁 + 2μ𝕂  ----------------------------
k, μ = symbols("k μ", positive = true)
ℂ = 3k * 𝕁 + 2μ * 𝕂

# --- 2.  Isotropic projection of a fully symbolic ℂᵢⱼₖₗ ----------------------
ℂ = Tens(SymmetricTensor{4, 3}((i, j, k, l) -> symbols("C$i$j$k$l", real = true)))

μ = simplify((ℂ ⊙ 𝕂) / 10)        # shear modulus
k = (ℂ ⊙ 𝕁) / 3                    # bulk modulus
λ = k - 2μ / 3                     # Lamé first parameter
