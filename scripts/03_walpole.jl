# # The Walpole basis end to end
#
# Everything the transversely isotropic machinery rests on, computed rather than
# quoted: the six ``\mathbb{W}_i`` as ``6\times6`` matrices, the multiplication
# table, the Gram matrix, the correct decomposition of ``\mathbb{I},\mathbb{J},
# \mathbb{K}``, the synthetic ``2\times2`` product and inverse rules, and the
# three parametrizations of a transversely isotropic material compared on the
# **same** material.
#
# Theory: [The Walpole basis](@ref th-walpole) and
# [TI parametrizations](@ref th-ti-parametrizations).

import Pkg                                                          #jl
Pkg.activate(joinpath(@__DIR__, "..", "docs"); io = devnull)        #jl
## The `docs` environment declares every dependency the tutorials use     #jl
## (SymPy, Symbolics, NLopt, BenchmarkTools, ...) and resolves TensND      #jl
## itself through `[sources]`. These four lines are stripped from the      #jl
## generated page and notebook, which already run inside that project.     #jl

using TensND
using LinearAlgebra
using Printf

n = (0.0, 0.0, 1.0)          # symmetry axis
W = walpole_basis(n)
arr(x) = get_array(x)

# ## The six tensors, in Kelvin–Mandel form
#
# Built from ``\boldsymbol{1}_n=\underline{n}\otimes\underline{n}`` and
# ``\boldsymbol{1}_T=\boldsymbol{1}-\boldsymbol{1}_n``, in the frame where
# ``\underline{n}=\underline{e}_3``.

function show_KM(label, t)
    M = KM(t)
    println(label)
    for r in 1:6
        println("   ", join([abs(M[r, c]) < 1.0e-12 ? "   .   " : @sprintf("%7.4f", M[r, c]) for c in 1:6], " "))
    end
    return println()
end

for i in 1:6
    show_KM("𝕎$i =", W[i])
end

# ``\mathbb{W}_3`` and ``\mathbb{W}_4`` are transposes of each other — the only
# two that are not major-symmetric individually:

norm(arr(W[3]) - permutedims(arr(W[4]), (3, 4, 1, 2)))

# ## The multiplication table
#
# The basis is closed under `⊡`. Decomposing each product back onto the basis
# shows the structure: ``(\mathbb{W}_1,\mathbb{W}_2,\mathbb{W}_3,\mathbb{W}_4)``
# behave exactly as the matrix units ``(E_{11},E_{22},E_{12},E_{21})`` of
# ``2\times2`` matrices, while ``\mathbb{W}_5`` and ``\mathbb{W}_6`` are
# orthogonal idempotents.

gram = [sum(arr(W[i]) .* arr(W[j])) for i in 1:6, j in 1:6]
coeffs(P) = [sum(P .* arr(W[k])) / gram[k, k] for k in 1:6]

println("   𝕎ᵢ : 𝕎ⱼ")
println("      ", join([lpad("𝕎$j", 5) for j in 1:6]))
for i in 1:6
    row = String[]
    for j in 1:6
        P = arr(W[i]) ⊡ arr(W[j])
        if norm(P) < 1.0e-12
            push!(row, lpad("0", 5))
        else
            c = coeffs(P)
            k = argmax(abs.(c))
            push!(row, lpad("𝕎$k", 5))
        end
    end
    println("  𝕎$i  ", join(row))
end

# ## The Gram matrix is diagonal
#
# ``\langle\mathbb{W}_i,\mathbb{W}_j\rangle=g_i\delta_{ij}`` with
# ``(g_i)=(1,1,1,1,2,2)``. This is the single fact that makes every TI
# projection closed-form: the normal equations decouple.

round.(gram, digits = 12)

# ## The isotropic tensors on this basis
#
# Easy to get wrong, so worth computing. ``\mathbb{I}`` is **not**
# ``\sum_i\mathbb{W}_i``, and ``\mathbb{J}`` is **not**
# ``\mathbb{W}_1+\mathbb{W}_2``.

𝕀, 𝕁, 𝕂 = ISO(Val(3), Val(Float64))
s2 = sqrt(2)

checks = [
    "𝕀 = 𝕎₁+𝕎₂+𝕎₅+𝕎₆" =>
        norm(arr(𝕀) - (arr(W[1]) + arr(W[2]) + arr(W[5]) + arr(W[6]))),
    "𝕁 = (𝕎₁+2𝕎₂+√2𝕎₃+√2𝕎₄)/3" =>
        norm(arr(𝕁) - (arr(W[1]) + 2arr(W[2]) + s2 * arr(W[3]) + s2 * arr(W[4])) / 3),
    "𝕂 = (2𝕎₁+𝕎₂−√2𝕎₃−√2𝕎₄)/3 + 𝕎₅ + 𝕎₆" =>
        norm(arr(𝕂) - ((2arr(W[1]) + arr(W[2]) - s2 * (arr(W[3]) + arr(W[4]))) / 3 + arr(W[5]) + arr(W[6]))),
    "— and the WRONG one: 𝕀 = Σᵢ𝕎ᵢ" =>
        norm(arr(𝕀) - sum(arr(W[i]) for i in 1:6)),
    "— and the WRONG one: 𝕁 = 𝕎₁+𝕎₂" =>
        norm(arr(𝕁) - (arr(W[1]) + arr(W[2]))),
]
for (name, r) in checks
    @printf "%-40s residual = %.3e\n" name r
end

# The last two are off by ``\sqrt2`` and ``1`` respectively — not round-off.
#
# Reading the Walpole coefficients of the three directly:

for (name, t) in (("𝕀", 𝕀), ("𝕁", 𝕁), ("𝕂", 𝕂))
    B, _, _ = proj_tens(:TI, arr(t), [0.0, 0.0, 1.0])
    println(rpad(name, 3), " → (ℓ₁,…,ℓ₆) = ", round.(get_ℓ(B), digits = 6))
end

# ## The synthetic algebra
#
# With ``L=\begin{pmatrix}\ell_1&\ell_3\\\ell_4&\ell_2\end{pmatrix}``, the double
# contraction is ``(L\,M,\ \ell_5m_5,\ \ell_6m_6)`` and the inverse is
# ``(L^{-1},\ 1/\ell_5,\ 1/\ell_6)``.

A = TensTI{4}(2.0, 1.0, 0.5, 0.9, 1.3, n)      # major-symmetric, N=5
B = TensTI{4}(3.0, 1.5, 0.2, 1.1, 0.4, n)
P = A ⊡ B

println("typeof(A)   = ", typeof(A))
println("typeof(A⊡B) = ", typeof(P))

# The product of two **major-symmetric** TI tensors is generally not
# major-symmetric — ``L\,M\neq{}^{t}(L\,M)`` unless the blocks commute — so the
# storage widens from `N=5` to `N=6`:

ℓA, ℓB, ℓP = get_ℓ(A), get_ℓ(B), get_ℓ(P)
LA = [ℓA[1] ℓA[3]; ℓA[4] ℓA[2]]
LB = [ℓB[1] ℓB[3]; ℓB[4] ℓB[2]]

println("L_A · L_B   = ", round.(LA * LB, digits = 10))
println("block of P  = ", round.([ℓP[1] ℓP[3]; ℓP[4] ℓP[2]], digits = 10))
println("ℓ₅ rule     : ", ℓA[5] * ℓB[5], " vs ", ℓP[5])
println("ℓ₆ rule     : ", ℓA[6] * ℓB[6], " vs ", ℓP[6])
println("ℓ₃ == ℓ₄ ?  ", isapprox(ℓP[3], ℓP[4]))

# Against the dense 81-component product:

norm(arr(P) - arr(A) ⊡ arr(B))

# Inversion, likewise closed-form:

norm(arr(A ⊡ inv(A)) - arr(𝕀))

# ## The symmetrized basis
#
# Merging the pair, ``\mathbb{W}^s_3=\mathbb{W}_3+\mathbb{W}_4``, spans exactly
# the five-dimensional major-symmetric subspace, with Gram ``(1,1,2,2,2)``:

Ws = walpole_basis_sym(n)
[
    norm(arr(Ws[1]) - arr(W[1])),
    norm(arr(Ws[2]) - arr(W[2])),
    norm(arr(Ws[3]) - (arr(W[3]) + arr(W[4]))),
    norm(arr(Ws[4]) - arr(W[5])),
    norm(arr(Ws[5]) - arr(W[6])),
]

#-

[sum(arr(Ws[i]) .* arr(Ws[i])) for i in 1:5]

# ## Three parametrizations, one material
#
# A unidirectional carbon/epoxy ply — a standard transversely isotropic
# composite, the fiber direction being the symmetry axis. Engineering constants
# first (moduli in GPa).

nv = [0.0, 0.0, 1.0]
E₁, E₃, ν₁₂, ν₃₁, G₃₁ = 9.0, 140.0, 0.40, 0.30, 4.6

𝕊 = tens_TI_eng(E₁, E₃, ν₁₂, ν₃₁, G₃₁, nv)     # compliance
ℂ = inv(𝕊)                                       # stiffness

println("engineering → back      : ", round.(arg_TI_eng(𝕊), digits = 8))
println("‖ℂ:𝕊 − 𝕀‖               : ", norm(arr(ℂ ⊡ 𝕊) - arr(𝕀)))

# The same compliance read in **component** form,
# ``(S_{1111},S_{1122},S_{1133},S_{3333},S_{2323})``:

round.(arg_TI(𝕊), sigdigits = 6)

# and in **Hoenig** form, the dimensionless ratios:

Eh, ν₁, ν₂, H, Γ = arg_TI_Hoenig(𝕊)
@printf "E  = %8.4f GPa\nν₁ = %8.4f\nν₂ = %8.4f\nH  = %8.4f   (axial/transverse modulus ratio)\nΓ  = %8.4f   (shear anisotropy)\n" Eh ν₁ ν₂ H Γ

# ``H`` and ``\Gamma`` quantify the anisotropy directly and dimensionlessly:
# here ``H\approx15.6``, so the fiber direction is some fifteen times stiffer
# than the transverse plane, while ``\Gamma\approx1.43`` says the axial shear
# modulus departs from its in-plane counterpart by rather less.
# Rebuilding from the Hoenig parameters returns the same tensor:

norm(arr(tens_TI_Hoenig(Eh, ν₁, ν₂, H, Γ, nv)) - arr(𝕊))

# An isotropic material sits at the point ``H=\Gamma=1`` with
# ``\nu_1=\nu_2`` — the Hoenig parameters measure departure from isotropy:

𝕊iso = inv(3 * 20.0 * 𝕁 + 2 * 8.0 * 𝕂)
Eᵢ, ν₁ᵢ, ν₂ᵢ, Hᵢ, Γᵢ = arg_TI_Hoenig(fromISO(𝕊iso, nv))
@printf "E = %.6f   ν₁ = %.6f   ν₂ = %.6f   H = %.10f   Γ = %.10f\n" Eᵢ ν₁ᵢ ν₂ᵢ Hᵢ Γᵢ

# ## The Walpole coefficients of the ply
#
# Stiffness and compliance are inverse in the synthetic algebra, so their
# ``2\times2`` blocks are inverse matrices and their ``\ell_5,\ell_6`` are
# reciprocal:

ℓC, ℓS = get_ℓ(ℂ), get_ℓ(𝕊)
LC = [ℓC[1] ℓC[3]; ℓC[4] ℓC[2]]
LS = [ℓS[1] ℓS[3]; ℓS[4] ℓS[2]]

println("L_ℂ · L_𝕊 = ", round.(LC * LS, digits = 12))
println("ℓ₅ᶜ·ℓ₅ˢ   = ", round(ℓC[5] * ℓS[5], digits = 12))
println("ℓ₆ᶜ·ℓ₆ˢ   = ", round(ℓC[6] * ℓS[6], digits = 12))
