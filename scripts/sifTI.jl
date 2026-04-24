# ============================================================================
#  Stress intensity factors — transversely-isotropic matrix
#
#  Legacy research sandbox that predates the in-library Walpole support.
#  Re-defines the Walpole basis locally (`Walpole_Basis`, `Walpole_Basis_sym`)
#  and the TI 4th-order constructor (`defTI`, `defcompTI`) to experiment
#  with SIF integrals for a TI medium.  The in-library equivalents are
#  `walpole_basis` / `walpole_basis_sym` in TensND and `tens_TI` /
#  `tens_TI_eng` / `tens_TI_Hoenig` — prefer those in new code.
# ============================================================================

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

using TensND, LinearAlgebra, SymPy, Tensors, OMEinsum, Rotations, Latexify
sympy.init_printing(use_unicode = true)

𝕀, 𝕁, 𝕂 = iso_projectors(Val(3), Val(Sym))
𝟏 = tens_Id2(Val(3), Val(Sym))

E, k, μ = symbols("E k μ", positive = true)
ν = symbols("ν", real = true)
k = E / (3(1-2ν)) ; μ = E / (2(1+ν))
λ = k -2μ/3

function Walpole_Basis(𝐧)
    T = eltype(𝐧)
    𝟏 = tens_Id2(3, T)
    𝐩 = 𝐧⊗𝐧 ; 𝐪=𝟏-𝐩
    return 𝐩⊗𝐩, 𝐪⊗𝐪/2, 𝐩⊗𝐪/√(T(2)), 𝐪⊗𝐩/√(T(2)), 𝐪⊠ˢ𝐪-𝐪⊗𝐪/2, 𝐪⊠ˢ𝐩+𝐩⊠ˢ𝐪
end

function Walpole_Basis_sym(𝐧)
    T = eltype(𝐧)
    𝟏 = tens_Id2(3, T)
    𝐩 = 𝐧⊗𝐧 ; 𝐪=𝟏-𝐩
    return 𝐩⊗𝐩, 𝐪⊗𝐪/2, 𝐩⊗𝐪/√(T(2))+𝐪⊗𝐩/√(T(2)), 𝐪⊠ˢ𝐪-𝐪⊗𝐪/2, 𝐪⊠ˢ𝐩+𝐩⊠ˢ𝐪
end

function Walpole(t::AbstractTens,𝐄)
    return ntuple(i->t⊙𝐄[i]/(𝐄[i]⊙𝐄[i]),length(𝐄))
end

function Walpole(t,𝐄)
    return sum([τ*𝐄[i] for (i,τ) ∈ enumerate(t)])
end

𝐧 = 𝐞(3)

𝐄ˢ = Walpole_Basis_sym(𝐧)

w = Walpole(𝕂,𝐄ˢ)

function defTI(T₁₁₁₁, T₁₁₂₂, T₁₁₃₃, T₃₃₃₃, T₂₃₂₃, ṉ = tens_basis(CanonicalBasis{3,typeof(T₁₁₁₁)}(), 3))
    T = eltype(ṉ)
    𝕎 = Walpole_Basis_sym(ṉ)
    𝕋 = T₃₃₃₃ * 𝕎[1] + (T₁₁₁₁ + T₁₁₂₂) * 𝕎[2] + √(T(2)) * T₁₁₃₃ * 𝕎[3] + (T₁₁₁₁ - T₁₁₂₂) * 𝕎[4] + 2 * T₂₃₂₃ * 𝕎[5]
    return 𝕋
end

function defcompTI(E₁, E₃, ν₁₂, ν₃₁, G₃₁, ṉ = tens_basis(CanonicalBasis{3,typeof(E₁)}(), 3))
    T = eltype(ṉ)
    𝕎 = Walpole_Basis_sym(ṉ)
    T₁₁₁₁ = inv(E₁)
    T₃₃₃₃ = inv(E₃)
    T₁₁₂₂ = -ν₁₂ / E₁
    T₁₁₃₃ = -ν₃₁ / E₃
    T₂₃₂₃ = inv(4G₃₁)
    𝕋 = T₃₃₃₃ * 𝕎[1] + (T₁₁₁₁ + T₁₁₂₂) * 𝕎[2] + √(T(2)) * T₁₁₃₃ * 𝕎[3] + (T₁₁₁₁ - T₁₁₂₂) * 𝕎[4] + 2 * T₂₃₂₃ * 𝕎[5]
    return 𝕋
end

function arg_TI(ℂ, ṉ = tens_basis(get_basis(ℂ), 3))
    C = tens_basis(get_basis(ℂ), 3) == ṉ ? ℂ : change_tens(ℂ, Basis(angles(components_canon(ṉ))...))
    return C[1, 1, 1, 1], C[1, 1, 2, 2], C[1, 1, 3, 3], C[3, 3, 3, 3], C[2, 3, 2, 3]
end

function argcompTI(𝕊, ṉ = tens_basis(get_basis(𝕊), 3))
    S = tens_basis(get_basis(𝕊), 3) == ṉ ? 𝕊 : change_tens(𝕊, Basis(angles(components_canon(ṉ))...))
    E₁ = inv(S[1, 1, 1, 1])
    E₃ = inv(S[3, 3, 3, 3])
    ν₁₂ = -E₁ * S[1, 1, 2, 2]
    ν₃₁ = -E₃ * S[1, 1, 3, 3]
    G₃₁ = inv(4S[2, 3, 2, 3])
    return E₁, E₃, ν₁₂, ν₃₁, G₃₁
end

function B_TI(η::Sym, ℬ, C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃)
    T = Sym
    σᵞ² = (C₁₁₁₁ * C₃₃₃₃ - C₁₁₃₃^2 - 2C₁₁₃₃ * C₂₃₂₃) / (C₂₃₂₃ * C₃₃₃₃)
    πᵞ = √(C₁₁₁₁ / C₃₃₃₃)
    σᵞ = symbols("σᵞ", real = true)
    σᵞ = √(σᵞ² + 2πᵞ)
    R₁₁₁₁ = -1 / σᵞ * (1 / C₂₃₂₃ + 1 / √(C₁₁₁₁ * C₃₃₃₃))
    R₁₁₃₃ = 1 / σᵞ * (C₁₁₃₃ + C₂₃₂₃) / (C₂₃₂₃ * C₃₃₃₃)
    R₃₃₃₃ = 1 / σᵞ * (√(C₁₁₁₁ / C₃₃₃₃) - C₁₁₃₃ * (C₁₁₃₃ + 2C₂₃₂₃) / (C₂₃₂₃ * C₃₃₃₃)) / C₃₃₃₃
    R₂₃₂₃ = 1 / (4 * √(T(2))) * √((C₁₁₁₁ - C₁₁₂₂) / C₂₃₂₃^3)
    R₃₁₃₁ = 1 / (4σᵞ) * (C₁₁₁₁ * C₃₃₃₃ - C₁₁₃₃^2) / (C₂₃₂₃^2 * C₃₃₃₃)
    η² = η^2
    k² = 1 - η²
    ℰ, 𝒦 = symbols("ℰ, 𝒦", real = true)
    ℰ = sympy.elliptic_e(k²)
    𝒦 = sympy.elliptic_k(k²)
    I₀ = ℰ / η²
    I₂ = k² == 0 ? PI / 4 : (𝒦 - ℰ) / k²
    Bₙₙ = 4 / (3η) / (I₀ * (C₁₁₃₃^2 * R₁₁₁₁ + 2C₁₁₃₃ * C₃₃₃₃ * R₁₁₃₃ + C₃₃₃₃^2 * R₃₃₃₃))
    Bₗₗ = 1 / (3η) / (C₂₃₂₃^2 * ((I₀ - I₂) * R₂₃₂₃ + I₂ * R₃₁₃₁))
    Bₘₘ = 1 / (3η) / (C₂₃₂₃^2 * ((I₀ - I₂) * R₃₁₃₁ + I₂ * R₂₃₂₃))
    𝐁 = Tens(Diagonal([Bₗₗ, Bₘₘ, Bₙₙ]), ℬ)
    ṉ = tens_basis(ℬ, 3)
    ℍ = 3 / 4 * ṉ ⊗ˢ 𝐁 ⊗ˢ ṉ
    return 𝐁, ℍ
end

C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃ = symbols("C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃", real = true)

η, a = symbols("η, a", positive = true)

ℬ = Basis()
𝐁, ℍ = B_TI(η::Sym, ℬ, C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃)

θ = symbols("θ", positive = true)
d = √(η^2*cos(θ)^2+sin(θ)^2)
𝛎 = (η*cos(θ) * 𝐞(1) + sin(θ) * 𝐞(2))/d
𝛕 = (-sin(θ) * 𝐞(1) + η*cos(θ) * 𝐞(2))/d

ℬ̄ = Basis(hcat(𝛎,𝛕,𝐧))

for 𝐯 ∈ (:𝛕,:𝛎) 
    @eval $(𝐯) = tsimplify(change_tens($(𝐯),ℬ̄))
end

Bₙₙᶜʸˡ = 3PI/8 * limit(𝐁[3,3]/η, η=>0)
Bₘₘᶜʸˡ = 3PI/8 * limit(𝐁[2,2]/η, η=>0)
Bₗₗᶜʸˡ = 3PI/8 * limit(𝐁[1,1]/η, η=>0)

𝐁ᶜʸˡ =  Bₗₗᶜʸˡ * 𝛕⊗𝛕 + Bₘₘᶜʸˡ * 𝛎⊗𝛎 + Bₙₙᶜʸˡ * 𝐧⊗𝐧 
𝐁ᶜʸˡ⁻¹ =  (1/Bₗₗᶜʸˡ) * 𝛕⊗𝛕 + (1/Bₘₘᶜʸˡ) * 𝛎⊗𝛎 + (1/Bₙₙᶜʸˡ) * 𝐧⊗𝐧 

p, q, ω = symbols("p, q, ω", real = true)
𝐓 = q * (cos(ω)*𝐞(1) + sin(ω)*𝐞(2)) + p * 𝐧

𝐊 = 3*PI^(3//2)/8 *√(a*d/η) * 𝐁ᶜʸˡ⁻¹⋅𝐁⋅𝐓

Kᴵ = 𝐊[3]

ℂ = defTI(C₁₁₁₁, C₁₁₂₂, C₁₁₃₃, C₃₃₃₃, C₂₃₂₃)

Ciso = 5. * 𝕁 + 5. * 𝕂 + 0.0001*𝐄ˢ[1]
valnum = (C₁₁₁₁=>Ciso[1,1,1,1], C₁₁₂₂=>Ciso[1,1,2,2], C₁₁₃₃=>Ciso[1,1,3,3], C₃₃₃₃=>Ciso[3,3,3,3], C₂₃₂₃=>Ciso[2,3,2,3],θ=>0*π/3,η=>0.999999,a=>1.,ω=>0*π/6,q=>1.)

valnum = (C₁₁₁₁=>5.8, C₁₁₂₂=>1., C₁₁₃₃=>2.1, C₃₃₃₃=>2.8, C₂₃₂₃=>1.3,θ=>π/5,η=>0.2,a=>0.5,ω=>π/3,p=>100.,q=>10000.)

𝐊ⁿᶸᵐ = tsimplify(tsubs(𝐊,valnum...).evalf())

A = C₁₁₁₁*C₂₃₂₃ ; B = C₁₁₃₃*(C₁₁₃₃+2C₂₃₂₃)-C₁₁₁₁*C₃₃₃₃ ; C = C₃₃₃₃*C₂₃₂₃
Δ = B^2-4A*C ; n₁ = (-B+√Δ)/(2A) ; n₂ = (-B-√Δ)/(2A)
m(n) = (C₁₁₁₁*n-C₂₃₂₃)/(C₁₁₃₃+C₂₃₂₃)
n₃ = 2C₂₃₂₃/(C₁₁₁₁-C₁₁₂₂)
mm(n) = (C₁₁₃₃+C₂₃₂₃)*n/(C₃₃₃₃-C₂₃₂₃*n)
m₁ = m(n₁) ; m₂ = m(n₂)

α₂ = (1+m₁)*(1+m₂)/(m₁-m₂)*(1/√n₁-1/√n₂) ; β₂ = -1/√n₃

b = a*η
η² = η^2
k² = 1 - η²
ℰ = sympy.elliptic_e(k²)
𝒦 = sympy.elliptic_k(k²)

A₁ = a*b^2*k²*q*cos(ω)/((k²-1+β₂/α₂)*ℰ+(1-β₂/α₂)*η²*𝒦)/(2α₂*C₂₃₂₃)
B₁ = a*b^2*k²*q*sin(ω)/((k²+(1-β₂/α₂)*η²)*ℰ-(1-β₂/α₂)*η²*𝒦)/(2α₂*C₂₃₂₃)

k₂ = 2α₂*C₂₃₂₃*(a*B₁*sin(θ)+b*A₁*cos(θ))/(a*b)^(3/2)/√(a*d)*√(PI) ;
k₃ = -2β₂*C₂₃₂₃*(a*A₁*sin(θ)-b*B₁*cos(θ))/(a*b)^(3/2)/√(a*d)*√(PI) ;
k₁ = p/ℰ*√(b*d)*√(PI) ;

real(simplify(subs(k₂,valnum...).evalf()))
real(simplify(subs(k₃,valnum...).evalf()))
simplify(subs(k₁,valnum...).evalf())

𝐊ⁿᶸᵐ
