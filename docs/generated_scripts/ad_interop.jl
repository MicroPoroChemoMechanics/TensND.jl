import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "docs"); io = devnull)
# The `docs` environment declares every dependency the tutorials use
# (SymPy, Symbolics, NLopt, BenchmarkTools, ...) and resolves TensND
# itself through `[sources]`. These four lines are stripped from the
# generated page and notebook, which already run inside that project.

using TensND
using LinearAlgebra
using ForwardDiff
using Printf
using Tensors

arr(x) = get_array(x)

for T in (Float64, Rational{Int})
    t = TensISO{3}(T(3), T(2))
    println(rpad(string(T), 18), " → ", typeof(t))
end

using SymPy
ks, μs = symbols("k μ", positive = true)
ℂˢ = TensISO{3}(3ks, 2μs)
typeof(ℂˢ)

using Symbolics
Symbolics.@variables kn μn
ℂⁿ = TensISO{3}(3kn, 2μn)
typeof(ℂⁿ)

(tsimplify(get_data(inv(ℂˢ))), tsimplify(get_data(inv(ℂⁿ))))

𝕀 = tens_Id4(Val(3), Val(Float64))
𝕁 = tens_J4(Val(3), Val(Float64))
𝕂 = tens_K4(Val(3), Val(Float64))

function dilute_k(k₁)
    ℂ₀ = TensISO{3}(3 * 20.0, 2 * 8.0)
    ℂ₁ = TensISO{3}(3 * k₁, 2 * 30.0)
    # Hill tensor of a sphere in an isotropic matrix
    k₀, μ₀ = 20.0, 8.0
    ℙ = TensISO{3}(1 / (3k₀ + 4μ₀), (3 * (k₀ + 2μ₀)) / (5μ₀ * (3k₀ + 4μ₀)))
    𝔸 = inv(𝕀 + ℙ ⊡ (ℂ₁ - ℂ₀))
    return get_data(ℂ₀ + 0.2 * ((ℂ₁ - ℂ₀) ⊡ 𝔸))[1] / 3
end

k_ad = ForwardDiff.derivative(dilute_k, 60.0)
k_fd = (dilute_k(60.0 + 1.0e-6) - dilute_k(60.0 - 1.0e-6)) / 2.0e-6
@printf "d k_eff / d k₁ : AD = %.12f   finite differences = %.12f\n" k_ad k_fd

n_fixed = [0.0, 0.0, 1.0]

function projected_ℓ₁(α)
    C = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, [sin(α), 0.0, cos(α)])
    B, _, _ = proj_tens(:TI, arr(C), n_fixed)
    return get_data(B)[1]
end

for α in (0.0, 0.3, 0.7)
    ad = ForwardDiff.derivative(projected_ℓ₁, α)
    fd = (projected_ℓ₁(α + 1.0e-6) - projected_ℓ₁(α - 1.0e-6)) / 2.0e-6
    @printf "  α = %.1f :  dℓ₁/dα  AD = %+.9f   FD = %+.9f\n" α ad fd
end

frame_rot = Basis(0.3, 0.8, 0.2)

function projected_C₁₁(α, frame)
    C = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, [sin(α), 0.0, cos(α)])
    B, _, _ = proj_tens(:ORTHO, arr(C), frame)
    return get_data(B)[1]
end

for (name, fr) in ("canonical" => CanonicalBasis{3, Float64}(), "rotated" => frame_rot)
    ad = ForwardDiff.derivative(a -> projected_C₁₁(a, fr), 0.5)
    fd = (projected_C₁₁(0.5 + 1.0e-6, fr) - projected_C₁₁(0.5 - 1.0e-6, fr)) / 2.0e-6
    @printf "  %-10s frame :  dC₁₁/dα  AD = %+.9f   FD = %+.9f\n" name ad fd
end

k₁s = symbols("k₁", positive = true)
f = 20 + 0.2 * (k₁s - 20) * (3 * 20 + 4 * 8) / (3 * k₁s + 4 * 8)
dfdk = diff(f, k₁s)

sym_value = Float64(dfdk.subs(k₁s, 60))
ad_value = ForwardDiff.derivative(k -> 20 + 0.2 * (k - 20) * (3 * 20 + 4 * 8) / (3k + 4 * 8), 60.0)
@printf "symbolic = %.12f   automatic = %.12f\n" sym_value ad_value

st = SymmetricTensor{2, 3}((i, j) -> Float64(i + j))
t = Tens(st)
(typeof(t), typeof(get_array(t)))

C4 = SymmetricTensor{4, 3}((i, j, k, l) -> Float64(i + j + k + l))
tc = Tens(C4)
norm(arr(inv_KM(KM(tc))) - arr(tc))

Bsym, dsym, _ = proj_tens(:TI, arr(TensISO{3}(3ks, 2μs)), [Sym(0), Sym(0), Sym(1)])
(typeof(Bsym), tsimplify(dsym))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
