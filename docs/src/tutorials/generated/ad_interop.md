```@meta
EditURL = "../../../../scripts/30_ad_interop.jl"
```

# Differentiation and interoperability

`TensND` is generic in its element type. The same code paths run on
`Float64`, on `ForwardDiff.Dual`, on SymPy's `Sym` and on Symbolics' `Num`,
which means derivatives can be taken by automatic differentiation *or*
symbolically, and results compared.

Theory: [Bases and variance](@ref th-bases-variance),
[Projection onto a symmetry class](@ref th-projection).

````@example ad_interop
using TensND
using LinearAlgebra
using ForwardDiff
using Printf
using Tensors

arr(x) = get_array(x)
````

## The four scalar worlds

The same tensor construction, four element types:

````@example ad_interop
for T in (Float64, Rational{Int})
    t = TensISO{3}(T(3), T(2))
    println(rpad(string(T), 18), " → ", typeof(t))
end
````

Symbolic, with SymPy:

````@example ad_interop
using SymPy
ks, μs = symbols("k μ", positive = true)
ℂˢ = TensISO{3}(3ks, 2μs)
typeof(ℂˢ)
````

and with Symbolics, the native Julia CAS:

````@example ad_interop
using Symbolics
Symbolics.@variables kn μn
ℂⁿ = TensISO{3}(3kn, 2μn)
typeof(ℂⁿ)
````

Both carry the same algebra. The inverse of an isotropic tensor is
``\tfrac{1}{3k}\mathbb{J}+\tfrac{1}{2\mu}\mathbb{K}`` in either:

````@example ad_interop
(tsimplify(get_data(inv(ℂˢ))), tsimplify(get_data(inv(ℂⁿ))))
````

## Derivatives through the tensor algebra

A quantity every homogenization scheme needs: the derivative of an effective
modulus with respect to a phase modulus. Here the simplest possible instance,
a dilute strain concentration
``\mathbb{A}=[\mathbb{I}+\mathbb{P}:(\mathbb{C}_1-\mathbb{C}_0)]^{-1}``, whose
whole computation is inversions and double contractions of `TensISO` objects.

````@example ad_interop
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
````

The derivative flows through `inv` and `⊡` on structured types without any
special handling: those methods are written on the stored coefficients, and
the coefficients are `Dual` numbers here.

## Differentiating a projection

More demanding: the derivative of a **projected** modulus with respect to an
Euler angle. Rotating a transversely isotropic tensor and projecting the
result back onto transverse isotropy about a *fixed* axis gives a quantity
that varies with the rotation.

````@example ad_interop
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
````

At ``\alpha=0`` the derivative vanishes: the tensor is already aligned with
the projection axis, so the fit is stationary — which is exactly the condition
the orientation optimizer of
[Projection onto a symmetry class](@ref th-projection) solves for.

The orthotropic projection differentiates just as well, in the canonical frame
and in a rotated one:

````@example ad_interop
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
````

## Symbolic and automatic differentiation agree

The same derivative, computed symbolically and numerically. Take the bulk
modulus of the dilute estimate as a function of the inclusion modulus:

````@example ad_interop
k₁s = symbols("k₁", positive = true)
f = 20 + 0.2 * (k₁s - 20) * (3 * 20 + 4 * 8) / (3 * k₁s + 4 * 8)
dfdk = diff(f, k₁s)

sym_value = Float64(dfdk.subs(k₁s, 60))
ad_value = ForwardDiff.derivative(k -> 20 + 0.2 * (k - 20) * (3 * 20 + 4 * 8) / (3k + 4 * 8), 60.0)
@printf "symbolic = %.12f   automatic = %.12f\n" sym_value ad_value
````

## Interoperability with `Tensors.jl`

`TensND` stores its data in `Tensors.jl` containers whenever the shape allows
it, so the two libraries compose. A `SymmetricTensor` goes in, and the
symmetry is preserved:

````@example ad_interop
st = SymmetricTensor{2, 3}((i, j) -> Float64(i + j))
t = Tens(st)
(typeof(t), typeof(get_array(t)))
````

and an order-4 minor-symmetric tensor round-trips through Kelvin–Mandel:

````@example ad_interop
C4 = SymmetricTensor{4, 3}((i, j, k, l) -> Float64(i + j + k + l))
tc = Tens(C4)
norm(arr(inv_KM(KM(tc))) - arr(tc))
````

## What is generic, and what is not

| Element type | tensor algebra | structured types | projections | symbolic operators | numerical operators |
|:--|:--|:--|:--|:--|:--|
| `Float64` | yes | yes | yes | — | yes |
| `ForwardDiff.Dual` | yes | yes | yes | — | yes |
| `SymPy.Sym` | yes | yes | fixed axis only | yes | — |
| `Symbolics.Num` | yes | yes | fixed axis only | partial | — |

The one genuine restriction is the **orientation search**: it is a numerical
optimization, so it needs a numeric element type. Projection onto a class with
a *given* axis or frame is closed-form and works symbolically.

````@example ad_interop
Bsym, dsym, _ = proj_tens(:TI, arr(TensISO{3}(3ks, 2μs)), [Sym(0), Sym(0), Sym(1)])
(typeof(Bsym), tsimplify(dsym))
````

An isotropic tensor is exactly transversely isotropic about any axis, so the
symbolic projection distance is identically zero.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

