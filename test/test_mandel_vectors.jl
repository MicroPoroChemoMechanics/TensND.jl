# Symmetric tensors as Kelvin–Mandel vectors: each operation against Tensors.jl's
# own on the symmetric tensor, exact cases (uniaxial tension, a change of
# volume), derivatives through them, symbolic components, inference.

using TensND, Test, LinearAlgebra, Random, ForwardDiff, StaticArrays
import Symbolics
const _MV = TensND.Tensors
_mv_tensor(v) = _MV.frommandel(_MV.SymmetricTensor{2, 3}, v)

@testset "Mandel vectors against Tensors.jl" begin
    rng = Xoshiro(20261030)
    @test identity_mandel() isa Mandel{Float64} && identity_mandel(Float32) isa Mandel{Float32}
    @test _mv_tensor(identity_mandel()) == one(_MV.SymmetricTensor{2, 3})
    for _ in 1:50
        v = Mandel{Float64}(randn(rng, 6))
        t = _mv_tensor(v)
        @test trace_mandel(v) ≈ tr(t)
        @test _mv_tensor(deviator(v)) ≈ _MV.dev(t)
        @test von_mises(v) ≈ sqrt(3 / 2 * (_MV.dev(t) ⊡ _MV.dev(t)))
        @test equivalent_strain(v) ≈ sqrt(2 / 3 * (_MV.dev(t) ⊡ _MV.dev(t)))
        @test dot(v, v) ≈ t ⊡ t                       # the basis is orthonormal
        @test sqnorm_deviator(v) ≈ _MV.dev(t) ⊡ _MV.dev(t)
    end
    # Uniaxial tension and incompressible uniaxial strain.
    @test von_mises(Mandel{Float64}(-3.0, 0, 0, 0, 0, 0)) == 3.0
    @test equivalent_strain(Mandel{Float64}(2.0, -1.0, -1.0, 0, 0, 0)) ≈ 2.0
    # A change of volume has no deviatoric part, to the bit, for any amount.
    for a in (0.1, 1 / 3, -7.3e-5)
        @test von_mises(a * identity_mandel()) == 0
        @test sqnorm_deviator(a * identity_mandel()) == 0
    end
end

@testset "Mandel vectors: derivatives, symbols, inference" begin
    v = Mandel{Float64}(1.0, 2.0, -0.5, 0.3, -0.2, 0.7)
    # d tr / dv = 1 (the identity), d σeq / dσ = 3 s / (2 σeq).
    @test ForwardDiff.gradient(trace_mandel, v) == identity_mandel()
    @test ForwardDiff.gradient(von_mises, v) ≈ 3 * deviator(v) / (2 * von_mises(v))
    @test ForwardDiff.jacobian(deviator, v) ≈ I - identity_mandel() * transpose(identity_mandel()) / 3
    # Symbolic components, evaluated back.
    Symbolics.@variables a b c d e f
    vs = Mandel{Symbolics.Num}(a, b, c, d, e, f)
    vals = Dict(a => v[1], b => v[2], c => v[3], d => v[4], e => v[5], f => v[6])
    num(x) = Symbolics.value(Symbolics.substitute(x, vals; fold = Val(true)))
    @test num(trace_mandel(vs)) ≈ trace_mandel(v)
    @test num.(deviator(vs)) ≈ deviator(v)
    @test num(von_mises(vs)) ≈ von_mises(v)
    @test num(equivalent_strain(vs)) ≈ equivalent_strain(v)
    @test (@inferred deviator(v)) isa Mandel{Float64}
    @test (@inferred von_mises(v)) isa Float64
    @test (@inferred identity_mandel(Float32)) isa Mandel{Float32}
end
