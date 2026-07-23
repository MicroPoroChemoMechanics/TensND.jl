using TensND
using Test
using TimerOutputs
using LinearAlgebra, SymPy, Tensors, OMEinsum, Rotations
using Random

# `test_bases.jl` and `test_tens_isotropic.jl` draw random bases and moduli.
# Seed once so a CI failure is always reproducible locally instead of
# depending on the draw.
Random.seed!(20260723)

macro testsection(str, block)
    return quote
        @timeit "$($(esc(str)))" begin
            @testset "$($(esc(str)))" begin
                $(esc(block))
            end
        end
    end
end

reset_timer!()

opequal(x, y) = x == y || x ≈ y


include("test_bases.jl")
include("test_tens.jl")
include("test_tens_isotropic.jl")
include("test_tens_walpole.jl")
include("test_tens_projection.jl")
include("test_special_tens.jl")
include("test_coorsystems.jl")
include("test_coorsystems_num.jl")

print_timer()
println()
