using TensND
using Test
using TimerOutputs
using LinearAlgebra, SymPy, Tensors, OMEinsum, Rotations
using ForwardDiff
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
include("test_tens_anisotropic.jl")
include("test_tens_projection.jl")
include("test_special_tens.jl")
include("test_coorsystems.jl")
include("test_coorsystems_num.jl")
include("test_submanifold.jl")
# Pins the conventions the documentation asserts — the Walpole identities, the
# Kelvin-Mandel congruence, the operator index placement, and so on — so that
# the code cannot drift away from the docs unnoticed. See
# `docs/src/developer/testing_conventions.md`.
include("test_conventions.jl")
# Must stay last: it loads NLopt, which `test_tens_projection.jl` requires to
# be absent.
include("test_nlopt_ext.jl")

print_timer()
println()
