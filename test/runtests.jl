using TensND
using Aqua
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

"""
    aqua_persistent_tasks(m::Module) -> Bool

`Aqua.test_persistent_tasks` in a child process with default bounds checking.

The check precompiles a synthetic package that depends on `m` and waits for a
sentinel written from that package's module body. Julia writes no
precompilation cache when `--check-bounds` is forced, so the body never runs,
the sentinel never appears, and the check fails for a reason that says nothing
about `m`. `Pkg.test()` forces the flag by default.

Running it in a child with `--check-bounds=auto` keeps both halves: the test
suite proper still runs with bounds checking on, and the check still runs.
"""
function aqua_persistent_tasks(m::Module)
    code = """
    using Aqua, $(nameof(m)), Test
    @testset "persistent_tasks" begin
        Aqua.test_persistent_tasks($(nameof(m)))
    end
    """
    cmd = `$(first(Base.julia_cmd())) --check-bounds=auto --startup-file=no --project=$(Base.active_project()) -e $code`
    return success(run(ignorestatus(cmd)))
end

@testsection "Aqua" begin
    Aqua.test_all(
        TensND;
        # `array_utils.jl` extends `Tensors` and `LinearAlgebra` functions to
        # `Tensors.Tensor`, `SymmetricTensor` and plain `AbstractArray`s. Both
        # the function and the type belong elsewhere, so it is type piracy by
        # definition, and it is deliberate: the float-tolerant predicates exist
        # because the upstream ones compare with `==`, useless on `Float64`,
        # `Dual`, `Num` or `Sym`, and the array-level operators are what let the
        # tensor algebra accept ordinary arrays. Declared rather than hidden.
        #
        # What the declaration does NOT fix: `Tensors.otimes(::AbstractArray,
        # ::AbstractArray)` applies to any array as soon as TensND is loaded,
        # including in code that never asked for it.
        piracies = (
            treat_as_own = [
                # `FourthOrderTensor` is an alias for a three-member union;
                # `is_foreign` on a union is an `any`, so every member has to be
                # named or the whole alias counts as foreign.
                Tensors.Tensor, Tensors.SymmetricTensor, Tensors.MixedTensor,
                Tensors.AbstractTensor,
                Tensors.ismajorsymmetric, Tensors.isminorsymmetric,
                Tensors.majortranspose, LinearAlgebra.issymmetric,
                Tensors.otimes, Tensors.otimesu, Tensors.otimesl,
                Tensors.dcontract, Tensors.dotdot,
                # `tens_isotropic.jl` also defines `+` and `*` between a
                # `UniformScaling` and a `Tensors.AbstractTensor`, which is what
                # makes `I + t` and `2I * t` work on tensors.
                LinearAlgebra.UniformScaling,
            ],
        ),
        # The twelve ambiguities left are all crossings with third-party
        # signatures: Symbolics' `/(::AbstractArray{<:Real}, ::Num)`,
        # MutableArithmetics' `dot` on mutable element types, and
        # `LinearAlgebra.dot` on a `Diagonal`. Resolving each would mean adding a
        # method that decides what, say, `dot(::TensISO{2}, ::Diagonal)` returns
        # — a modeling question, not a dispatch one. The package's own
        # ambiguities, seventeen of them, were fixed rather than excluded.
        ambiguities = (exclude = [LinearAlgebra.dot, LinearAlgebra.normalize, /],),
        # Run apart, see `aqua_persistent_tasks` above.
        persistent_tasks = false,
    )
    @test aqua_persistent_tasks(TensND)
end

opequal(x, y) = x == y || x ≈ y


include("test_bases.jl")
include("test_tens.jl")
include("test_tensor_products.jl")
include("test_tens_isotropic.jl")
include("test_tens_anisotropic.jl")
include("test_tens_projection.jl")
include("test_tens_rotational_average.jl")
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
