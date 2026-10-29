module TensND

import Base: eltype
import LinearAlgebra: normalize, dot, tr

using LinearAlgebra, Tensors, Rotations
# `divergence` is Tensors.jl's function, extended to the coordinate systems: one
# function, not two that the same `using` would make ambiguous.
import Tensors: divergence
using ForwardDiff, StaticArrays
using Random: Random, AbstractRNG
using SymPyCore: SymPyCore, Sym
using Symbolics

include("contraction_kernels.jl")
include("array_utils.jl")
include("bases.jl")
include("tens.jl")
include("tens_display.jl")
include("tens_components.jl")
include("tensors_bridge.jl")
include("kelvin_mandel.jl")
include("tens_algebra.jl")
include("tens_isotropic.jl")
include("tens_anisotropic.jl")
include("tens_cubic.jl")
include("structured_tens_ops.jl")
include("structured_tens_promotion.jl")
include("symmetric_parts.jl")
include("tens_projection.jl")
include("orientation_search.jl")
include("symmetry_search.jl")
include("tens_rotational_average.jl")
include("km_blocks.jl")
include("special_tens.jl")
include("random_tens.jl")
include("coorsystems.jl")
include("coorsystems_num.jl")
include("submanifold.jl")
include("deprecated.jl")

# The constructors whose default coordinates are SymPy symbols get them from the
# extension that SymPyPythonCall loads. Called without coordinates while it is
# not loaded, they would end in a bare `MethodError`: the hint says what to do.
const _SYMPY_DEFAULTS = (
    coorsys_cartesian, coorsys_polar, coorsys_cylindrical, coorsys_spherical, coorsys_spheroidal,
    init_cartesian, init_polar, init_cylindrical, init_spherical, init_rotated,
)

function _sympy_defaults_hint(io, exc, argtypes, kwargs)
    exc.f in _SYMPY_DEFAULTS || return nothing
    Base.get_extension(TensND, :TensNDSymPyPythonCallExt) === nothing || return nothing
    print(
        io, "\nThe default coordinates of `", exc.f, "` are SymPy symbols: load SymPyPythonCall ",
        "(`using SymPyPythonCall`), or pass the coordinates (for instance from Symbolics' `@variables`)."
    )
    return nothing
end

__init__() = Base.Experimental.register_error_hint(_sympy_defaults_hint, MethodError)

end # module
