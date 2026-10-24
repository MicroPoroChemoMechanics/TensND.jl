module TensND

import Base: eltype
import LinearAlgebra: normalize, dot, tr

using LinearAlgebra, Tensors, Rotations
# `divergence` is Tensors.jl's function, extended to the coordinate systems: one
# function, not two that the same `using` would make ambiguous.
import Tensors: divergence
using ForwardDiff, StaticArrays
using SymPy
using Symbolics

include("sympy_matrix.jl")
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
include("coorsystems.jl")
include("coorsystems_num.jl")
include("submanifold.jl")
include("deprecated.jl")

end # module
