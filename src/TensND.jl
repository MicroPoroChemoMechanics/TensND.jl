module TensND

import Base: @pure, eltype
import LinearAlgebra: normalize, dot, tr

using LinearAlgebra, Tensors, OMEinsum, Rotations
using ForwardDiff, StaticArrays
using SymPy
using Symbolics

include("array_utils.jl")
include("bases.jl")
include("tens.jl")
include("tens_isotropic.jl")
include("tens_walpole.jl")
include("structured_tens_ops.jl")
include("structured_tens_promotion.jl")
include("tens_projection.jl")
include("special_tens.jl")
include("coorsystems.jl")
include("coorsystems_num.jl")
include("submanifold.jl")

end # module
