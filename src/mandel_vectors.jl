#########################################
# Symmetric tensors as Kelvin–Mandel vectors
#########################################

# Constitutive laws are written on the six Kelvin–Mandel components of a
# symmetric second-order tensor of dimension 3, a static vector, rather than on
# a tensor object: the components of the strain, the stress, and their
# operations in a few flops. These are the operations of that representation,
# in the order of TensND and Tensors.jl, (11, 22, 33, √2·23, √2·13, √2·12).

"""
    Mandel{T}

A symmetric second-order tensor of dimension 3 in **Kelvin–Mandel** components,
an `SVector{6, T}` ordered

```math
(v_{11},\\ v_{22},\\ v_{33},\\ \\sqrt2\\,v_{23},\\ \\sqrt2\\,v_{13},\\ \\sqrt2\\,v_{12}),
```

the order of [`mandel`](@ref) and of Tensors.jl's `tomandel`. The basis is
orthonormal: ``\\boldsymbol{\\sigma}:\\boldsymbol{\\varepsilon} =
\\boldsymbol{\\sigma}\\cdot\\boldsymbol{\\varepsilon}`` (`dot(σ, ε)`), the norm of
the tensor is `norm(v)`, and a fourth-order tensor with the minor symmetries is
an `SMatrix{6, 6}` acting by the matrix product.

See also [`identity_mandel`](@ref), [`trace_mandel`](@ref), [`deviator`](@ref),
[`sqnorm_deviator`](@ref), [`von_mises`](@ref), [`equivalent_strain`](@ref).
"""
const Mandel{T} = SVector{6, T}

"""
    identity_mandel(T = Float64) -> Mandel{T}

The second-order identity ``\\boldsymbol{1}``, in Kelvin–Mandel components.
"""
identity_mandel(::Type{T} = Float64) where {T} = Mandel{T}(one(T), one(T), one(T), zero(T), zero(T), zero(T))

"""
    trace_mandel(v::Mandel)

The trace ``\\operatorname{tr}\\boldsymbol{v} = \\boldsymbol{v} : \\boldsymbol{1}``.
"""
trace_mandel(v::Mandel) = v[1] + v[2] + v[3]

"""
    deviator(v::Mandel) -> Mandel

The deviatoric part ``\\boldsymbol{v} - \\tfrac13 (\\operatorname{tr}\\boldsymbol{v})\\,\\boldsymbol{1}``.
"""
deviator(v::Mandel) = v - (trace_mandel(v) / 3) * identity_mandel(eltype(v))

"""
    von_mises(σ::Mandel)

The von Mises equivalent stress ``\\sigma_{eq} = \\sqrt{\\tfrac32\\,\\boldsymbol{s}:\\boldsymbol{s}}``,
``\\boldsymbol{s}`` being the deviator of ``\\boldsymbol{\\sigma}``: equal to
``|\\sigma|`` in uniaxial tension.

!!! note "Not differentiable at zero"
    The square root has no derivative at a vanishing deviator. A law whose
    response must be differentiated there is written with ``\\sigma_{eq}^2``
    instead.
"""
von_mises(σ::Mandel) = sqrt(3 * sqnorm_deviator(σ) / 2)

"""
    equivalent_strain(ε::Mandel)

The equivalent strain ``\\varepsilon_{eq} = \\sqrt{\\tfrac23\\,\\boldsymbol{e}:\\boldsymbol{e}}``,
``\\boldsymbol{e}`` being the deviator of ``\\boldsymbol{\\varepsilon}``: work-conjugate
to [`von_mises`](@ref), and equal to ``|\\varepsilon|`` for an incompressible
uniaxial strain.
"""
equivalent_strain(ε::Mandel) = sqrt(2 * sqnorm_deviator(ε) / 3)

"""
    sqnorm_deviator(v::Mandel)

The squared norm of the deviator, ``\\boldsymbol{s}:\\boldsymbol{s}``, without a
square root: the form a law differentiates where the deviator vanishes
(``\\sigma_{eq}^2 = \\tfrac32\\,\\boldsymbol{s}:\\boldsymbol{s}``). It is formed from
the differences of the normal components, so that it is exactly zero for a
change of volume, whatever the rounding, where subtracting the mean would leave
a few ulps (and a law with a power of the deviatoric strain would take them for
a deviatoric strain, or produce 0 × ∞ in its derivatives).
"""
sqnorm_deviator(v::Mandel) =
    ((v[1] - v[2])^2 + (v[2] - v[3])^2 + (v[3] - v[1])^2) / 3 + v[4]^2 + v[5]^2 + v[6]^2

export Mandel, identity_mandel, trace_mandel, deviator, sqnorm_deviator, von_mises, equivalent_strain
