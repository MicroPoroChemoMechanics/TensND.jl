# [Special tensors and rotations](@id api-special)

Levi-Civita symbol, rotation tensors, predefined unit vectors and random
tensors. Theory: [Rotations and Euler angles](@ref th-rotations).

```@docs
levi_civita
rotation_matrix
rotation_tensor4
𝐞
𝐞ᵖ
𝐞ᶜ
𝐞ˢ
```

## Random tensors and orientations

Manual: [Random tensors](@ref man-random-tensors). Theory:
[Uniform random rotations](@ref th-random-rotations).

```@docs
Base.rand(::Type{<:TensND._RandomClass})
random_rotation
random_axis
random_basis
```
