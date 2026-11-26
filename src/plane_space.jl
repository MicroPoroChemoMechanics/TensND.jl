# ─────────────────────────────────────────────────────────────────────────────
#  plane_space.jl — a tensor of space in one of its planes, a tensor of the
#  plane in space.
#
#  The plane normal to the canonical axis e_k is spanned by the two other
#  canonical axes in the cyclic order, (e_i, e_j) = (e₁, e₂) for k = 3,
#  (e₂, e₃) for k = 1, (e₃, e₁) for k = 2 (`_plane_axes` of
#  `tens_plane_classes.jl`), which are the axes 1 and 2 of the plane.
#
#  - `plane_block`: the restriction of a tensor of space to the plane, its
#    components with every index in the plane: for a stiffness, the
#    plane-strain stiffness of any anisotropy.
#  - `plane_stress_block`: the stiffness of the plane under a zero stress on
#    the facets of normal e_k, the Schur complement of the out-of-plane block
#    of the Kelvin-Mandel matrix.
#  - `space_embedding`: the tensor of space whose restriction is a given tensor
#    of the plane, uncoupled from the axis, with given entries along it.
#
#  `plane_strain_stiffness` and `plane_stress_stiffness` give the same tensors
#  for the structured classes, as structured tensors in their material frame;
#  the functions here take any tensor and read its canonical components, so
#  that the basis a tensor is written in does not change the result.
#
#  The Kelvin-Mandel slots come from `mandel_pairs` (`kelvin_mandel.jl`). The
#  storage of a result follows the type of the input, not its values: a tensor
#  whose type guarantees its (minor) symmetries gives a `SymmetricTensor`.
# ─────────────────────────────────────────────────────────────────────────────

# ── Indices and slots of the plane normal to e_k ─────────────────────────────

# The index of the plane carried by each index of space, 0 for the normal one.
function _plane_local(k::Integer)
    i, j = _plane_axes(k)
    return ntuple(m -> m == i ? 1 : (m == j ? 2 : 0), Val(3))
end

# The Kelvin-Mandel slots of space of the plane normal to e_k: `inplane[s]` is
# the slot of the slot `s` of the plane (11, 22, √2·12 of the plane), and
# `outplane` holds the slot of e_k⊗e_k, then those of the two out-of-plane
# shears. The weights are those of the plane: a diagonal pair stays diagonal.
function _plane_km_slots_of(k::Int)
    P = mandel_pairs(Val(3))
    ax = _plane_axes(k)
    slot(a, b) = findfirst(==(minmax(a, b)), P)
    inplane = map(((a, b),) -> slot(ax[a], ax[b]), mandel_pairs(Val(2)))
    normal = slot(k, k)
    shears = Tuple(s for s in eachindex(P) if s ∉ inplane && s != normal)
    return (inplane = SVector(inplane), outplane = SVector(normal, shears...))
end

# Computed once for the three planes: the kernels read them, never search.
const _PLANE_KM_SLOTS = ntuple(_plane_km_slots_of, 3)

_plane_km_slots(k::Integer) = (_plane_axes(k); @inbounds _PLANE_KM_SLOTS[k])

# ── Storage decided by the type ──────────────────────────────────────────────

# Whether the type of a tensor guarantees its minor symmetries (order 4) or its
# symmetry (order 2). Any other tensor is stored as `Tens` stores a full
# Tensors.jl tensor, on the values of its components.
_symmetric_by_type(::AbstractTens) = Val(false)
_symmetric_by_type(::TensArray{<:Any, <:Any, <:Any, <:Tensors.SymmetricTensor}) = Val(true)
_symmetric_by_type(::Union{TensISO, TensOrtho, TensCubic, TensTI{4}, TensTI{2, <:Any, 2}}) = Val(true)

_storage(::Val{true}, ::Val{order}, ::Val{dim}, ::Type{T}, f::F) where {order, dim, T, F} = SymmetricTensor{order, dim, T}(f)
_storage(::Val{false}, ::Val{order}, ::Val{dim}, ::Type{T}, f::F) where {order, dim, T, F} = Tensor{order, dim, T}(f)

# ── Restriction to a plane ───────────────────────────────────────────────────

"""
    plane_block(t::AbstractTens{4,3}, k = 3) → Tens{4,2}
    plane_block(t::AbstractTens{2,3}, k = 3) → Tens{2,2}

The restriction of a tensor of space to the plane normal to the canonical axis
``\\underline{e}_k``: its canonical components with every index in the plane,

```math
t'_{abcd} = t_{\\pi(a)\\pi(b)\\pi(c)\\pi(d)} , \\qquad t'_{ab} = t_{\\pi(a)\\pi(b)} ,
\\qquad a, b, c, d \\in \\{1, 2\\} ,
```

where ``(\\pi(1), \\pi(2))`` are the two other canonical axes in the cyclic
order, `(e₁, e₂)` for `k = 3`, `(e₂, e₃)` for `k = 1`, `(e₃, e₁)` for `k = 2`.
The tensor is first read in canonical components, so that the basis it is
written in does not change the result. The result is a general tensor on the
canonical basis of the plane, stored as a `SymmetricTensor` when the type of
`t` guarantees its (minor) symmetries.

For a stiffness it is the plane-strain stiffness in that plane, whatever the
anisotropy: the strains that carry the index `k` vanish, and the stresses of
the plane are those of the block. [`plane_strain_stiffness`](@ref) gives the
same tensor for the orthotropic, cubic, transversely isotropic and isotropic
classes, as a tensor of the plane of the same class.

See also [`plane_stress_block`](@ref), [`space_embedding`](@ref).
"""
plane_block(t::AbstractTens{4, 3}, k::Integer = 3) = _plane_block(_symmetric_by_type(t), canonical_components(t), _plane_axes(k))
plane_block(t::AbstractTens{2, 3}, k::Integer = 3) = _plane_block(_symmetric_by_type(t), canonical_components(t), _plane_axes(k))

_plane_block(s::Val, C::AbstractArray{T, 4}, ax::NTuple{2, Integer}) where {T} =
    Tens(_storage(s, Val(4), Val(2), T, (a, b, c, d) -> C[ax[a], ax[b], ax[c], ax[d]]))
_plane_block(s::Val, C::AbstractArray{T, 2}, ax::NTuple{2, Integer}) where {T} =
    Tens(_storage(s, Val(2), Val(2), T, (a, b) -> C[ax[a], ax[b]]))

"""
    plane_stress_block(ℂ::AbstractTens{4,3}, k = 3) → Tens{4,2}

The plane-stress stiffness of `ℂ` in the plane normal to the canonical axis
``\\underline{e}_k``, with the plane axes of [`plane_block`](@ref): the stress
vector on the facets of normal ``\\underline{e}_k`` vanishes,
``\\boldsymbol{\\sigma}\\cdot\\underline{e}_k = 0``. With the Kelvin-Mandel matrix
``\\boldsymbol{M}`` of the minor-symmetric part of `ℂ` split into the three slots
of the plane, ``p`` (the two normal strains of the plane and its shear), and
the three slots that carry the index ``k``, ``o`` (the normal one and the two
out-of-plane shears), the stresses ``\\sigma_o = 0`` give
``\\varepsilon_o = -\\boldsymbol{M}_{oo}^{-1}\\,\\boldsymbol{M}_{op}\\,\\varepsilon_p`` and

```math
\\boldsymbol{M}' = \\boldsymbol{M}_{pp} - \\boldsymbol{M}_{po}\\,\\boldsymbol{M}_{oo}^{-1}\\,\\boldsymbol{M}_{op} ,
```

the inverse of the block ``pp`` of the compliance. No major symmetry is
assumed. The 3×3 inverse is taken by cofactors ([`inv3`](@ref)), with no pivot,
so that duals and symbols go through.

For an orthotropic tensor whose frame is canonical the out-of-plane shears are
uncoupled, and the formula is that of [`plane_stress_stiffness`](@ref),
``C'_{ab} = C_{ab} - C_{ak}\\,C_{kb}/C_{kk}``.
"""
function plane_stress_block(ℂ::AbstractTens{4, 3}, k::Integer = 3)
    M = SMatrix{6, 6}(_mandel_minor(canonical_components(ℂ)))
    p, o = _plane_km_slots(k)
    return from_mandel(M[p, p] - M[p, o] * inv3(M[o, o]) * M[o, p])
end

# ── Embedding of a tensor of the plane ───────────────────────────────────────

"""
    space_embedding(t::AbstractTens{2,2}, k = 3; axial = 0) → Tens{2,3}
    space_embedding(t::AbstractTens{4,2}, k = 3; axial = 0, shear = 0) → Tens{4,3}

The tensor of space whose restriction to the plane normal to the canonical
axis ``\\underline{e}_k`` ([`plane_block`](@ref), with the same plane axes) is
`t`, with no coupling between the plane and the axis.

- Order 2: the components of `t` on the plane, `axial` on
  ``\\underline{e}_k\\otimes\\underline{e}_k``, zero elsewhere. The result is
  stored as a `SymmetricTensor` when the type of `t` guarantees its symmetry.
- Order 4: a tensor with the minor symmetries, whose Kelvin-Mandel matrix holds
  that of `t` on the three slots of the plane, `axial` on the diagonal slot of
  ``\\underline{e}_k\\otimes\\underline{e}_k`` (``C_{kkkk}``), `shear` on the
  diagonal slots of the two out-of-plane shears (Kelvin-Mandel entries, that is
  ``2\\,C_{akak}`` for each axis ``a`` of the plane), and zero elsewhere. The
  components of `t` read are those of its minor-symmetric part.

`plane_block(space_embedding(t, k; kw...), k) == t` whatever the keywords, for
a tensor of order 4 with the minor symmetries.
"""
function space_embedding(t::AbstractTens{2, 2}, k::Integer = 3; axial = zero(eltype(t)))
    A = canonical_components(t)
    T = promote_type(eltype(A), typeof(axial))
    return _space_embedding(_symmetric_by_type(t), A, T, _plane_local(k), convert(T, axial))
end

function space_embedding(t::AbstractTens{4, 2}, k::Integer = 3; axial = zero(eltype(t)), shear = zero(eltype(t)))
    A = canonical_components(t)
    T = halves_eltype(promote_type(eltype(A), typeof(axial), typeof(shear)))
    # The Kelvin-Mandel entry of the shear slot (a, k) is 2C_akak.
    return _space_embedding(A, T, _plane_local(k), convert(T, axial), convert(T, shear) / 2)
end

# The element type is a static parameter here, so that the closures that fill
# the components convert to a type known at compile time. An index of space is
# in the plane when it carries an index of the plane (`loc[i] != 0`).
function _space_embedding(s::Val, A::AbstractMatrix, ::Type{T}, loc::NTuple{3, Int}, axial::T) where {T}
    entry(i, j) = loc[i] != 0 && loc[j] != 0 ? convert(T, A[loc[i], loc[j]]) : (i == j ? axial : zero(T))
    return Tens(_storage(s, Val(2), Val(3), T, entry))
end

# Read on the independent entries: (i, j) and (l, m) are each the pair of a
# slot. A shear entry pairs the normal index with the same axis of the plane on
# both sides.
function _space_embedding(A::AbstractArray{<:Any, 4}, ::Type{T}, loc::NTuple{3, Int}, axial::T, h::T) where {T}
    function entry(i, j, l, m)
        a, b, c, d = loc[i], loc[j], loc[l], loc[m]
        a != 0 && b != 0 && c != 0 && d != 0 && return convert(T, _minor_part(A, a, b, c, d))
        a == 0 && b == 0 && c == 0 && d == 0 && return axial
        ((a == 0) != (b == 0) && (c == 0) != (d == 0) && a + b == c + d) && return h
        return zero(T)
    end
    return Tens(SymmetricTensor{4, 3, T}(entry))
end

export plane_block, plane_stress_block, space_embedding
