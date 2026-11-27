# =============================================================================
#  reference_contractions.jl — contractions written index by index, the
#  independent references of the tests of the algebra: they share no code with
#  the kernels of TensND. Included once by `runtests.jl`.
# =============================================================================

"Components of A in the frame whose vectors are the columns of R: A′ᵢⱼ… = Rᵢₐ Rⱼᵦ … Aₐᵦ…."
ref_rotate(R, A::AbstractMatrix) = [sum(R[i, a] * R[j, b] * A[a, b] for a in axes(A, 1), b in axes(A, 2)) for i in axes(R, 1), j in axes(R, 1)]
ref_rotate(R, A::AbstractArray{<:Any, 4}) = [
    sum(R[i, a] * R[j, b] * R[k, c] * R[l, d] * A[a, b, c, d] for a in axes(A, 1), b in axes(A, 2), c in axes(A, 3), d in axes(A, 4))
        for i in axes(R, 1), j in axes(R, 1), k in axes(R, 1), l in axes(R, 1)
]

"The converse: A′ₐᵦ… = Rᵢₐ Rⱼᵦ … Aᵢⱼ…."
ref_rotate_back(R, A) = ref_rotate(transpose(R), A)

"The tensor product (a ⊗ b)ᵢⱼₖₗ = aᵢⱼ bₖₗ."
ref_outer(a::AbstractMatrix, b::AbstractMatrix) = [a[i, j] * b[k, l] for i in axes(a, 1), j in axes(a, 2), k in axes(b, 1), l in axes(b, 2)]

"Double contraction on the two inner indices: (L : a)ᵢⱼ, (a : L)ₖₗ and (A : B)ᵢⱼₘₙ."
ref_dc(L::AbstractArray{<:Any, 4}, a::AbstractMatrix) = [sum(L[i, j, k, l] * a[k, l] for k in axes(L, 3), l in axes(L, 4)) for i in axes(L, 1), j in axes(L, 2)]
ref_dc(a::AbstractMatrix, L::AbstractArray{<:Any, 4}) = [sum(a[i, j] * L[i, j, k, l] for i in axes(L, 1), j in axes(L, 2)) for k in axes(L, 3), l in axes(L, 4)]
ref_dc(A::AbstractArray{<:Any, 4}, B::AbstractArray{<:Any, 4}) = [
    sum(A[i, j, k, l] * B[k, l, m, n] for k in axes(A, 3), l in axes(A, 4))
        for i in axes(A, 1), j in axes(A, 2), m in axes(B, 3), n in axes(B, 4)
]

"Full contraction Aᵢⱼₖₗ Bᵢⱼₖₗ, without conjugation."
ref_full(A, B) = sum(A[I] * B[I] for I in eachindex(A, B))

"The acoustic contraction aₖ Aₖᵢⱼₗ bₗ."
ref_acoustic(a, A, b) = [sum(a[k] * A[k, i, j, l] * b[l] for k in axes(A, 1), l in axes(A, 4)) for i in axes(A, 2), j in axes(A, 3)]
