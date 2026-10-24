# The paths that are type-stable by construction must stay so: a regression to
# dynamic dispatch costs no test its result, only time, and nothing else would
# notice it. Each case is inferred, and those that build only isbits results
# allocate nothing.
@testsection "Hot paths: inference and allocations" begin
    _allocs(f, args...) = (f(args...); @allocated f(args...))

    ISO = TensISO{3}(3.0, 2.0)
    TI = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, [0.0, 0.0, 1.0])
    OR = TensOrtho(250.0, 180.0, 150.0, 90.0, 80.0, 70.0, 60.0, 50.0, 40.0, CanonicalBasis{3, Float64}())
    # Structured algebra: closed forms on the coefficients.
    for (f, args) in ((inv, (ISO,)), (⊡, (ISO, ISO)), (inv, (TI,)), (inv, (OR,)))
        @inferred f(args...)
        @test _allocs(f, args...) == 0
    end
    # Two TI tensors are compared for a common axis at run time, with the
    # generic product as the fallback, so the product is not inferred; the
    # coaxial case still allocates nothing.
    @test _allocs(⊡, TI, TI) == 0
    # An isotropic operand is lifted onto the axis of the TI one, so the product,
    # the sum and the difference skip that comparison and are inferred; so is the
    # chain of a dilute or Mori-Tanaka concentration, `C₁ ⊡ (𝕀 + ℙ ⊡ δℂ)⁻¹`.
    TI6, TI8 = TI ⊡ ISO, TensND._lift_walpole_N8(TI)
    @test _allocs(⊡, TI8, TI8) == 0 && TI8 ⊡ TI8 ≈ generic_tens(TI8) ⊡ generic_tens(TI8)
    for t in (TI, TI6, TI8), (f, args) in ((⊡, (ISO, t)), (⊡, (t, ISO)), (+, (ISO, t)), (-, (t, ISO)))
        @inferred f(args...)
        @test _allocs(f, args...) == 0
    end
    concentration(P, δC, C₁, 𝕀) = C₁ ⊡ inv(𝕀 + P ⊡ δC)
    @inferred concentration(TI, ISO, ISO, TensISO{3}(1.0, 1.0))
    # Kelvin-Mandel matrices of structured tensors, at the dimension of the type.
    for t in (ISO, TI, OR, tens_cubic(250.0, 120.0, 90.0, CanonicalBasis{3, Float64}()))
        @test (@inferred KM(t)) isa Matrix{Float64}
    end
    # The product of a dual tensor with a plain one, both canonical: no change
    # of basis is needed, and none is dispatched at run time.
    TId = tens_TI(ForwardDiff.Dual{Nothing}(10.0, 1.0), 3.0, 2.5, 12.0, 2.0, [0.6, 0.0, 0.8])
    Gd, Ga = TensND.generic_tens(TId), TensND.generic_tens(TI)
    @inferred Gd ⊡ Ga
    @test _allocs(⊡, Gd, Ga) == 0
    # A static Kelvin-Mandel array carries its size in its type: the tensor
    # type comes by dispatch, and the result is the one of the plain array.
    # Inferred when the size gives a symmetric tensor; a full one is stored
    # symmetric or not according to its components, a value.
    SA = TensND.StaticArrays
    static(v) = v isa AbstractMatrix ? SA.SMatrix{size(v)...}(v) : SA.SVector{length(v)}(v)
    ℬ3, ℬ2 = RotatedBasis(0.3, 0.2, 0.1), RotatedBasis(0.3)
    M6 = [10.0 3 1 0 0 0; 3 8 0.5 0 0 0; 1 0.5 4 0 0 0; 0 0 0 2 0 0; 0 0 0 0 3 0; 0 0 0 0 0 1.5]
    for (v, ℬ) in ((M6, ℬ3), (M6[1:3, 1:3], ℬ2), (M6[:, 1], ℬ3), (M6[1:3, 1], ℬ2))
        @test get_array(@inferred inv_KM(static(v))) == get_array(inv_KM(v))
        @test get_array(@inferred inv_KM(static(v), ℬ)) == get_array(inv_KM(v, ℬ))
    end
    M9 = [0.1 * (i + 2j) + (i == j) for i in 1:9, j in 1:9]
    M4 = [1.0 2 0 0; 2 5 0 0; 0 0 3 1; 0 0 1 4]
    for (v, ℬ) in ((M9, ℬ3), (M9[:, 1], ℬ3), (M4, ℬ2), (M4[:, 1], ℬ2))
        @test get_array(inv_KM(static(v))) == get_array(inv_KM(v))
        @test get_array(inv_KM(static(v), ℬ)) == get_array(inv_KM(v, ℬ))
    end
    # The parent of an orthogonal basis is in its type.
    ℬo = Basis(RotatedBasis(0.3, 0.2, 0.1), [2.0, 1.0, 0.5])
    @test (@inferred TensND.relevant_OrthonormalBasis(ℬo)) isa RotatedBasis{3, Float64}
    # The contraction kernels.
    A4 = [0.1 * (i + 2j) * (k + l) + (i == k) for i in 1:3, j in 1:3, k in 1:3, l in 1:3]
    A2 = [1.0 2.0 0.5; 0.3 4.0 1.0; 0.2 0.1 3.0]
    @inferred TensND._contract_index(A4, A2, 2)
    @inferred Tensors.dcontract(A4, A2)

    # Change of basis, Kelvin-Mandel, products and inverses on the four
    # structured classes and on generic tensors, canonical or rotated.
    ℬr = RotatedBasis(0.7, -0.4, 1.1)
    CU = tens_cubic(60.0, 16.0, 10.0, ℬr)
    ORr = TensOrtho(160.0, 120.0, 95.0, 55.0, 48.0, 42.0, 33.0, 29.0, 27.0, ℬr)
    G = Tens(SymmetricTensor{4, 3}((i, j, k, l) -> Float64(i + j + k + l)))
    S2 = Tens(SymmetricTensor{2, 3}((i, j) -> Float64(i * j)))
    GR, S2R = change_tens(G, ℬr), change_tens(S2, ℬr)
    for t in (ISO, TI, ORr, CU)
        @test (@inferred change_tens(t, ℬr)) isa TensND.TensRotated{4, 3, Float64, SymmetricTensor{4, 3, Float64, 36}}
        @inferred KM(t, ℬr)
        @inferred inv(t)
        # (TI + TI checks for a common axis at run time, like TI ⊡ TI above.)
        t isa TensTI || @inferred t + t
    end
    # An isotropic tensor in a rotated basis: the same components as before,
    # now at the static dimension of its type.
    @test get_array(change_tens(ISO, ℬr)) == get_array(Tens(get_array(ISO), ℬr))
    @inferred change_tens(TensISO{3}(2.0), ℬr)   # order 2: the constructor takes the dimension
    for (f, args) in (
            (KM, (G,)), (KM, (GR,)), (KM, (GR, CanonicalBasis{3, Float64}())),
            (⊗, (S2, S2)), (⊗, (S2R, S2R)), (⊡, (G, S2)), (⊡, (GR, S2R)),
            (inv, (G,)), (inv, (GR,)), (+, (GR, GR)), (isotropify, (ORr,)), (KM_rotation, (Matrix(ℬr),)),
        )
        @inferred f(args...)
    end
    # Small unions decided by a value, by design: a tensor already in the
    # target basis is returned as it is, a basis that is the identity is the
    # canonical one, and the storage is symmetric only if the components are.
    @inferred Union{typeof(G), typeof(GR)} change_tens(G, ℬr)
    # Two orthotropic tensors on the same frame multiply in the twelve-constant
    # form, on different frames by the generic route.
    @inferred Union{TensOrtho{4, Float64, 12, typeof(ℬr)}, TensND.TensCanonical{4, 3, Float64, SymmetricTensor{4, 3, Float64, 36}}} ORr ⊡ ORr
    # The product is inlined: a caller splits the union and allocates nothing
    # (only the union itself, returned from a non-inlined call, is boxed).
    @test _allocs((a, b) -> get_data(a ⊡ b), ORr, ORr) == 0
    @inferred Union{typeof(G), typeof(GR)} change_tens(GR, CanonicalBasis{3, Float64}())
    @inferred Union{Array{Float64, 4}, SymmetricTensor{4, 3, Float64, 36}} components_canon(GR)
    @inferred Union{RotatedBasis{<:Any, Float64}, CanonicalBasis{<:Any, Float64}} frame_from_axis((0.3, -0.5, 0.8))
    # An isotropic tensor combined with a generic one: the result is at the
    # static dimension of the types, its storage decided by its components.
    GRs = Union{TensND.TensRotated{4, 3, Float64, SymmetricTensor{4, 3, Float64, 36}, typeof(ℬr)}, TensND.TensRotated{4, 3, Float64, Tensor{4, 3, Float64, 81}, typeof(ℬr)}}
    for (f, args) in ((+, (ISO, GR)), (-, (GR, ISO)), (⊡, (ISO, GR)), (⊡, (GR, ISO)))
        @inferred GRs f(args...)
    end
end
