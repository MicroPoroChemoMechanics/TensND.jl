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
    # Kelvin-Mandel matrices of structured tensors, at the dimension of the type.
    for t in (ISO, TI, OR, tens_cubic(250.0, 120.0, 90.0, CanonicalBasis{3, Float64}()))
        @test (@inferred KM(t)) isa Matrix{Float64}
    end
    # The product of a dual tensor with a plain one, both canonical: no change
    # of basis is needed, and none is dispatched at run time.
    TId = tens_TI(ForwardDiff.Dual{Nothing}(10.0, 1.0), 3.0, 2.5, 12.0, 2.0, [0.6, 0.0, 0.8])
    Gd, Ga = TensND._generic_tens(TId), TensND._generic_tens(TI)
    @inferred Gd ⊡ Ga
    @test _allocs(⊡, Gd, Ga) == 0
    # A static Kelvin-Mandel array carries its size in its type: the tensor
    # type comes by dispatch, and the result is the one of the plain array.
    # Inferred when the size gives a symmetric tensor; a full one is stored
    # symmetric or not according to its components, a value.
    SA = TensND.StaticArrays
    static(v) = v isa AbstractMatrix ? SA.SMatrix{size(v)...}(v) : SA.SVector{length(v)}(v)
    M6 = [10.0 3 1 0 0 0; 3 8 0.5 0 0 0; 1 0.5 4 0 0 0; 0 0 0 2 0 0; 0 0 0 0 3 0; 0 0 0 0 0 1.5]
    for v in (M6, M6[1:3, 1:3], M6[:, 1], M6[1:3, 1])
        @test get_array(@inferred inv_KM(static(v))) == get_array(inv_KM(v))
        @test get_array(@inferred inv_KM(static(v), RotatedBasis(0.3))) == get_array(inv_KM(v, RotatedBasis(0.3)))
    end
    for v in ([1.0 2 0 0; 2 5 0 0; 0 0 3 1; 0 0 1 4], [1.0, 2, 3, 4])
        @test get_array(inv_KM(static(v))) == get_array(inv_KM(v))
        @test get_array(inv_KM(static(v), RotatedBasis(0.3))) == get_array(inv_KM(v, RotatedBasis(0.3)))
    end
    # The parent of an orthogonal basis is in its type.
    ℬo = Basis(RotatedBasis(0.3, 0.2, 0.1), [2.0, 1.0, 0.5])
    @test (@inferred TensND.relevant_OrthonormalBasis(ℬo)) isa RotatedBasis{3, Float64}
    # The contraction kernels.
    A4 = [0.1 * (i + 2j) * (k + l) + (i == k) for i in 1:3, j in 1:3, k in 1:3, l in 1:3]
    A2 = [1.0 2.0 0.5; 0.3 4.0 1.0; 0.2 0.1 3.0]
    @inferred TensND._contract_index(A4, A2, 2)
    @inferred Tensors.dcontract(A4, A2)
end
