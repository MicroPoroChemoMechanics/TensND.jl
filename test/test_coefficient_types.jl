import Symbolics

@testsection "Element types of the structured tensors" begin

    # The components of the structured tensors carry 1/2, 1/3 (the isotropic and
    # cubic projectors, 𝕎₂, 𝕎₅) or 1/√2 (𝕎₃, 𝕎₄). Integer coefficients were
    # stored as they were, which left tensors whose components could not be
    # formed (`InexactError` on display, `get_array` or a product), and rational
    # ones rounded 1/√2 into a rational silently. Each tensor is held to the same
    # tensor built from floating-point coefficients.
    can = CanonicalBasis{3, Float64}()

    @testsection "Walpole basis on an integer axis" begin
        for (W, Wf) in zip(walpole_basis((0, 0, 1)), walpole_basis((0.0, 0.0, 1.0)))
            @test eltype(W) === Float64
            @test get_array(W) == get_array(Wf)
            @test !isempty(sprint(show, MIME"text/plain"(), W))
        end
        @test get_array(walpole_basis((0, 0, 1))[3])[3, 3, 1, 1] == 1 / sqrt(2.0)
    end

    @testsection "transversely isotropic, order 4" begin
        for args in ((1, 2, 0, 1, 1), (1, 2, 1, 1, 1, 1), (1, 2, 0, 0, 1, 1, 1, 0))
            t = TensTI{4}(args..., (0, 0, 1))
            tf = TensTI{4}(float.(args)..., (0.0, 0.0, 1.0))
            @test eltype(t) === Float64
            @test get_array(t) == get_array(tf)
            @test get_array(inv(t)) ≈ get_array(inv(tf))
            @test get_array(t ⊡ t) == get_array(tf ⊡ tf)
        end
        # 1/√2 is not a rational: the coefficients become floats, not a rational
        # rounding of the float.
        t = TensTI{4}(1 // 1, 2 // 1, 1 // 1, 1 // 1, 1 // 1, (0 // 1, 0 // 1, 1 // 1))
        @test eltype(t) === Float64
        @test get_array(t) == get_array(TensTI{4}(1.0, 2.0, 1.0, 1.0, 1.0, (0.0, 0.0, 1.0)))
        @test eltype(TensTI{4}(big(1), big(2), big(1), big(1), big(1), (big(0), big(0), big(1)))) === BigFloat
        @test eltype(@inferred TensTI{4}(1, 2, 0, 1, 1, (0, 0, 1))) === Float64
        @test (@allocated TensTI{4}(1.0, 2.0, 0.0, 1.0, 1.0, (0.0, 0.0, 1.0))) == 0
        @test_throws ArgumentError TensTI{4, Int, 5}((1, 2, 0, 1, 1), (0, 0, 1))
        @test_throws ArgumentError TensTI{4, Rational{Int}, 5}((1 // 1, 2 // 1, 0 // 1, 1 // 1, 1 // 1), (0 // 1, 0 // 1, 1 // 1))
    end

    @testsection "transversely isotropic, order 2: integer components, float inverse" begin
        t = TensTI{2}(2, 4, (0, 0, 1))
        @test eltype(t) === Int
        @test get_array(t) == [2 0 0; 0 2 0; 0 0 4]
        @test inv(t) isa TensTI{2, Float64, 2}
        @test get_array(inv(t)) == [0.5 0 0; 0 0.5 0; 0 0 0.25]
        t3 = TensTI{2}(1, 2, 1, (0, 0, 1))
        @test eltype(inv(t3)) === Float64
        @test get_array(inv(t3)) ≈ inv(Float64.(get_array(t3)))
    end

    @testsection "isotropic, order 4" begin
        for dim in (2, 3)
            t = TensISO{dim}(3, 2)
            @test eltype(t) === Float64
            @test get_array(t) == get_array(TensISO{dim}(3.0, 2.0))
        end
        # 1/3 is a rational: rational coefficients stay exact.
        tq = TensISO{3}(3 // 1, 2 // 1)
        @test eltype(tq) === Rational{Int}
        @test get_array(tq)[1, 1, 1, 1] == 3 // 3 + 2 * 2 // 3
        @test eltype(TensISO{3}(big(3), big(2))) === BigFloat
        @test eltype(TensISO{3}(2)) === Int                 # order 2: integer components
        @test get_array(TensISO{4, 3, Float64}()) == get_array(TensISO{3}(1.0, 1.0))
        @test_throws ArgumentError TensISO{4, 3, Int}()
    end

    @testsection "cubic" begin
        frame_int = CanonicalBasis{3, Int}()
        t = TensCubic(5, 2, 1, frame_int)
        @test eltype(t) === Float64
        @test get_array(t) == get_array(TensCubic(5.0, 2.0, 1.0, can))
        @test eltype(TensCubic(5 // 1, 2 // 1, 1 // 1, CanonicalBasis{3, Rational{Int}}())) === Rational{Int}
        @test_throws ArgumentError TensCubic{Int, typeof(frame_int)}((5, 2, 1), frame_int)
    end

    @testsection "symbolic, dual and other floating types are kept" begin
        x = SymPyPythonCall.symbols("x", real = true)
        Symbolics.@variables y
        d = ForwardDiff.Dual(1.0, 1.0)
        @test eltype(TensTI{4}(x, x, x, x, x, (Sym(0), Sym(0), Sym(1)))) <: Sym
        @test eltype(TensISO{3}(x, 2x)) <: Sym
        @test eltype(TensTI{4}(y, y, y, y, y, (0, 0, 1))) === Symbolics.Num
        @test eltype(TensISO{3}(y, y)) === Symbolics.Num
        @test eltype(TensTI{4}(d, d, d, d, d, (0.0, 0.0, 1.0))) === typeof(d)
        @test eltype(TensTI{4}(1.0f0, 2.0f0, 0.0f0, 1.0f0, 1.0f0, (0.0f0, 0.0f0, 1.0f0))) === Float32
        @test eltype(TensISO{3}(1 + 2im, 3)) === ComplexF64
    end
end
