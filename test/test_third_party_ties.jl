# Calls that matched a TensND method and a method of a package TensND loads
# equally well, and were refused as ambiguous: Symbolics divides any real array
# by a `Num`, LinearAlgebra contracts a `Diagonal` with any matrix, and
# `normalize` takes a number as a second argument. Aqua now checks that no such
# tie is left; this file checks what the calls return.

import Symbolics

@testsection "Ties with Symbolics and LinearAlgebra" begin
    Symbolics.@variables x
    # The structured tensors are compared on their own coefficients, which are
    # what the division acts on (the components of a `TensTI` carry a symbolic
    # √2 that substitution does not fold).
    coefficients(t::Union{TensISO, TensTI, TensOrtho, TensCubic}) = collect(get_data(t))
    coefficients(t) = collect(get_array(t))
    at(t) = Float64.(Symbolics.value.(Symbolics.substitute.(coefficients(t), Ref(Dict(x => 2.0)))))
    n3 = 𝐞(Val(3), Val(3), Val(Float64))
    can = CanonicalBasis{3, Float64}()

    @testsection "a tensor divided by a Num keeps its structure" begin
        for A in (
                TensISO{3}(3.0, 2.0), TensISO{3}(2.0),
                TensTI{4}(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, n3),
                TensOrtho(10.0, 8.0, 9.0, 3.0, 2.0, 4.0, 2.5, 3.0, 1.5, can),
                cubic_from_constants(10.0, 4.0, 2.0, can),
                Tens([1.0 2.0 0.0; 2.0 3.0 1.0; 0.0 1.0 4.0]),
            )
            B = A / x
            @test nameof(typeof(B)) == nameof(typeof(A)) || A isa Tens
            @test eltype(B) == Symbolics.Num
            @test at(B) ≈ coefficients(A) / 2
        end
        # Entries that are themselves symbolic, the call that MeanFieldHomogenization
        # makes when it closes an average by a symbolic volume fraction.
        Symbolics.@variables k μ
        A = TensISO{3}(3k, 2μ) / x
        @test A isa TensISO{4, 3, Symbolics.Num}
        @test all(isequal.(Symbolics.simplify.(collect(get_data(A)) .* x .- [3k, 2μ]), 0))
    end

    @testsection "an isotropic tensor contracted with a Diagonal" begin
        D = Diagonal([1.0, 2.0, 3.0])
        λ = TensISO{3}(2.0)
        @test dot(λ, D) == 2.0 * D
        @test dot(D, λ) == D * 2.0
        @test λ ⋅ D == 2.0 * D
    end

    @testsection "normalize takes a variance, LinearAlgebra a p-norm" begin
        ℬ = Basis(Float64[2 0 0; 0 3 0; 0 0 4])
        @test vecbasis(normalize(ℬ), :cov) ≈ Matrix(1.0I, 3, 3)
        @test vecbasis(normalize(ℬ, :cov), :cov) ≈ Matrix(1.0I, 3, 3)
        @test which(normalize, Tuple{typeof(ℬ), Float64}).module === LinearAlgebra
    end
end
