# The `sympy.Matrix` that every SymPy operation on a `Matrix{Sym}` builds is
# assembled row by row in the SymPyPythonCall extension of TensND, so that no Julia array of
# Python objects exists (the one-piece conversion crashes the compiler of Julia
# 1.15-DEV, JuliaLang/julia#63672). Each operation that goes through it is
# checked against its closed form, so a wrong assembly (a transposed matrix,
# a row vector instead of a column) fails here.

@testsection "SymPy matrices assembled row by row" begin
    @syms a b c d e
    M = [a b; c d]
    # The object handed to SymPy: its shape, for a matrix and for a vector (a
    # Python tuple, converted to compare).
    shape(x) = SymPyPythonCall.PythonCall.pyconvert(Tuple{Int, Int}, SymPyPythonCall.SymPyCore.:↓(x).shape)
    @test shape(M) == (2, 2)
    @test shape([a, b, c]) == (3, 1)
    # A vector that is not a `Vector` (a view) takes the `AbstractVector`
    # method, a `Vector` its own: both are columns.
    @test shape(view([a, b, c, d], 2:4)) == (3, 1)
    @test shape(view(Sym[1, 2, 3, 4], 1:2)) == (2, 1)
    @test shape(Sym[1, 2]) == (2, 1)
    @test shape([a b c]) == (1, 3)
    # `Sym` itself is not a concrete type: `Sym[...]` arrays take the same route.
    ext = Base.get_extension(TensND, :TensNDSymPyPythonCallExt)
    @test ext isa Module
    for A in (Sym[1 0 0; 0 1 0; 0 1 1], Sym[1, 2, 3], view(Sym[1 0; 0 1], :, 1:2))
        @test parentmodule(which(SymPyPythonCall.SymPyCore.:↓, Tuple{typeof(A)})) === ext
    end
    @test parentmodule(which(SymPyPythonCall.SymPyCore.:↓, Tuple{typeof(M)})) === ext
    @test shape(Sym[1, 2, 3]) == (3, 1)
    @test transpose(Sym[1 0 0; 0 1 0; 0 1 1]) == Sym[1 0 0; 0 1 1; 0 0 1]
    @test inv(Sym[1 0 0; 0 1 0; 0 1 1]) == Sym[1 0 0; 0 1 0; 0 -1 1]

    # The operations that go through it, against their closed forms.
    @test transpose(M) == [a c; b d]
    @test adjoint(M) == conj.([a c; b d])
    @test simplify(det(M) - (a * d - b * c)) == 0
    @test all(simplify.(inv(M) .- [d -b; -c a] ./ (a * d - b * c)) .== 0)
    @test norm([a, b]) == sqrt(abs(a)^2 + abs(b)^2)
    @test M.rank() == 2
    # A view goes through the same assembly.
    N = [a b e; c d e; e e e]
    @test simplify(det(view(N, 1:2, 1:2)) - (a * d - b * c)) == 0
    @test transpose(view(N, 1:2, 2:3)) == [b d; e e]
end
