# =============================================================================
#  sympy_matrix.jl — how a Julia matrix of `Sym` becomes a `sympy.Matrix`.
#
#  Every SymPy operation on a `Matrix{Sym}` or a `Vector{Sym}` (`transpose`,
#  `adjoint`, `inv`, `det`, `norm`, `eigen`, `\`, `M.method()`, LaTeX display)
#  hands the array to Python as one `sympy.Matrix`, built by SymPyCore's `↓`.
#  The Julia code that `↓` compiles to crashes the compiler of Julia 1.15-DEV,
#  and of 1.14 once the change is backported (JuliaLang/julia#63672, a
#  segmentation fault in `LowerGCAllocMemory`): so does any Julia array of
#  `PyObject`s handed to Python in one piece, PyCall's own conversion included.
#
#  The same `sympy.Matrix` is built here row by row, each entry appended to a
#  Python list, so that no Julia array of `PyObject`s exists at any point. The
#  object reaching SymPy is the one it received before, and every operation
#  above returns exactly what it returned: only the way the rows are assembled
#  changes. The methods are more specific than SymPyCore's (an exact element
#  type), so they take precedence without overwriting them. Both element types
#  occur: `Sym{PyObject}` (`[x y; 0 x]`), and `Sym` itself, which is not concrete
#  (`Sym[1 0; 0 1]`).
# =============================================================================

function _sympy_matrix(M::AbstractMatrix)
    PyCall = SymPy.PyCall
    list = PyCall.pybuiltin("list")
    rows = PyCall.pycall(list, PyCall.PyObject)
    for i in axes(M, 1)
        row = PyCall.pycall(list, PyCall.PyObject)
        for j in axes(M, 2)
            PyCall.pycall(row."append", PyCall.PyObject, SymPy.SymPyCore.:↓(M[i, j]))
        end
        PyCall.pycall(rows."append", PyCall.PyObject, row)
    end
    return PyCall.pycall(SymPy.sympy.o."Matrix", PyCall.PyObject, rows)
end

for E in (SymPy.Sym, SymPy.Sym{SymPy.PyCall.PyObject})
    @eval begin
        SymPy.SymPyCore.:↓(M::AbstractMatrix{$E}) = _sympy_matrix(M)
        # A vector is a column, as SymPyCore makes it. `Vector` as well as
        # `AbstractVector`: SymPyCore's method is on `Vector{<:Sym}`.
        SymPy.SymPyCore.:↓(v::AbstractVector{$E}) = _sympy_matrix(reshape(v, length(v), 1))
        SymPy.SymPyCore.:↓(v::Vector{$E}) = _sympy_matrix(reshape(v, length(v), 1))
    end
end
