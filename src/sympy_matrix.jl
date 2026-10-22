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
#  changes. The methods are more specific than SymPyCore's (a concrete element
#  type), so they take precedence without overwriting them.
# =============================================================================

const _SymPyObject = SymPy.Sym{SymPy.PyCall.PyObject}

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

SymPy.SymPyCore.:↓(M::AbstractMatrix{_SymPyObject}) = _sympy_matrix(M)
# A vector is a column, as SymPyCore makes it.
SymPy.SymPyCore.:↓(v::AbstractVector{_SymPyObject}) = _sympy_matrix(reshape(v, length(v), 1))
