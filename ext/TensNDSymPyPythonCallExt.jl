# =============================================================================
#  TensNDSymPyPythonCallExt — what TensND does with SymPy that calls Python.
#
#  TensND dispatches on the SymPy type `Sym` through SymPyCore, which is pure
#  Julia and light; everything that reaches Python lives here, loaded with
#  SymPyPythonCall, so that `using TensND` no longer loads PythonCall (2 s of
#  the 5 s that loading MeanFieldHomogenization took). Three pieces: the
#  elementwise symbolic operations on a `Sym` (`tsimplify`, `tlimit`, …), the
#  default coordinates of the coordinate systems (SymPy symbols), and the
#  conversion of a Julia matrix of `Sym` to a `sympy.Matrix`, below.
# =============================================================================

module TensNDSymPyPythonCallExt

using TensND
using SymPyPythonCall
using SymPyPythonCall: Sym, sympy, symbols
import TensND: tsimplify, tfactor, tsubs, tdiff, ttrigsimp, texpand_trig, tlimit
import TensND: coorsys_cartesian, coorsys_polar, coorsys_cylindrical, coorsys_spherical, coorsys_spheroidal
import TensND: init_cartesian, init_polar, init_cylindrical, init_spherical, init_rotated

# ── Elementwise symbolic operations on a `Sym` ──────────────────────────────

tlimit(x::T, args...; kwargs...) where {T <: Sym} = sympy.limit(x, args...; kwargs...)
for OP in (:(simplify), :(factor), :(subs), :(diff))
    @eval $(Symbol("t", OP))(x::T, args...; kwargs...) where {T <: Sym} = SymPyPythonCall.$OP(x, args...; kwargs...)
end
for OP in (:(trigsimp), :(expand_trig))
    @eval $(Symbol("t", OP))(x::T, args...; kwargs...) where {T <: Sym} = sympy.$OP(x, args...; kwargs...)
end

# ── Default coordinates: SymPy symbols ──────────────────────────────────────
#
# Each set once, read by every method that defaults to it; the coordinates
# keep the assumptions they always had (a radius positive, angles real).

_cartesian_symbols(::Val{3}) = symbols("x y z", real = true)
_cartesian_symbols(::Val{2}) = symbols("x y", real = true)
_polar_system_symbols() = (symbols("r", positive = true), symbols("θ", real = true))
_cylindrical_symbols() = (symbols("r", positive = true), symbols("θ", real = true), symbols("z", real = true))
_spherical_symbols() = (symbols("θ", real = true), symbols("ϕ", real = true), symbols("r", positive = true))
_spheroidal_symbols() = (symbols("ϕ", real = true), symbols("p", real = true), symbols("q", positive = true))
_spheroidal_focal() = symbols("c", positive = true)
_spheroidal_auxiliary() = symbols("p̄ q̄", positive = true)

coorsys_cartesian() = coorsys_cartesian(_cartesian_symbols(Val(3)))
coorsys_cartesian(v::Val{<:Sym}) = coorsys_cartesian(v, _cartesian_symbols(Val(3)))
coorsys_polar(; canonical = false) = coorsys_polar(_polar_system_symbols(); canonical)
coorsys_polar(v::Val{<:Sym}; canonical = false) = coorsys_polar(v, _polar_system_symbols(); canonical)
coorsys_cylindrical(; canonical = false) = coorsys_cylindrical(_cylindrical_symbols(); canonical)
coorsys_cylindrical(v::Val{<:Sym}; canonical = false) = coorsys_cylindrical(v, _cylindrical_symbols(); canonical)
coorsys_spherical(; canonical = false) = coorsys_spherical(_spherical_symbols(); canonical)
coorsys_spherical(v::Val{<:Sym}; canonical = false) = coorsys_spherical(v, _spherical_symbols(); canonical)
coorsys_spheroidal(coords::NTuple{3} = _spheroidal_symbols(), c = _spheroidal_focal()) =
    coorsys_spheroidal(coords, c, _spheroidal_auxiliary())
coorsys_spheroidal(v::Val{<:Sym}, coords = _spheroidal_symbols(), c = _spheroidal_focal()) =
    coorsys_spheroidal(v, coords, c, _spheroidal_auxiliary())

init_cartesian() = init_cartesian(_cartesian_symbols(Val(3)))
init_cartesian(v::Union{Val{2}, Val{3}}) = init_cartesian(_cartesian_symbols(v))
init_polar(; canonical = false) = init_polar(symbols("r θ", real = true); canonical)
init_cylindrical(; canonical = false) = init_cylindrical((symbols("r", positive = true), symbols("θ z", real = true)...); canonical)
init_spherical(; canonical = false) = init_spherical((symbols("θ ϕ", real = true)..., symbols("r", positive = true)); canonical)
init_rotated(; canonical = false) = init_rotated(symbols("θ ϕ ψ", real = true); canonical)

# ── A Julia matrix of `Sym` handed to SymPy ─────────────────────────────────
#
# Every SymPy operation on a `Matrix{Sym}` or a `Vector{Sym}` (`transpose`,
# `adjoint`, `inv`, `det`, `norm`, `eigen`, `\`, `M.method()`, LaTeX display)
# hands the array to Python as one `sympy.Matrix`, built by SymPyCore's `↓`.
# The Julia code that `↓` compiles to crashes the compiler of Julia 1.15-DEV,
# and of 1.14 once the change is backported (JuliaLang/julia#63672, a
# segmentation fault in `LowerGCAllocMemory`): so does any Julia array of
# Python objects (`PythonCall.Py`) handed to Python in one piece.
# The same `sympy.Matrix` is built here row by row, each entry appended to a
# Python list, so that no Julia array of Python objects exists at any point. The
# object reaching SymPy is the one it received before, and every operation
# above returns exactly what it returned: only the way the rows are assembled
# changes. The methods are more specific than SymPyCore's (an exact element
# type), so they take precedence without overwriting them. Both element types
# occur: `Sym{Py}` (`[x y; 0 x]`), and `Sym` itself, which is not concrete
# (`Sym[1 0; 0 1]`).

function _sympy_matrix(M::AbstractMatrix)
    PythonCall = SymPyPythonCall.PythonCall
    rows = PythonCall.pylist()
    for i in axes(M, 1)
        row = PythonCall.pylist()
        for j in axes(M, 2)
            row.append(SymPyPythonCall.SymPyCore.:↓(M[i, j]))
        end
        rows.append(row)
    end
    return SymPyPythonCall.sympy.o.Matrix(rows)
end

_sympy_column(v::AbstractVector) = _sympy_matrix(reshape(v, length(v), 1))

for E in (SymPyPythonCall.Sym, SymPyPythonCall.Sym{SymPyPythonCall.PythonCall.Py})
    # Written out, one `return` per line: the coverage of a one-line method
    # generated by `@eval` was attributed to one line or the other from run to
    # run.
    @eval begin
        function SymPyPythonCall.SymPyCore.:↓(M::AbstractMatrix{$E})
            return _sympy_matrix(M)
        end
        # A vector is a column, as SymPyCore makes it. `Vector` as well as
        # `AbstractVector`: SymPyCore's method is on `Vector{<:Sym}`.
        function SymPyPythonCall.SymPyCore.:↓(v::AbstractVector{$E})
            return _sympy_column(v)
        end
        function SymPyPythonCall.SymPyCore.:↓(v::Vector{$E})
            return _sympy_column(v)
        end
    end
end

end # module
