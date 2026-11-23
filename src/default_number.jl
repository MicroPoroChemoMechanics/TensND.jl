# =============================================================================
#  default_number.jl — the element type of what is built without one.
#
#  A few constructors take their element type as an optional argument:
#  `tens_Id2`, `tens_Id4`, `tens_J4`, `tens_K4`, `iso_projectors`, `𝐞`,
#  `levi_civita`, `CanonicalBasis()`. Their default is SymPy's `Sym` when the
#  package that makes `Sym` usable, SymPyPythonCall, is loaded (it brings the
#  Python backend of SymPyCore), and `Float64` otherwise, so that these calls
#  work in a session without Python.
# =============================================================================

"""
    default_number_type() -> Type

The element type of the tensors, vectors and bases built without one
(`tens_Id4()`, `tens_J4()`, `𝐞(1)`, `levi_civita()`, `CanonicalBasis()`, …):
`Sym` once SymPyPythonCall is loaded, `Float64` otherwise.

Loading SymPyPythonCall later in a session changes it. Code that needs a given
type passes it, `tens_J4(3, Float64)`.

# Examples
```jldoctest
julia> default_number_type()    # SymPyPythonCall is loaded in this session
Sym{PythonCall.Py}
```
"""
default_number_type() = _default_number(nothing)

# `Float64` here; the extension of SymPyPythonCall adds the more specific method
# on `Nothing`, which returns the concrete type of SymPy's symbols. A method
# rather than a global value, so that the type is a constant the compiler
# folds: the calls without an element type stay inferred and cost nothing more.
_default_number(::Any) = Float64

export default_number_type
