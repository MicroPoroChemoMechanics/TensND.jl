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

The choice is made at each call: loading SymPyPythonCall later in a session
changes it. Code that needs a given type passes it, `tens_J4(3, Float64)`,
which is also type-stable.

# Examples
```jldoctest
julia> default_number_type()    # SymPyPythonCall is loaded in this session
Sym{PythonCall.Py}
```
"""
default_number_type() = _DEFAULT_NUMBER[]

# Set by the extension of SymPyPythonCall when it is loaded, to the concrete
# type of its symbols (`Sym` with its Python object), so that the tensors built
# by default have a concrete element type.
const _DEFAULT_NUMBER = Ref{DataType}(Float64)

export default_number_type
