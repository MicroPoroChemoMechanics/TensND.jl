# ── Displaying tensors ───────────────────────────────────────────────────────

for OP in (:show, :print, :display)
    @eval begin
        function Base.$OP(t::AbstractTens)
            $OP(typeof(t))
            if ndims(t) == 4
                print("→ KM: ")
                $OP(KM(t))
            else
                print("→ array: ")
                $OP(get_array(t))
            end
            print("→ basis: ")
            $OP(vecbasis(get_basis(t)))
            print("→ var: ")
            return $OP(get_var(t))
        end
        function Base.$OP(t::TensOrthonormal)
            $OP(typeof(t))
            if ndims(t) == 4
                print("→ KM: ")
                $OP(KM(t))
            else
                print("→ array: ")
                $OP(get_array(t))
            end
            print("→ basis: ")
            return $OP(vecbasis(get_basis(t)))
        end

        # Base.$OP(t::AbstractTens{order,dim,T}; vec = '𝐞', coords = ntuple(i -> i, dim)) where {order,dim,T} = pprint(t; vec= vec, coords = coords)
    end
end

# Prevent SymPy's ImmutableMatrix from receiving 4D Julia arrays as elements,
# which triggers a SymPyDeprecationWarning for non-Expr objects in a Matrix.
# Redirect to the 2D Kelvin-Mandel representation instead.
#
# Restricted to symbolic element types **on purpose**: the fallback exists only
# because SymPy supplies a `text/latex` show for its own matrices. For a numeric
# tensor `KM(t)` is an ordinary `Matrix{Float64}`, which has no such method, so
# an unrestricted signature turned every numeric order-4 tensor into a
# `MethodError` in any front-end that asks for `text/latex` — Documenter and
# Jupyter both do. Numeric tensors now fall through to the plain-text display.
function Base.show(io::IO, M::MIME"text/latex", t::AbstractTens{4, dim, T}) where {dim, T <: SymType}
    return show(io, M, KM(t))
end

"""
    pprint(t::AbstractTens; vec = '𝐞', coords = 1:dim)
    pprint(t::AbstractTens, CS::AbstractCoorSystem; vec = '𝐞')

Print `t` **expanded on its basis** — as a sum of basis products such as
`(σʳʳ)𝐞ʳ⊗𝐞ʳ + (σᶿᶿ)𝐞ᶿ⊗𝐞ᶿ` — instead of as an array of components.

This is a display helper: it writes to `stdout` and returns `nothing`. For
symbolic results it is far more readable than the component array, and it is the
natural way to report a field computed in a curvilinear frame. Zero components
are omitted, and a vanishing tensor prints as `0`.

Passing a coordinate system names the basis vectors after its coordinates;
[`@set_coorsys`](@ref) makes that the default, so `pprint(t)` alone then
prints `𝐞ʳ`, `𝐞ᶿ`, … rather than `𝐞₁`, `𝐞₂`, ….

!!! note "Formerly `intrinsic`, then `print_tensor`"
    `intrinsic` was misleading: what is printed is explicitly *basis dependent*
    — it changes with the chart — whereas "intrinsic" denotes the basis-free
    form. `print_tensor` was misleading in turn, because the fallback below
    accepts a **scalar**, and `pprint` is the name MPCM already uses for a
    generic pretty-printer (`ChemistryLab.pprint`). Both former names still work
    and forward here, with a deprecation warning.

# Examples
```julia
julia> Spherical = coorsys_spherical() ; 𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ = unitvec(Spherical) ;

julia> @set_coorsys Spherical

julia> pprint(𝐞ʳ ⊗ 𝐞ʳ + 2 * 𝐞ᶿ ⊗ 𝐞ᶿ)   # spherical coordinates are ordered (θ, ϕ, r)
(2)𝐞ᶿ⊗𝐞ᶿ + 𝐞ʳ⊗𝐞ʳ
```

Applied to anything that is **not** a tensor — a scalar field, typically, so
that `LAPLACE(f(r)) |> pprint` works in a pipeline — it falls back to the
rich `text/plain` display, which for a symbolic expression is SymPy's
two-dimensional form:

```julia
julia> LAPLACE(SymFunction("f", real = true)(r)) |> pprint
              d
 2          2⋅──(f(r))
d             dr
───(f(r)) + ──────────
  2             r
dr
```

Tensor coefficients are deliberately *not* printed that way: a multi-line
coefficient interpolated into `(…)𝐞ʳ⊗𝐞ʳ + (…)𝐞ᶿ⊗𝐞ᶿ` would destroy the layout
that makes the expanded form readable in the first place.

See also [`get_array`](@ref), [`components`](@ref).
"""
function pprint(t)
    # The three-argument `show`, not `println`: SymPy's two-dimensional form is
    # reachable only through `MIME"text/plain"`, which is what the REPL calls.
    # `println` and `string` go through the two-argument `show` and give the
    # flat form. `sympy.init_printing` does not change any of this — it is a
    # no-op for all three paths in current SymPy.jl.
    show(stdout, MIME"text/plain"(), t)
    return println()
end

function pprint(t::AbstractTens{order, dim, T}; vec = '𝐞', coords = ntuple(i -> i, dim)) where {order, dim, T}
    return println(_pprint_string(t; vec = vec, coords = coords))
end

# Formatting split from printing, so the layout can be tested against a string
# instead of capturing `stdout`.
function _pprint_string(t::AbstractTens{order, dim, T}; vec = '𝐞', coords = ntuple(i -> i, dim)) where {order, dim, T}
    ind = CartesianIndices(t)
    ℬ = get_basis(t)
    firstprint = true
    s = ""
    for i in ind
        x = t[i]
        if !iszero(x)
            if !firstprint
                s *= " + "
            end
            if !isone(x)
                s *= "(" * string(x) * ")"
            end
            j = Tuple(i)
            for k in 1:order
                s *= strvecbasis(ℬ, coords[j[k]], invvar(get_var(t, k)); vec = vec)
                if k < order
                    s *= "⊗"
                end
            end
            firstprint = false
        end
    end
    return length(s) > 0 ? s : "0"
end


@deprecate intrinsic(t; kwargs...) pprint(t; kwargs...)
@deprecate intrinsic(t, CS; kwargs...) pprint(t, CS; kwargs...)
@deprecate print_tensor(t; kwargs...) pprint(t; kwargs...)
@deprecate print_tensor(t, CS; kwargs...) pprint(t, CS; kwargs...)

export pprint
