# ── Searching the symmetry class of a tensor ───────────────────────────────────
#
# `best_sym_tens` tries the symmetry classes from the most to the least
# restrictive with the projections of `tens_projection.jl`, and returns the
# first one the tensor belongs to within a tolerance.

project(A::AbstractArray, sym::Symbol) = project(A, Val(sym))

# The axis a TI container carries, the frame an orthotropic or cubic one does;
# for any other tensor, the candidate read off the eigenstructure of its array
# `A` — and, for the cube, the canonical frame (see `best_sym_tens`).
_axis_or_candidate(t::TensTI, A) = axis(t)
_axis_or_candidate(t::AbstractTens, A) = _candidate_TI_axis(A)
_frame_or_candidate(t::Union{TensOrtho, TensCubic}, A) = frame(t)
_frame_or_candidate(t::AbstractTens, A) = _candidate_ORTHO_frame(A)
_cube_frame_or_canonical(t::Union{TensOrtho, TensCubic}, A) = frame(t)
_cube_frame_or_canonical(t::AbstractTens, A) = CanonicalBasis{3, eltype(A)}()

# The classes a fixed orientation can serve: a frame serves all but `:TI`, an
# axis all but `:ORTHO` and `:CUBIC` (see `best_sym_tens(t, n_or_frame)`).
_classes_served(::OrthonormalBasis{3}, proj) = filter(!=(:TI), proj)
_classes_served(n, proj) = filter(sym -> sym ∉ (:ORTHO, :CUBIC), proj)

# Internal helper: loop over symmetry projections and return the first match.
# `proj_fn(sym)` must return `(projected, d, drel)`.
function _best_sym_loop(newt::AbstractTens{order, dim, T}, proj, ε, proj_fn) where {order, dim, T}
    for sym in proj
        (projt, d, drel) = proj_fn(sym)
        if iszero(d) || drel < ε
            return projt, d, drel, sym
        end
    end
    return newt, zero(T), zero(T), :ANISO
end

for order in (2, 4)
    @eval begin
        """
            best_sym_tens(t; proj=(:ISO, :CUBIC, :TI, :ORTHO), ε=1e-6, optimize_angles=false)

        Find the best (most restrictive) symmetry of tensor `t` by trying each
        symmetry class in `proj` (from most to least symmetric) and accepting the
        first whose relative projection error is below `ε`.

        - `optimize_angles=false` (**default**, cheap path, no optimization): the
          `:ISO` projection is closed-form; for `:TI` the symmetry axis is taken
          from `t` itself (if it is a structured TI container) or derived from the
          Kelvin-Mandel eigenstructure (otherwise); for `:ORTHO` the material frame
          is taken from `t` or derived likewise.
        - `optimize_angles=true`: the `:TI` axis and the `:ORTHO` and `:CUBIC`
          frames are searched by [`project`](@ref) (a screening of fixed
          orientations, then Newton's method on the best ones).

        Returns `(projected, d, drel, sym)` where
        `sym ∈ {:ISO, :CUBIC, :TI, :ORTHO, :ANISO}`.

        # Examples
        ```julia
        julia> n = [0., 0., 1.];

        julia> C = ti_from_constants(10., 3., 2.5, 12., 2., n);

        julia> _, _, _, sym = best_sym_tens(C);

        julia> sym === :TI
        true
        ```

        See also [`best_sym_tens(t, n_or_frame)`](@ref) for fixed-axis/frame use,
        [`project`](@ref).
        """
        function best_sym_tens(
                t::AbstractTens{$order, dim, T};
                proj = (:ISO, :CUBIC, :TI, :ORTHO),
                ε = 1.0e-6,
                optimize_angles::Bool = false,
            ) where {dim, T}
            basis = relevant_OrthonormalBasis(get_basis(t))
            newt = express(t, basis)
            A = Array(get_array(newt))
            # `A` is read in the orthonormal frame `basis`; `Q` takes that frame
            # to the canonical one (`nothing` when it is canonical). Candidates
            # and projections live in that frame, and only the axis or frame of
            # the result is brought back: the 81 components are never rotated.
            Q = _rotation_to_canon(basis)
            if optimize_angles
                return _best_sym_loop(newt, proj, ε, sym -> _from_frame(Q, project(A, sym)))
            else
                # Cheap path: the axis or frame a structured container carries,
                # otherwise a candidate from the eigendecomposition of the trace
                # tensor of `A`, the array already materialized above — each
                # computed only if `proj` asks for the symmetry that reads it.
                n_default = :TI in proj ? _axis_or_candidate(t, A) : nothing
                frame_default = :ORTHO in proj ? _frame_or_candidate(t, A) : nothing
                # `:CUBIC` cannot use the orthotropic eigencandidate: every
                # second-order contraction of a cubic tensor is isotropic, so no
                # eigenframe says where the cube points, and that candidate
                # returns an arbitrary frame. A structured tensor carries its
                # own; anything else falls back to the canonical frame, which
                # finds an axis-aligned cube and nothing else. Detecting a
                # *rotated* cubic tensor needs `optimize_angles = true`.
                cube_default = :CUBIC in proj ? _to_frame(Q, _cube_frame_or_canonical(t, A)) : nothing
                return _best_sym_loop(
                    newt, proj, ε,
                    function (sym)
                        if sym === :ISO
                            return project(A, sym)
                        elseif sym === :TI
                            return _from_frame(Q, project(A, sym, n_default))
                        elseif sym === :CUBIC
                            return _from_frame(Q, project(A, sym, cube_default))
                        else    # :ORTHO
                            return _from_frame(Q, project(A, sym, frame_default))
                        end
                    end,
                )
            end
        end

        """
            best_sym_tens(t, n_or_frame; proj=(:ISO, :CUBIC, :TI, :ORTHO), ε=1e-6)

        Find the best symmetry of tensor `t` with a **fixed** symmetry axis `n`
        (for TI) or material frame `frame` (for ORTHO).  No rotation optimization
        is performed.

        - `n_or_frame`: a vector (axis for TI) or `OrthonormalBasis{3}` (frame for
          ORTHO and CUBIC). For ISO projection the extra argument is ignored.

        Because there is one argument for two kinds of orientation, the classes
        it cannot serve are **dropped from `proj`**: an axis leaves `:ISO` and
        `:TI`, a frame leaves `:ISO`, `:CUBIC` and `:ORTHO`. Pass the kind that
        matches the classes you mean to test.

        Returns `(projected, d, drel, sym)`.

        # Examples
        ```julia
        julia> n = [0., 0., 1.];

        julia> C = ti_from_constants(10., 3., 2.5, 12., 2., n);

        julia> _, _, drel, sym = best_sym_tens(C, n);

        julia> sym == :TI && drel < 1e-12
        true
        ```
        """
        function best_sym_tens(
                t::AbstractTens{$order, dim, T},
                n_or_frame;
                proj = (:ISO, :CUBIC, :TI, :ORTHO),
                ε = 1.0e-6,
            ) where {dim, T}
            basis = relevant_OrthonormalBasis(get_basis(t))
            newt = express(t, basis)
            A = Array(get_array(newt))
            # `A` is read in the orthonormal frame `basis`; `Q` takes that frame
            # to the canonical one (`nothing` when it is canonical). Candidates
            # and projections live in that frame, and only the axis or frame of
            # the result is brought back: the 81 components are never rotated.
            Q = _rotation_to_canon(basis)
            # `:TI` wants an axis, `:ORTHO` and `:CUBIC` a frame, and there is
            # one argument for both. Rather than let the mismatch surface as a
            # `MethodError` from deep inside `project` — which is what
            # happened before `:CUBIC` existed, whenever `:ORTHO` was reached
            # with an axis — the classes the argument cannot serve are dropped,
            # and the docstring says so.
            usable = _classes_served(n_or_frame, proj)
            local_arg = _to_frame(Q, n_or_frame)
            return _best_sym_loop(
                newt, usable, ε,
                sym -> sym == :ISO ? project(A, sym) :
                    _from_frame(Q, project(A, sym, local_arg))
            )
        end
    end
end

export project, best_sym_tens
