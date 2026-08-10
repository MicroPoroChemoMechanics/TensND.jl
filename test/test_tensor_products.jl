# =============================================================================
#  test_tensor_products.jl — the outer products, pinned against their definition.
#
#  `otimes`, `otimesu`, `otimesl` and `sotimes` on `AbstractArray` were rewritten
#  from `einsum` calls to single broadcasts (v0.3.3). They contract nothing —
#  their output indices are the union of the input indices — so the rewrite is
#  purely about how elements are placed, and a wrong placement is a silent
#  numerical error, not an exception.
#
#  Three independent references, because one is not enough:
#
#  1. the **definition** written out as explicit loops, which owes nothing to
#     any implementation;
#  2. the **einsum** form that was replaced, kept here verbatim as the historical
#     reference;
#  3. **algebraic identities** the products must satisfy whatever the code.
#
#  And, deliberately, coverage over the *types the callers actually pass*.
#  The v0.3.3 regression that reached the documentation build — `transpose` of a
#  `Tensors.Vec` is discontinued upstream — passed every array-based check and
#  the whole suite, because nothing here exercised a first-order `Tensors`
#  operand. That gap is the reason for the "Tensors element types" section.
# =============================================================================

# ── Reference 1: the definitions, as explicit loops ──────────────────────────
#
# `out[…] = t1[…] * t2[…]` with each operand's m-th index at output position
# `ec[m]`. Nothing clever: this is the specification.

function _ref_placed(t1, t2, ec1, ec2)
    n = length(ec1) + length(ec2)
    sz = ntuple(n) do k
        m = findfirst(isequal(k), ec1)
        m === nothing ? size(t2, findfirst(isequal(k), ec2)) : size(t1, m)
    end
    out = Array{promote_type(eltype(t1), eltype(t2))}(undef, sz)
    for I in CartesianIndices(out)
        i1 = ntuple(m -> I[ec1[m]], length(ec1))
        i2 = ntuple(m -> I[ec2[m]], length(ec2))
        out[I] = t1[i1...] * t2[i2...]
    end
    return out
end

_ec_otimes(o1, o2) = (ntuple(i -> i, o1), ntuple(i -> o1 + i, o2))
_ec_otimesu(o1, o2) = ((ntuple(i -> i, o1 - 1)..., o1 + 1), (o1, ntuple(i -> o1 + 1 + i, o2 - 1)...))
_ec_otimesl(o1, o2) = ((ntuple(i -> i, o1 - 1)..., o1 + 2), (o1, o1 + 1, ntuple(i -> o1 + 2 + i, o2 - 2)...))

_ref_otimes(a, b) = _ref_placed(a, b, _ec_otimes(ndims(a), ndims(b))...)
_ref_otimesu(a, b) = _ref_placed(a, b, _ec_otimesu(ndims(a), ndims(b))...)
_ref_otimesl(a, b) = _ref_placed(a, b, _ec_otimesl(ndims(a), ndims(b))...)
_ref_sotimes(a, b) = (_ref_otimes(a, b) .+ _ref_otimesu(a, b)) ./ 2

# ── Reference 2: the einsum implementations that were replaced ───────────────

function _einsum_placed(t1::AbstractArray{A}, t2::AbstractArray{B}, ec1, ec2) where {A, B}
    ec3 = ntuple(i -> i, length(ec1) + length(ec2))
    return einsum(
        EinCode((ec1, ec2), ec3), (AbstractArray{A}(t1), AbstractArray{B}(t2))
    )
end

_old_otimes(a, b) = _einsum_placed(a, b, _ec_otimes(ndims(a), ndims(b))...)
_old_otimesu(a, b) = _einsum_placed(a, b, _ec_otimesu(ndims(a), ndims(b))...)
_old_otimesl(a, b) = _einsum_placed(a, b, _ec_otimesl(ndims(a), ndims(b))...)
_old_sotimes(a, b) = (_old_otimes(a, b) + _old_otimesu(a, b)) / 2

@testsection "Tensor products (otimes family)" begin

    @testset "the convention, spelled out for order 2 ⊗ 2" begin
        # The four products differ only in which index of which operand goes
        # where. Writing that out once, by hand, is what stops a refactor from
        # silently permuting them.
        a = rand(3, 3)
        b = rand(3, 3)
        o, u, l = otimes(a, b), otimesu(a, b), otimesl(a, b)
        s = TensND.sotimes(a, b)
        for i in 1:3, j in 1:3, k in 1:3, m in 1:3
            @test o[i, j, k, m] == a[i, j] * b[k, m]
            @test u[i, j, k, m] == a[i, k] * b[j, m]
            @test l[i, j, k, m] == a[i, m] * b[j, k]
            @test s[i, j, k, m] == (a[i, j] * b[k, m] + a[i, k] * b[j, m]) / 2
        end
    end

    @testset "against the definition, over shapes and orders" begin
        shapes = (
            ((3, 3), (3, 3)), ((2, 2), (2, 2)), ((4, 4), (4, 4)),
            ((2, 3), (3, 2)), ((3, 3), (3, 3, 3)), ((3, 3, 3), (3, 3)),
            ((3, 3, 3), (3, 3, 3)), ((2, 3, 4), (4, 3)),
        )
        for (s1, s2) in shapes
            a = rand(s1...)
            b = rand(s2...)
            # `===` and not `≈`: the rewrite must not even reassociate, it
            # places the same products in the same order.
            @test all(otimes(a, b) .=== _ref_otimes(a, b))
            @test all(otimesu(a, b) .=== _ref_otimesu(a, b))
            @test all(otimesl(a, b) .=== _ref_otimesl(a, b))
            @test TensND.sotimes(a, b) ≈ _ref_sotimes(a, b)
            @test size(otimes(a, b)) == (s1..., s2...)
        end

        # `otimes` and `otimesu` also accept a first-order second operand;
        # `otimesl` needs two indices there to displace.
        for (s1, s2) in (((3, 3), (3,)), ((3, 3, 3), (3,)), ((3,), (3,)))
            a = rand(s1...)
            b = rand(s2...)
            @test all(otimes(a, b) .=== _ref_otimes(a, b))
            @test all(otimesu(a, b) .=== _ref_otimesu(a, b))
        end
    end

    @testset "against the einsum implementation it replaced" begin
        for (s1, s2) in (
                ((3, 3), (3, 3)), ((2, 3), (3, 2)), ((3, 3), (3, 3, 3)),
                ((3, 3, 3), (3, 3)), ((4, 4), (4, 4)),
            )
            a = rand(s1...)
            b = rand(s2...)
            @test all(otimes(a, b) .=== _old_otimes(a, b))
            @test all(otimesu(a, b) .=== _old_otimesu(a, b))
            @test all(otimesl(a, b) .=== _old_otimesl(a, b))
            @test TensND.sotimes(a, b) ≈ _old_sotimes(a, b)
        end
    end

    @testset "Tensors element types — the first-order operand" begin
        # This is the case that broke the documentation build in v0.3.3 while
        # every array-based test stayed green: `transpose` of a `Tensors.Vec`
        # is discontinued upstream, so an implementation that transposed raised
        # here and nowhere else.
        v2 = Tensors.Vec{2}((1.0, 2.0))
        v3 = Tensors.Vec{3}((1.0, 2.0, 3.0))
        t2 = Tensors.Tensor{2, 3}((i, j) -> 1.0i + 3.0j)
        s2 = Tensors.SymmetricTensor{2, 3}((i, j) -> 1.0i + 1.0j)
        arr = rand(3, 3)

        for (x, y) in (
                (v2, v2), (v3, v3), (t2, t2), (s2, s2),
                (arr, v3), (v3, arr), (t2, arr),
            )
            p = otimes(x, y)
            @test size(p) == (size(x)..., size(y)...)
            for I in CartesianIndices(p)
                i1 = ntuple(m -> I[m], ndims(x))
                i2 = ntuple(m -> I[ndims(x) + m], ndims(y))
                @test p[I] ≈ x[i1...] * y[i2...]
            end
        end

        # …and the interleaving products on second-order Tensors arrays.
        for (x, y) in ((t2, t2), (t2, arr), (arr, t2))
            @test otimesu(x, y) ≈ _ref_otimesu(collect(x), collect(y))
            @test otimesl(x, y) ≈ _ref_otimesl(collect(x), collect(y))
        end
    end

    @testset "exact and generic element types" begin
        # Rational arithmetic is exact, so any misplacement shows as an exact
        # inequality rather than as a tolerance question.
        a = [1 // 2 1 // 3; 1 // 5 1 // 7]
        b = [2 // 3 3 // 5; 5 // 7 7 // 11]
        @test all(otimes(a, b) .== _ref_otimes(a, b))
        @test all(otimesu(a, b) .== _ref_otimesu(a, b))
        @test all(otimesl(a, b) .== _ref_otimesl(a, b))
        @test eltype(otimes(a, b)) <: Rational

        # Integers must not be silently promoted to Float64.
        ai = [1 2; 3 4]
        bi = [5 6; 7 8]
        @test eltype(otimes(ai, bi)) === Int
        @test all(otimesu(ai, bi) .== _ref_otimesu(ai, bi))

        # Forward-mode duals must survive: these products sit under every
        # sensitivity computation in MeanFieldHomogenization.
        d = ForwardDiff.Dual(2.0, 1.0)
        ad = fill(d, 2, 2)
        bd = rand(2, 2)
        pd = otimes(ad, bd)
        @test eltype(pd) <: ForwardDiff.Dual
        @test ForwardDiff.value(pd[1, 1, 1, 1]) ≈ 2.0 * bd[1, 1]
        @test ForwardDiff.partials(pd[1, 1, 1, 1], 1) ≈ bd[1, 1]
        @test all(otimesu(ad, bd) .≈ _ref_otimesu(ad, bd))

        # And symbolic elements, where no BLAS path could ever be taken.
        @syms x::real y::real
        as = [x 0; 0 y]
        bs = [y x; x y]
        @test simplify(otimes(as, bs)[1, 1, 1, 2] - x * x) == 0
        @test simplify(otimesu(as, bs)[1, 1, 1, 2] - x * x) == 0
        # `sotimes` divides by 2 rather than multiplying by 0.5, which keeps a
        # symbolic result readable; both are bit-identical in Float64.
        @test simplify(TensND.sotimes(as, bs)[1, 1, 1, 1] - x * y) == 0
    end

    @testset "algebraic identities" begin
        𝟏 = Matrix(1.0I, 3, 3)
        m = [1.0 0.4 0.2; 0.1 2.0 0.3; 0.5 0.6 3.0]     # deliberately not symmetric

        # Documented in the `otimesu` docstring: 𝟏 ⊠ 𝟏 is the identity of ALL
        # order-2 tensors, where ⊠ˢ gives the identity of the symmetric ones.
        @test (𝟏 ⊠ 𝟏) ⊡ m ≈ m
        @test (𝟏 ⊠ˢ 𝟏) ⊡ m ≈ (m + m') / 2

        # Also documented: ⊠ inverts termwise, unlike its symmetrized partner.
        a = rand(3, 3) + 3I
        b = rand(3, 3) + 3I
        @test (a ⊠ b) ⊡ (inv(a) ⊠ inv(b)) ≈ 𝟏 ⊠ 𝟏

        # The symmetrized product of two vectors is symmetric by construction.
        u = Tensors.Vec{3}((1.0, 2.0, 3.0))
        w = Tensors.Vec{3}((4.0, 5.0, 6.0))
        su = TensND.sotimes(u, w)
        @test su ≈ (otimes(u, w) + otimes(w, u)) / 2
        @test su[1, 2] ≈ su[2, 1]

        # ⊠ and ⊠ˡ differ by a transposition of the second operand's indices.
        x = rand(3, 3)
        y = rand(3, 3)
        @test otimesl(x, y) ≈ permutedims(otimesu(x, y), (1, 2, 4, 3))
    end

    @testset "no allocation beyond the result" begin
        # The rewrite's point: one broadcast, one array. A future change that
        # reintroduces an intermediate shows up here rather than in a profile.
        a = rand(3, 3)
        b = rand(3, 3)
        result_bytes = sizeof(Float64) * 81
        for f in (otimes, otimesu, otimesl)
            @test (@allocated f(a, b)) < 3 * result_bytes
        end
        # `sotimes` sums two placements; still a single output array.
        @test (@allocated TensND.sotimes(a, b)) < 3 * result_bytes
    end
end

# =============================================================================
#  best_sym_tens: restricting `proj` must not change the answers it does give.
#
#  The candidate axis and frame are now derived only when `proj` asks for the
#  symmetry that reads them, and from the array already materialized rather
#  than from two fresh copies. Laziness that changed a result would be a silent
#  behavior change, so the restricted calls are compared with the full one.
#
#  Note the nesting ISO ⊂ TI ⊂ ORTHO: an isotropic tensor is also transversely
#  isotropic, so restricting `proj` legitimately changes *which* symmetry is
#  reported. What must not change are the numbers for a given symmetry.
# =============================================================================

@testsection "best_sym_tens — restricting proj" begin
    Ciso = TensISO{3}(3 * 10.0, 2 * 4.0)
    Cti = tens_TI(10.0, 3.0, 2.5, 12.0, 2.0, [0.0, 0.0, 1.0])

    @testset "the most restrictive symmetry wins when all are offered" begin
        @test best_sym_tens(Ciso)[4] === :ISO
        @test best_sym_tens(Cti)[4] === :TI
    end

    @testset "a restricted proj gives the same numbers for that symmetry" begin
        for C in (Ciso, Cti), sym in (:TI, :ORTHO)
            # Both tensors are at least TI, so TI and ORTHO always fit.
            r = best_sym_tens(C; proj = (sym,))
            @test r[4] === sym
            @test r[3] < 1.0e-12
            # Asking for the same symmetry twice must be deterministic, which
            # is what the lazy candidate could have broken.
            r2 = best_sym_tens(C; proj = (sym,))
            @test r[1] ≈ r2[1]
            @test r[2] == r2[2]
            @test r[3] == r2[3]
        end
        # An isotropic tensor is isotropic; a TI one is not.
        @test best_sym_tens(Ciso; proj = (:ISO,))[4] === :ISO
        @test best_sym_tens(Cti; proj = (:ISO,))[4] === :ANISO
    end

    @testset "the full call agrees with the restricted one it selected" begin
        for C in (Ciso, Cti)
            full = best_sym_tens(C)
            restricted = best_sym_tens(C; proj = (full[4],))
            @test restricted[4] === full[4]
            @test restricted[1] ≈ full[1]
            @test restricted[3] ≈ full[3]
        end
    end
end

# =============================================================================
#  Arbitrary orders, and variance.
#
#  TensND is an N-dimensional package: these products must work for any pair of
#  orders, which is what the `einsum` calls provided generically. The broadcast
#  that replaced them relies on a property that has to hold at every order —
#  that both index lists are increasing, so neither operand's own axes need
#  reordering — so that property is checked here rather than assumed.
#
#  The contraction paths (`contract`, `dcontract`, `dotdot`, `qcontract`) still
#  go through `einsum` and were deliberately left alone; they are exercised too,
#  so that a future attempt to give them the same treatment cannot land
#  unnoticed.
# =============================================================================

@testsection "Arbitrary orders and variance" begin

    @testset "outer products at every order pair" begin
        for o1 in 1:5, o2 in 1:5
            o1 + o2 > 8 && continue          # keeps the arrays a sane size
            a = rand(ntuple(_ -> 3, o1)...)
            b = rand(ntuple(_ -> 3, o2)...)

            @test all(otimes(a, b) .=== _ref_placed(a, b, _ec_otimes(o1, o2)...))
            @test all(otimesu(a, b) .=== _ref_placed(a, b, _ec_otimesu(o1, o2)...))
            @test ndims(otimes(a, b)) == o1 + o2
            if o2 ≥ 2                         # otimesl needs two indices to displace
                @test all(otimesl(a, b) .=== _ref_placed(a, b, _ec_otimesl(o1, o2)...))
            end

            # The property the implementation rests on: with both index lists
            # increasing, singleton axes suffice and no permutation is needed.
            increasing(t) = all(t[i] < t[i + 1] for i in 1:(length(t) - 1))
            @test all(increasing, _ec_otimesu(o1, o2))
            o2 ≥ 2 && @test all(increasing, _ec_otimesl(o1, o2))
        end
    end

    @testset "contractions still reduce the order correctly" begin
        for (o1, o2) in ((5, 3), (3, 5), (4, 2), (2, 4), (5, 2), (3, 3), (4, 4))
            a = rand(ntuple(_ -> 3, o1)...)
            b = rand(ntuple(_ -> 3, o2)...)
            @test ndims(dcontract(a, b)) == o1 + o2 - 4
        end
        @test ndims(TensND.contract(rand(3, 3, 3, 3, 3), 2, 4)) == 3
    end

    @testset "variance is carried, on a non-orthonormal basis" begin
        # A non-orthonormal basis is what makes co- and contravariance
        # observable at all: on an orthonormal one the two coincide.
        ℬ = Basis(rand(3, 3) + 2I)
        for (o1, o2) in ((2, 1), (2, 2), (3, 2), (5, 3), (2, 3))
            v1 = ntuple(_ -> rand((:cov, :cont)), o1)
            v2 = ntuple(_ -> rand((:cov, :cont)), o2)
            t1 = Tens(rand(ntuple(_ -> 3, o1)...), ℬ, v1)
            t2 = Tens(rand(ntuple(_ -> 3, o2)...), ℬ, v2)

            p = t1 ⊗ t2
            @test get_var(p) == (v1..., v2...)
            # Basis-independent invariant: the canonical components of the
            # product are the product of the canonical components.
            @test Array(components_canon(p)) ≈ _ref_placed(
                Array(components_canon(t1)), Array(components_canon(t2)),
                _ec_otimes(o1, o2)...
            )

            # ⊠ interleaves the indices, so it must interleave the variance the
            # same way — a mismatch there would be silent on an orthonormal
            # basis and wrong on any other.
            pu = t1 ⊠ t2
            @test get_var(pu) == (v1[1:(end - 1)]..., v2[1], v1[end], v2[2:end]...)
            @test Array(components_canon(pu)) ≈ _ref_placed(
                Array(components_canon(t1)), Array(components_canon(t2)),
                _ec_otimesu(o1, o2)...
            )
        end
    end
end
