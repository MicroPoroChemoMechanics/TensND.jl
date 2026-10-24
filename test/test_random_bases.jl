# =============================================================================
#  test_random_bases.jl — every operation is the same whatever basis its
#  operands are written in.
#
#  `get_array(t)` and `t[i, j, …]` read a tensor in its own basis, while an axis
#  or a frame is given in canonical coordinates. Code tested only with
#  canonical or structured tensors (which answer in the canonical basis) looks
#  right and is wrong for a `Tens` written in a rotated basis. Round angles hide
#  it, since a component read in the wrong frame can coincide with the right one
#  there, so every frame below is drawn at random, and each operand gets its own.
#
#  The reference of each check is the same tensor written in the canonical
#  basis, its components rotated by a comprehension written here (not by
#  `change_tens`), or an explicit sum over indices.
# =============================================================================

const RB_SEED = 20261021
const RB_DRAWS = 5

"A uniform random rotation (QR of a Gaussian matrix, det fixed to +1)."
function _rb_rotation(rng)
    F = qr(randn(rng, 3, 3))
    Q = Matrix(F.Q) * Diagonal(sign.(diag(F.R)))
    return det(Q) < 0 ? Q * Diagonal([1.0, 1.0, -1.0]) : Q
end

"A random unit vector."
_rb_unit(rng) = normalize(randn(rng, 3))

"A random minor- and major-symmetric 4th-order array (from a symmetric 6×6)."
function _rb_sym4(rng)
    M = randn(rng, 6, 6)
    return TensND.array_from_mandel66(M + transpose(M) + 12I)
end

"Components in the frame whose vectors are the columns of `R`."
_rb_in_frame(A::AbstractArray{<:Any, 4}, R) = [
    sum(R[i, a] * R[j, b] * R[k, c] * R[l, d] * A[i, j, k, l] for i in 1:3, j in 1:3, k in 1:3, l in 1:3)
        for a in 1:3, b in 1:3, c in 1:3, d in 1:3
]
_rb_in_frame(A::AbstractArray{<:Any, 2}, R) =
    [sum(R[i, a] * R[j, b] * A[i, j] for i in 1:3, j in 1:3) for a in 1:3, b in 1:3]
_rb_in_frame(v::AbstractVector, R) = transpose(R) * v

"The tensor of canonical components `A`, written in the rotated basis of `R`."
_rb_rotated(A, R) = Tens(_rb_in_frame(A, R), RotatedBasis(R))

# `AbstractTens <: AbstractArray`: the tensor method must be the more specific
# one, or `Array(t)` would read `t` in its own basis, the very defect tested.
_rb_canon(t::TensND.AbstractTens) = Array(components_canon(t))
_rb_canon(A::AbstractArray) = Array(A)

@testset "random bases — draw $k" for k in 1:RB_DRAWS
    rng = Random.Xoshiro(RB_SEED + k)
    R = _rb_rotation(rng)                  # basis of the tensor
    n = _rb_unit(rng)                      # an axis, canonical coordinates
    F = RotatedBasis(_rb_rotation(rng))    # a frame, canonical coordinates

    C = _rb_sym4(rng)
    tc = Tens(C)                           # canonical basis
    tr_ = _rb_rotated(C, R)                # the same tensor, rotated basis
    @test _rb_canon(tr_) ≈ C               # premise

    @testset "fixed-axis and fixed-frame projections" begin
        for (args, sym) in (((), :ISO), ((n,), :TI), ((F,), :ORTHO), ((F,), :CUBIC))
            Pc, dc, rc = proj_tens(sym, tc, args...)
            Pr, dr, rr = proj_tens(sym, tr_, args...)
            @test Pr ≈ Pc
            @test dr ≈ dc
            @test rr ≈ rc
        end
        # Order 2, with an antisymmetric part.
        a = randn(rng, 3, 3)
        sc, sr = Tens(a), _rb_rotated(a, R)
        for (args, sym) in (((), :ISO), ((n,), :TI), ((F,), :ORTHO))
            Pc, dc, _ = proj_tens(sym, sc, args...)
            Pr, dr, _ = proj_tens(sym, sr, args...)
            @test Pr ≈ Pc
            @test dr ≈ dc
        end
    end

    @testset "an exactly TI tensor is recognized about its own axis" begin
        ti = Array(get_array(tens_ti(10.0, 3.0, 2.5, 12.0, 2.0, n)))
        t = _rb_rotated(ti, R)
        @test istransverselyisotropic(t, n)
        P, d, drel = proj_tens(:TI, t, n)
        @test drel < 1.0e-12
        @test _rb_canon(P) ≈ ti
        # The search, with no axis: the same tensor, whatever axis sign it finds.
        Ps, _, drels = proj_tens(:TI, t)
        @test drels < 1.0e-8
        @test _rb_canon(Ps) ≈ ti rtol = 1.0e-8
        # best_sym_tens returns the tensor itself, not one about another axis.
        Pb, _, _, sym = best_sym_tens(t)
        @test sym === :TI
        @test _rb_canon(Pb) ≈ ti rtol = 1.0e-8
        Pf, _, _, symf = best_sym_tens(t, n)
        @test symf === :TI
        @test _rb_canon(Pf) ≈ ti
    end

    @testset "an isotropic tensor written in a rotated basis" begin
        iso = Array(get_array(TensISO{3}(3 * 7.0, 2 * 3.0)))
        t = _rb_rotated(iso, R)
        # The angle search is the path that brings each candidate back from the
        # tensor's frame, the isotropic one included.
        Pb, _, _, sym = best_sym_tens(t; optimize_angles = true)
        @test sym === :ISO
        @test Pb isa TensISO
        @test _rb_canon(Pb) ≈ iso
    end

    @testset "an exactly orthotropic tensor about a random frame" begin
        o = Array(get_array(TensOrtho(10.0, 8.0, 12.0, 3.0, 2.5, 1.5, 2.0, 3.0, 3.5, F)))
        t = _rb_rotated(o, R)
        @test isorthotropic(t, F)
        P, _, drel = proj_tens(:ORTHO, t, F)
        @test drel < 1.0e-12
        @test _rb_canon(P) ≈ o
        Pb, _, _, sym = best_sym_tens(t, F; proj = (:ISO, :CUBIC, :ORTHO))
        @test sym === :ORTHO
        @test _rb_canon(Pb) ≈ o
    end

    @testset "rotational averages" begin
        @test isotropify(tr_) ≈ isotropify(tc)
        @test transverse_isotropify(tr_, n) ≈ transverse_isotropify(tc, n)
        a = randn(rng, 3, 3)
        @test transverse_isotropify(_rb_rotated(a, R), n) ≈
            transverse_isotropify(Tens(a), n)
        # A left-handed frame: the antisymmetric part `c` of the order-2 average
        # comes from a cross product, which a reflection turns over.
        L = R * Diagonal([1.0, 1.0, -1.0])
        @test transverse_isotropify(_rb_rotated(a, L), n) ≈
            transverse_isotropify(Tens(a), n)
        @test transverse_isotropify(_rb_rotated(C, L), n) ≈ transverse_isotropify(tc, n)
        # A non-orthonormal basis: the traces need the metric.
        M = R + 0.3 * randn(rng, 3, 3)
        tg = change_tens(tc, Basis(M, Val(:cov)))
        @test isotropify(tg) ≈ isotropify(tc)
    end

    @testset "dotdot, each operand in its own basis" begin
        R1, R2 = _rb_rotation(rng), _rb_rotation(rng)
        u, v = randn(rng, 3), randn(rng, 3)
        ref = [sum(u[i] * C[i, j, k, l] * v[l] for i in 1:3, l in 1:3) for j in 1:3, k in 1:3]
        U = Tens(_rb_in_frame(u, R1), RotatedBasis(R1))
        V = Tens(_rb_in_frame(v, R2), RotatedBasis(R2))
        @test _rb_canon(dotdot(U, tr_, V)) ≈ ref
        @test _rb_canon(dotdot(U, tc, V)) ≈ ref
        # A structured operand: its data is a coefficient tuple, not components.
        ti = tens_ti(10.0, 3.0, 2.5, 12.0, 2.0, n)
        T4 = Array(get_array(ti))
        refti = [sum(u[i] * T4[i, j, k, l] * v[l] for i in 1:3, l in 1:3) for j in 1:3, k in 1:3]
        @test _rb_canon(dotdot(U, ti, V)) ≈ refti
    end

    @testset "(α𝕁 + β𝕂) ⊡ B for a non-symmetric B" begin
        α, β = 3.0 + rand(rng), 1.0 + rand(rng)
        A = TensISO{3}(α, β)
        AA = Array(get_array(A))
        b = randn(rng, 3, 3)
        ref = [sum(AA[i, j, k, l] * b[k, l] for k in 1:3, l in 1:3) for i in 1:3, j in 1:3]
        refl = [sum(b[i, j] * AA[i, j, k, l] for i in 1:3, j in 1:3) for k in 1:3, l in 1:3]
        for B in (Tens(b), _rb_rotated(b, R))
            @test _rb_canon(A ⊡ B) ≈ ref
            @test _rb_canon(B ⊡ A) ≈ refl
        end
    end

    @testset "an axis given as a vector in a rotated basis" begin
        nv = Tens(_rb_in_frame(n, R), RotatedBasis(R))
        @test collect(axis(tens_ti(10.0, 3.0, 2.5, 12.0, 2.0, nv))) ≈ n
    end
end

@testset "random bases — the canonical path is unchanged" begin
    rng = Random.Xoshiro(RB_SEED)
    C = _rb_sym4(rng)
    n = _rb_unit(rng)
    F = RotatedBasis(_rb_rotation(rng))
    t = Tens(C)
    # A tensor in the canonical basis is projected on its own stored components,
    # with no change of basis: the result is the array kernel's on `get_array(t)`
    # (what the wrappers passed before), to the bit.
    for (args, sym) in (((n,), :TI), ((F,), :ORTHO), ((F,), :CUBIC))
        Pt, dt, rt = proj_tens(sym, t, args...)
        Pa, da, ra = proj_tens(sym, get_array(t), args...)
        @test get_data(Pt) == get_data(Pa)
        @test (dt, rt) == (da, ra)
    end
    @test TensND._orthonormal_components(t)[2] === nothing
end

@testset "Kelvin-Mandel of a slightly non-minor-symmetric array" begin
    rng = Random.Xoshiro(RB_SEED + 100)
    C = _rb_sym4(rng)
    C[1, 2, 1, 1] += 1.0e-5                  # beyond the symmetry tolerance
    @test size(TensND._KM_of_array(C)) == (6, 6)
    @test size(TensND._KM_of_array(_rb_sym4(rng))) == (6, 6)
    n = _rb_unit(rng)
    _, _, drel = proj_tens(:TI, C, n)       # was a DimensionMismatch
    @test isfinite(drel)
end

@testset "RotatedBasis refuses a non-orthogonal matrix" begin
    @test_throws ArgumentError RotatedBasis([1.0 0.5 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
    rng = Random.Xoshiro(RB_SEED + 200)
    Q = _rb_rotation(rng)
    @test vecbasis(RotatedBasis(Q)) ≈ Q
end

@testset "mandel_rotation: the Kelvin-Mandel matrix in a rotated basis" begin
    rng = Random.Xoshiro(RB_SEED + 300)
    for _ in 1:RB_DRAWS
        R = _rb_rotation(rng)
        Q = mandel_rotation(R)
        @test transpose(Q) * Q ≈ I atol = 1.0e-12
        M = TensND.mandel66_minor(_rb_sym4(rng))
        # The components in the basis of `R` are Qᵀ M Q, as `mandel` in that basis says.
        @test transpose(Q) * M * Q ≈ Matrix(mandel(Tens(from_mandel(M)), RotatedBasis(R))) atol = 1.0e-10
    end
    # QᵀQ = I along a path of rotations, so its derivative vanishes.
    dQQ = ForwardDiff.derivative(θ -> (Q = mandel_rotation(θ, 0.4, 1.3); transpose(Q) * Q), 0.7)
    @test norm(dQQ) ≤ 1.0e-12
end

@testset "primal and carries_partials, at every nesting level" begin
    x = ForwardDiff.Dual{:a}(2.0, 1.0)
    y = ForwardDiff.Dual{:b}(x, ForwardDiff.Dual{:a}(0.0, 0.0))     # nested, outer partial zero
    @test primal(x) === 2.0
    @test primal(y) === 2.0
    @test primal(3) === 3
    @test carries_partials(x)
    @test carries_partials(y)                    # the inner level carries one
    @test !carries_partials(ForwardDiff.Dual{:a}(2.0, 0.0))
    @test !carries_partials(1.5)
    @syms s
    @test !carries_partials(s) && primal(s) === s
    # The float type under the duals, the depth, the channels.
    @test primal_type(typeof(y)) === Float64
    @test primal_type(ForwardDiff.Dual{:a, Float32, 1}) === Float32
    @test primal_type(Int) === Float64 && primal_type(Rational{Int}) === Float64
    @test primal_type(typeof(s)) === typeof(s)
    @test dual_depth(typeof(y)) == 2 && dual_depth(typeof(x)) == 1 && dual_depth(Float64) == 0
    # Value, then each partial recursively: x = 2 + 1ε, y = x + (0 + 0ε) ε'.
    @test dual_channels(1.5) == (1.5,)
    @test dual_channels(x) == (2.0, 1.0)
    @test dual_channels(y) == (2.0, 1.0, 0.0, 0.0)
    @test length(dual_channels(ForwardDiff.Dual{:c}(1.0, 2.0, 3.0))) == 3
    @test float_like(Int) === Float64 && float_like(Float32) === Float32
    @test float_like(typeof(s)) === typeof(s) && float_like(typeof(x)) === typeof(x)
end

@testset "isparallel, numeric, dual and symbolic" begin
    rng = Random.Xoshiro(RB_SEED + 400)
    n = _rb_unit(rng)
    @test isparallel(n, -3.0 .* n)
    @test isparallel(Tuple(2.0 .* n), n)
    @test !isparallel(n, _rb_unit(rng))
    @test !isparallel(n, zeros(3))
    # Read in canonical components, whatever basis the vector is written in.
    R = _rb_rotation(rng)
    @test isparallel(Tens(transpose(R) * n, RotatedBasis(R)), n)
    # A dual direction, compared on its value.
    @test isparallel(ForwardDiff.Dual{:t}.(n, 1.0), n)
    # Symbolic: structural, exact, never a false positive.
    @syms θ::real
    a = [sin(θ), cos(θ), 0]
    @test isparallel(a, [2sin(θ), 2cos(θ), 0])
    @test !isparallel(a, [cos(θ), sin(θ), 0])
end
