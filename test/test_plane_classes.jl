# =============================================================================
#  test_plane_classes.jl — the orthotropic and the square classes of the
#  plane: `TensOrtho{4, 2}` of four or five constants, `TensOrtho{2, 2}`, and
#  `TensCubic{2}`, with their algebra, their symmetry blocks, the promotions
#  from the isotropic class and the reductions of a tensor of space to one of
#  its planes of symmetry.
#
#  References, independent of the code under test: the canonical components
#  rebuilt from the Kelvin-Mandel matrix in the material frame by explicit
#  sums over the indices (`reference_contractions.jl`), matrix algebra on these Kelvin-Mandel
#  matrices, the sub-matrices of the 6×6 Kelvin-Mandel matrix of a tensor of
#  space (plane strain) and of its inverse (plane stress), the Kelvin-Mandel
#  rotation of the plane written out, closed forms of the rotated entries,
#  and the rotational average computed by quadrature over the angle.
# =============================================================================

using StaticArrays: SMatrix, SVector

const PLC_SEED = 20261128
const PLC_PAIRS = ((1, 1), (2, 2), (1, 2))      # Kelvin-Mandel components (11, 22, √2·12)
const PLC_W = (1.0, 1.0, sqrt(2.0))
# A normal block without major symmetry and three shear moduli, for the tensors of space.
const PLC_N3 = [160.0 55.0 48.0; 31.0 120.0 42.0; 22.0 17.0 95.0]
const PLC_S3 = (33.0, 29.0, 27.0)

"An angle away from the round values, where `Basis(θ)` would be canonical."
plc_angle(rng) = π * (2 * rand(rng) - 1)

"Rotation of angle `θ`: its columns are the vectors of `Basis(θ)`."
plc_rot(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

"Five constants (C₁₁, C₂₂, C₁₂, C₆₆, C₂₁) of a positive definite normal block."
plc_draw(rng) = (10 + rand(rng), 6 + rand(rng), 2 + rand(rng), 1.5 + rand(rng), 1 + rand(rng))

"Kelvin-Mandel matrices in the material frame, written out."
plc_ortho_km(C11, C22, C12, C66, C21 = C12) = [C11 C12 0; C21 C22 0; 0 0 2C66]
plc_square_km(C11, C12, C44) = [C11 C12 0; C12 C11 0; 0 0 2C44]
plc_km6(N, s) = [N zeros(3, 3); zeros(3, 3) Diagonal(2 .* collect(s))]

plc_index(i, j) = i == j ? i : 3

"Array in the material frame of the Kelvin-Mandel matrix `M`, by the definition of its components."
plc_material_array(M) =
    [M[plc_index(i, j), plc_index(k, l)] / (PLC_W[plc_index(i, j)] * PLC_W[plc_index(k, l)]) for i in 1:2, j in 1:2, k in 1:2, l in 1:2]

"Canonical components of the tensor of material Kelvin-Mandel matrix `M` on the frame of columns `R`."
plc_array(M, R) = ref_rotate(R, plc_material_array(M))
plc_array2(a, R) = R * Diagonal(collect(a)) * transpose(R)

"Kelvin-Mandel matrix of an array of order 4 of the plane, by definition."
plc_km(A) = [PLC_W[I] * PLC_W[J] * A[PLC_PAIRS[I]..., PLC_PAIRS[J]...] for I in 1:3, J in 1:3]

"Kelvin-Mandel rotation of the plane, written out: the canonical matrix is `B M Bᵀ`."
function plc_bond(θ)
    c, s = cos(θ), sin(θ)
    r = sqrt(oftype(c, 2))
    return [c^2 s^2 -r * c * s; s^2 c^2 r * c * s; r * c * s -r * c * s c^2 - s^2]
end

# Canonical components, whatever the basis a generic result is stored in.
plc_arr(t::AbstractTens) = Array(canonical_components(t))
plc_arr(x::AbstractArray) = collect(x)
plc_cmandel(t) = plc_km(plc_arr(t))

plc_inferred_inv(t) = @inferred inv(t)
plc_inferred_material(t) = @inferred mandel_material(t)
plc_inferred_transpose(t) = @inferred transpose(t)
plc_inferred_constants(t) = @inferred cubic_constants(t)
plc_inferred_strain(t) = @inferred plane_strain_stiffness(t)
plc_inferred_blocks(t) = @inferred symmetry_blocks(t)

@testsection "Orthotropic and square classes of the plane" begin
    rng = Xoshiro(PLC_SEED)
    RB = RotatedBasis{2, Float64}
    ℬc = CanonicalBasis{2, Float64}()

    @testsection "construction and types" begin
        ℬ = Basis(plc_angle(rng))
        O = TensOrtho(10.0, 6.0, 3.0, 2.0, ℬ)
        G = TensOrtho(10.0, 6.0, 3.0, 2.0, 2.5, ℬ)
        Q = cubic_from_constants(10.0, 4.0, 2.0, ℬ)
        a = TensOrtho{2}(2.0, 5.0, ℬ)
        @test ℬ isa RB
        @test O isa TensOrtho{4, 2, Float64, 4, RB}
        @test G isa TensOrtho{4, 2, Float64, 5, RB}
        @test Q isa TensCubic{2, Float64, RB}
        @test a isa TensOrtho{2, 2, Float64, 2, RB}
        @test get_data(O) == (10.0, 6.0, 3.0, 2.0)
        @test get_data(G) == (10.0, 6.0, 3.0, 2.0, 2.5)
        @test get_data(a) == (2.0, 5.0)
        # α = C₁₁ + C₁₂, β = C₁₁ − C₁₂, γ = 2C₄₄.
        @test get_data(Q) == (14.0, 6.0, 4.0)
        @test cubic_constants(Q) == (10.0, 4.0, 2.0)
        @test get_data(TensCubic(14.0, 6.0, 4.0, ℬ)) == get_data(Q)
        for t in (O, G, Q)
            @test frame(t) === ℬ
            @test size(t) == (2, 2, 2, 2) && length(t) == 16
            @test get_dim(t) == 2 && eltype(t) === Float64
        end
        @test size(a) == (2, 2) && length(a) == 4 && get_dim(a) == 2 && frame(a) === ℬ
        # Mixed argument types promote.
        @test TensOrtho(10, 6.0, 3, 2, ℬc) isa TensOrtho{4, 2, Float64, 4}
        @test TensOrtho(10, 6, 3, 2, 1.5, ℬc) isa TensOrtho{4, 2, Float64, 5}
        @test TensOrtho{2}(1, 2.0, ℬc) isa TensOrtho{2, 2, Float64, 2}
        @test cubic_from_constants(10, 4, 2, ℬc) isa TensCubic{2, Float64}
        # From the Kelvin-Mandel matrix in the material frame.
        MO, MG, MQ = plc_ortho_km(10.0, 6.0, 3.0, 2.0), plc_ortho_km(10.0, 6.0, 3.0, 2.0, 2.5), plc_square_km(10.0, 4.0, 2.0)
        @test get_data(TensOrtho(MO, ℬ)) == get_data(O)
        @test get_data(TensOrtho(MO, ℬ, Val(4))) == get_data(O)
        @test get_data(TensOrtho(MG, ℬ, Val(5))) == get_data(G)
        @test get_data(TensCubic(MQ, ℬ)) == get_data(Q)
        # Refusals: a wrong number of coefficients, four constants on a frame of space.
        @test_throws ArgumentError TensOrtho{4, 2, Float64}((1.0, 2.0, 3.0), ℬc)
        @test_throws MethodError TensOrtho(1.0, 2.0, 3.0, 4.0, CanonicalBasis{3, Float64}())
    end

    @testsection "components against explicit index sums" begin
        for _ in 1:3
            θ = plc_angle(rng)
            ℬ, R = Basis(θ), plc_rot(θ)
            c = plc_draw(rng)
            cases = (
                (TensOrtho(c[1:4]..., ℬ), plc_ortho_km(c[1:4]...)),
                (TensOrtho(c..., ℬ), plc_ortho_km(c...)),
                (cubic_from_constants(c[1], c[3], c[4], ℬ), plc_square_km(c[1], c[3], c[4])),
            )
            θ′ = plc_angle(rng)
            R′ = plc_rot(θ′)
            for (t, M) in cases
                A = plc_array(M, R)
                At = Array(get_array(t))
                @test At ≈ A
                @test all(t[i, j, k, l] ≈ A[i, j, k, l] for i in 1:2, j in 1:2, k in 1:2, l in 1:2)
                @test canonical_components(t) ≈ A
                @test mandel(t) ≈ plc_km(A)
                @test mandel_material(t) ≈ M
                @test At ≈ permutedims(At, (2, 1, 3, 4)) && At ≈ permutedims(At, (1, 2, 4, 3))
                # In another frame: the components of the reference brought into it.
                @test Array(get_array(express(t, Basis(θ′)))) ≈ ref_rotate_back(R′, A)
                # In its own frame: the material components.
                @test Array(get_array(express(t, ℬ))) ≈ plc_material_array(M)
            end
            a = TensOrtho{2}(c[1], c[2], ℬ)
            Aa = plc_array2((c[1], c[2]), R)
            @test Array(get_array(a)) ≈ Aa
            @test all(a[i, j] ≈ Aa[i, j] for i in 1:2, j in 1:2)
            @test canonical_components(a) ≈ Aa
            @test Array(get_array(express(a, ℬ))) ≈ Diagonal([c[1], c[2]])
        end
    end

    @testsection "algebra on a shared frame" begin
        ℬ = Basis(plc_angle(rng))
        c, d = plc_draw(rng), plc_draw(rng)
        O, G = TensOrtho(c[1:4]..., ℬ), TensOrtho(d..., ℬ)
        Q = cubic_from_constants(c[2], d[3], d[4], ℬ)
        L = iso_stiffness(Val(2), 4.0 + rand(rng), 2.0 + rand(rng))
        MO, MG, MQ = plc_ortho_km(c[1:4]...), plc_ortho_km(d...), plc_square_km(c[2], d[3], d[4])
        # The inverse keeps the shape.
        @test inv(O) isa TensOrtho{4, 2, Float64, 4}
        @test inv(G) isa TensOrtho{4, 2, Float64, 5}
        @test inv(Q) isa TensCubic{2, Float64}
        for (t, M) in ((O, MO), (G, MG), (Q, MQ))
            @test mandel_material(inv(t)) ≈ inv(M)
            @test mandel_material(t^-1) ≈ inv(M)
            @test mandel_material(-t) ≈ -M
            @test mandel_material(2t) ≈ 2M && mandel_material(t * 3) ≈ 3M
            @test mandel_material(t / 4) ≈ M / 4
            @test t ⊙ t ≈ sum(abs2, plc_arr(t))
        end
        # Double contraction: the product of the Kelvin-Mandel matrices.
        @test O ⊡ O isa TensOrtho{4, 2, Float64, 5}
        @test mandel_material(O ⊡ G) ≈ MO * MG
        @test mandel_material(G ⊡ O) ≈ MG * MO
        @test Q ⊡ Q isa TensCubic{2}
        @test mandel_material(Q ⊡ Q) ≈ MQ * MQ
        for (l, r) in ((Q, O), (O, Q), (Q, G), (G, Q), (L, O), (O, L), (L, G))
            @test l ⊡ r isa TensOrtho{4, 2}
            @test mandel(l ⊡ r) ≈ mandel(l) * mandel(r)
        end
        for (l, r) in ((L, Q), (Q, L))
            @test l ⊡ r isa TensCubic{2}
            @test mandel(l ⊡ r) ≈ mandel(l) * mandel(r)
        end
        # Sums: the smallest shape that holds the result.
        @test O + O isa TensOrtho{4, 2, Float64, 4}
        @test G - G isa TensOrtho{4, 2, Float64, 5}
        for (l, r) in ((O, G), (G, O))
            @test l + r isa TensOrtho{4, 2, Float64, 5}
            @test mandel(l + r) ≈ mandel(l) + mandel(r)
            @test mandel(l - r) ≈ mandel(l) - mandel(r)
        end
        for (l, r) in ((Q, O), (O, Q), (Q, G), (L, O), (O, L), (L, G))
            @test l + r isa TensOrtho{4, 2}
            @test l - r isa TensOrtho{4, 2}
            @test mandel(l + r) ≈ mandel(l) + mandel(r)
            @test mandel(l - r) ≈ mandel(l) - mandel(r)
        end
        for (l, r) in ((L, Q), (Q, L))
            @test l + r isa TensCubic{2} && l - r isa TensCubic{2}
            @test mandel(l + r) ≈ mandel(l) + mandel(r)
            @test mandel(l - r) ≈ mandel(l) - mandel(r)
        end
        # Transpose and major-symmetric part.
        @test transpose(O) === O
        Gt = transpose(G)
        @test get_data(Gt) == (d[1], d[2], d[5], d[4], d[3])
        @test mandel_material(Gt) ≈ transpose(MG)
        @test plc_arr(Gt) ≈ permutedims(plc_arr(G), (3, 4, 1, 2))
        S = majorsymmetric(G)
        @test S isa TensOrtho{4, 2, Float64, 4}
        @test collect(get_data(S)) ≈ [d[1], d[2], (d[3] + d[5]) / 2, d[4]]
        @test majorsymmetric(O) === O
        @test issymmetric(O) && ismajorsymmetric(O) && isminorsymmetric(O)
        @test !issymmetric(G) && !ismajorsymmetric(G) && isminorsymmetric(G)
        @test issymmetric(TensOrtho(1.0, 2.0, 0.5, 0.7, 0.5, ℬ))
        @test issymmetric(Q) && ismajorsymmetric(Q)
        # zero and one of the square class stay in the class and on the frame.
        @test zero(Q) isa TensCubic{2} && all(iszero, get_data(zero(Q))) && frame(zero(Q)) === ℬ
        @test get_data(one(Q)) == (1.0, 1.0, 1.0)
        @test mandel(one(Q)) ≈ Matrix(1.0I, 3, 3)
        @test all(iszero, plc_arr(zero(O)))
        # Class predicates.
        @test symmetry(O) == symmetry(G) == :ORTHO
        @test symmetry(Q) == :CUBIC
        @test isorthotropic(O) && isorthotropic(G) && !isisotropic(O)
        @test iscubic(Q) && !isorthotropic(Q) && !isisotropic(Q)
    end

    @testsection "operands on different frames" begin
        c, d = plc_draw(rng), plc_draw(rng)
        θ₁, θ₂ = plc_angle(rng), plc_angle(rng)
        A, B = TensOrtho(c[1:4]..., Basis(θ₁)), TensOrtho(d..., Basis(θ₂))
        @test !(A + B isa TensOrtho)
        @test plc_cmandel(A + B) ≈ mandel(A) + mandel(B)
        @test plc_cmandel(A - B) ≈ mandel(A) - mandel(B)
        @test plc_cmandel(A ⊡ B) ≈ mandel(A) * mandel(B)
        @test plc_cmandel(B ⊡ A) ≈ mandel(B) * mandel(A)
        # A quarter turn exchanges the axes: the rectangle is the same, the
        # constants are relabeled and the sum stays in the class.
        B′ = TensOrtho(d[1:4]..., Basis(θ₁ + π / 2))
        for A′ in (TensOrtho(c[1:4]..., Basis(θ₁)), TensOrtho(c..., Basis(θ₁)))
            @test A′ + B′ isa typeof(A′)
            @test mandel(A′ + B′) ≈ mandel(A′) + mandel(B′)
            @test canonical_components(A′ + B′) ≈ canonical_components(A′) + canonical_components(B′)
        end
        Q₁ = cubic_from_constants(c[1], c[3], c[4], Basis(θ₁))
        Q₂ = cubic_from_constants(d[1], d[3], d[4], Basis(θ₁ + π / 2))
        @test Q₁ + Q₂ isa TensCubic{2} && Q₁ ⊡ Q₂ isa TensCubic{2}
        @test mandel(Q₁ + Q₂) ≈ mandel(Q₁) + mandel(Q₂)
        @test mandel(Q₁ ⊡ Q₂) ≈ mandel(Q₁) * mandel(Q₂)
    end

    @testsection "order 2" begin
        θ = plc_angle(rng)
        ℬ, R = Basis(θ), plc_rot(θ)
        a, b = TensOrtho{2}(2.0 + rand(rng), 5.0 + rand(rng), ℬ), TensOrtho{2}(1.0 + rand(rng), 3.0 + rand(rng), ℬ)
        Aa, Ab = plc_array2(get_data(a), R), plc_array2(get_data(b), R)
        @test tr(a) ≈ tr(Aa)
        @test inv(a) isa TensOrtho{2, 2, Float64, 2}
        @test plc_arr(inv(a)) ≈ inv(Aa)
        @test plc_arr(a ⋅ b) ≈ Aa * Ab
        λ = TensISO{2}(1.5)
        @test plc_arr(a ⋅ λ) ≈ 1.5 * Aa && plc_arr(λ ⋅ a) ≈ 1.5 * Aa
        v = [0.3, -1.2]
        @test plc_arr(a ⋅ v) ≈ Aa * v
        @test plc_arr(a ⋅ Tens(v)) ≈ Aa * v
        # With order-4 tensors of the same frame.
        c = plc_draw(rng)
        for t in (TensOrtho(c[1:4]..., ℬ), TensOrtho(c..., ℬ), cubic_from_constants(c[1], c[3], c[4], ℬ), iso_stiffness(Val(2), 3.0, 1.5))
            T4 = plc_arr(t)
            @test plc_arr(t ⊡ a) ≈ ref_dc(T4, Aa)
            @test plc_arr(a ⊡ t) ≈ ref_dc(Aa, T4)
        end
        # Tensor products: Σₘₙ aₘbₙ Pₘ⊗Pₙ. The self product of one argument is
        # major-symmetric, four constants; the product of two operands takes the
        # general shape, whatever they are.
        @test Tensors.otimes(a) isa TensOrtho{4, 2, Float64, 4}
        @test plc_arr(Tensors.otimes(a)) ≈ ref_outer(Aa, Aa)
        @test a ⊗ a isa TensOrtho{4, 2, Float64, 5}
        @test a ⊗ b isa TensOrtho{4, 2, Float64, 5}
        @test plc_arr(a ⊗ b) ≈ ref_outer(Aa, Ab)
        @test plc_arr(a ⊗ a) ≈ ref_outer(Aa, Aa)
        # Frames a quarter turn apart describe the same rectangle: the sum stays
        # orthotropic, the coefficients of the second operand exchanged.
        a⊥ = TensOrtho{2}(0.7, -1.3, Basis(θ + π / 2))
        @test a + a⊥ isa TensOrtho{2, 2}
        @test plc_arr(a + a⊥) ≈ Aa + plc_arr(a⊥)
        @test collect(get_data(a + a⊥)) ≈ collect(get_data(a)) .+ [-1.3, 0.7]
        # Sums with an isotropic tensor, and on two unrelated frames.
        @test λ + a isa TensOrtho{2, 2}
        @test plc_arr(λ + a) ≈ Aa + 1.5I && plc_arr(a - λ) ≈ Aa - 1.5I
        θ′ = plc_angle(rng)
        e = TensOrtho{2}(1.0, 3.0, Basis(θ′))
        s = a + e
        @test !(s isa TensOrtho)
        @test plc_arr(s) ≈ Aa + plc_array2((1.0, 3.0), plc_rot(θ′))
        @test plc_arr(s) ≈ transpose(plc_arr(s))
        @test isisotropic(TensOrtho{2}(2.0, 2.0, ℬ)) && !isisotropic(a)
    end

    @testsection "symmetry blocks" begin
        ℬ = Basis(plc_angle(rng))
        c = plc_draw(rng)
        O, G = TensOrtho(c[1:4]..., ℬ), TensOrtho(c..., ℬ)
        Q = cubic_from_constants(c[1], c[3], c[4], ℬ)
        @test symmetry_layout(O) == ((SVector(1, 2),), (SVector(3),))
        @test symmetry_layout(Q) == ((SVector(1),), (SVector(2),), (SVector(3),))
        @test symmetry_blocks(G)[1] ≈ [c[1] c[3]; c[5] c[2]]
        @test symmetry_blocks(O)[2][1] ≈ 2c[4]
        @test collect(only.(symmetry_blocks(Q))) ≈ [c[1] + c[3], c[1] - c[3], 2c[4]]
        for t in (O, G, Q)
            b, B = symmetry_blocks(t), symmetry_basis(t)
            @test B' * B ≈ I
            # Qᵀ M Q is block diagonal, each block on its copies.
            D = zeros(3, 3)
            for (bk, copies) in zip(b, symmetry_layout(t)), k in copies
                D[k, k] .= bk
            end
            @test B' * mandel(t) * B ≈ D
            @test symmetry_matrix(t, b) ≈ mandel(t)
            @test plc_arr(from_symmetry_blocks(t, b)) ≈ plc_arr(t)
            @test all(map((x, y) -> x ≈ y, symmetry_blocks(t, mandel(t)), b))
        end
        a = TensOrtho{2}(c[1], c[2], ℬ)
        b, B = symmetry_blocks(a), symmetry_basis(a)
        @test symmetry_layout(a) == ((SVector(1),), (SVector(2),))
        @test collect(only.(b)) == [c[1], c[2]]
        @test B' * B ≈ I && B' * plc_arr(a) * B ≈ Diagonal([c[1], c[2]])
        @test symmetry_matrix(a, b) ≈ plc_arr(a)
        @test plc_arr(from_symmetry_blocks(a, b)) ≈ plc_arr(a)
    end

    @testsection "promotions and the rotational average" begin
        k, μ = 4.0 + rand(rng), 1.5 + rand(rng)
        L = iso_stiffness(Val(2), k, μ)
        @test collect(get_data(L)) ≈ [2k, 2μ]
        ℬ = Basis(plc_angle(rng))
        # α𝕁 + β𝕂 on any frame: C₁₁ = C₂₂ = k + μ, C₁₂ = k − μ, C₆₆ = μ.
        Lo = iso_to_ortho(L, ℬ)
        @test Lo isa TensOrtho{4, 2, Float64, 4}
        @test collect(get_data(Lo)) ≈ [k + μ, k + μ, k - μ, μ]
        @test mandel(Lo) ≈ mandel(L)
        Lc = iso_to_cubic(L, ℬ)
        @test Lc isa TensCubic{2, Float64}
        @test mandel(Lc) ≈ mandel(L)
        @test abs(cubic_anisotropy(Lc)) < 1.0e-14
        λ = TensISO{2}(1.7)
        @test iso_to_ortho(λ, ℬ) isa TensOrtho{2, 2, Float64, 2}
        @test plc_arr(iso_to_ortho(λ, ℬ)) ≈ 1.7 * Matrix(1.0I, 2, 2)
        C11, C12, C44 = 10.0 + rand(rng), 3.0 + rand(rng), 2.0 + rand(rng)
        Q = cubic_from_constants(C11, C12, C44, ℬ)
        Qo = cubic_to_ortho(Q)
        @test Qo isa TensOrtho{4, 2, Float64, 4} && frame(Qo) === ℬ
        @test collect(get_data(Qo)) ≈ [C11, C11, C12, C44]
        @test mandel(Qo) ≈ mandel(Q)
        @test cubic_anisotropy(Q) ≈ (C11 - C12 - 2C44) / C11
        # The average over the rotations of the plane: the entries are
        # trigonometric polynomials of degree 4 in θ, of period π/2 for the
        # square, so a uniform grid of 64 angles integrates them exactly.
        α, β, γ = get_data(Q)
        Qa = isotropify(Q)
        @test Qa isa TensISO{4, 2}
        @test collect(get_data(Qa)) ≈ [α, (β + γ) / 2]
        n = 64
        avg = sum(mandel(cubic_from_constants(C11, C12, C44, Basis((m - 0.5) * π / (2n)))) for m in 1:n) / n
        @test avg ≈ mandel(Qa)
    end

    @testsection "reductions from space" begin
        ℬ3 = random_basis(rng)
        Ns = (PLC_N3 + transpose(PLC_N3)) / 2
        ℂ9 = TensOrtho(Ns[1, 1], Ns[2, 2], Ns[3, 3], Ns[1, 2], Ns[1, 3], Ns[2, 3], PLC_S3..., ℬ3)
        N = PLC_N3
        ℂ12 = TensOrtho(N[1, 1], N[2, 2], N[3, 3], N[1, 2], N[1, 3], N[2, 3], PLC_S3..., N[2, 1], N[3, 1], N[3, 2], ℬ3)
        CB = CanonicalBasis{2, Float64}
        for (ℂ, Nb, len) in ((ℂ9, Ns, 4), (ℂ12, N, 5))
            # Material frame of space, components 11, 22, 33, 23, 13, 12.
            M6 = plc_km6(Nb, PLC_S3)
            for k in 1:3
                # The plane normal to e_k, of axes (e_i, e_j) in the cyclic order,
                # and its shear component 3 + k.
                i, j = mod1(k + 1, 3), mod1(k + 2, 3)
                idx = [i, j, 3 + k]
                ε = plane_strain_stiffness(ℂ, k)
                σ = plane_stress_stiffness(ℂ, k)
                @test ε isa TensOrtho{4, 2, Float64, len, CB}
                @test σ isa TensOrtho{4, 2, Float64, len, CB}
                @test mandel(ε) ≈ M6[idx, idx]
                @test mandel(σ) ≈ inv(inv(M6)[idx, idx])
            end
            @test get_data(plane_strain_stiffness(ℂ)) == get_data(plane_strain_stiffness(ℂ, 3))
            @test_throws ArgumentError plane_strain_stiffness(ℂ, 4)
            @test_throws ArgumentError plane_stress_stiffness(ℂ, 0)
        end
        idx = [1, 2, 6]
        C11, C12, C44 = 170.0, 120.0, 75.0
        Cc = cubic_from_constants(C11, C12, C44, ℬ3)
        M6c = plc_km6([C11 C12 C12; C12 C11 C12; C12 C12 C11], (C44, C44, C44))
        @test plane_strain_stiffness(Cc) isa TensCubic{2, Float64}
        @test plane_stress_stiffness(Cc) isa TensCubic{2, Float64}
        @test mandel(plane_strain_stiffness(Cc)) ≈ M6c[idx, idx]
        @test mandel(plane_stress_stiffness(Cc)) ≈ inv(inv(M6c)[idx, idx])
        # Transversely isotropic about e₁, in the plane (e₁, e₂) that contains the axis.
        ℂti = ti_from_constants(10.0, 4.0, 3.0, 20.0, 2.5, [1.0, 0.0, 0.0]; form = :hill)
        M6t = Matrix(mandel(ℂti))
        @test plane_strain_stiffness(ℂti, Val(:axial)) isa TensOrtho{4, 2, Float64, 4}
        @test mandel(plane_strain_stiffness(ℂti, Val(:axial))) ≈ M6t[idx, idx]
        @test mandel(plane_stress_stiffness(ℂti, Val(:axial))) ≈ inv(inv(M6t)[idx, idx])
        # In the plane normal to the axis the result stays isotropic.
        @test plane_strain_stiffness(ℂti) isa TensISO{4, 2}
        @test plane_stress_stiffness(ℂti) isa TensISO{4, 2}
        @test mandel(plane_strain_stiffness(ℂti)) ≈ M6t[[2, 3, 4], [2, 3, 4]]
        @test plane_strain_stiffness(TensISO{3}(30.0, 8.0)) isa TensISO{4, 2}
    end

    @testsection "random draws" begin
        @test rand(rng, TensOrtho{4, 2}) isa TensOrtho{4, 2, Float64, 4, RB}
        @test rand(rng, TensOrtho{4, 2, Float64, 5}) isa TensOrtho{4, 2, Float64, 5, RB}
        @test rand(rng, TensOrtho{2, 2}) isa TensOrtho{2, 2, Float64, 2, RB}
        @test rand(rng, TensCubic{2}) isa TensCubic{2, Float64, RB}
        for Tt in (TensOrtho{4, 2}, TensCubic{2})
            t = rand(rng, Tt; definite = true)
            @test all(>(0), eigvals(Symmetric(Matrix(mandel(t)))))
        end
        @test all(>(0), get_data(rand(rng, TensOrtho{2, 2}; definite = true)))
        ℬ = Basis(0.3)
        for Tt in (TensOrtho{4, 2}, TensOrtho{4, 2, Float64, 5}, TensOrtho{2, 2}, TensCubic{2})
            @test frame(rand(rng, Tt; frame = ℬ)) === ℬ
        end
        @test rand(rng, TensOrtho{4, 2, Float32, 4}) isa TensOrtho{4, 2, Float32, 4}
        @test rand(rng, TensOrtho{2, 2, Float32}) isa TensOrtho{2, 2, Float32}
        @test rand(rng, TensCubic{2, Float32}) isa TensCubic{2, Float32}
        @test_throws ArgumentError rand(rng, TensOrtho{4, 2, Float64, 5}; definite = true)
        # Without a dimension the draws stay in space.
        @test rand(rng, TensOrtho{4}) isa TensOrtho{4, 3}
        @test rand(rng, TensCubic) isa TensCubic{3}
        @test rand(rng, TensOrtho{2}) isa TensOrtho{2, 3}
    end

    @testsection "inference" begin
        ℬ = Basis(plc_angle(rng))
        O, G = TensOrtho(10.0, 6.0, 3.0, 2.0, ℬ), TensOrtho(10.0, 6.0, 3.0, 2.0, 2.5, ℬ)
        Q = cubic_from_constants(10.0, 4.0, 2.0, ℬ)
        for t in (O, G, Q)
            @test plc_inferred_inv(t) isa AbstractTens{4, 2}
            @test plc_inferred_material(t) isa SMatrix{3, 3, Float64}
            @test plc_inferred_blocks(t) isa Tuple
        end
        @test plc_inferred_transpose(O) === O
        @test plc_inferred_transpose(G) isa TensOrtho{4, 2, Float64, 5}
        @test plc_inferred_constants(Q) == (10.0, 4.0, 2.0)
        @test plc_inferred_strain(cubic_from_constants(10.0, 4.0, 2.0, random_basis(rng))) isa TensCubic{2}
        @test plc_inferred_inv(TensOrtho(ForwardDiff.Dual(10.0, 1.0), 6.0, 3.0, 2.0, ℬ)) isa TensOrtho{4, 2}
    end

    @testsection "derivatives" begin
        θ = plc_angle(rng)
        ℬ = Basis(θ)
        c11, c22, c12, c66 = 10.0, 6.0, 3.0, 2.0
        # Independent route: the canonical Kelvin-Mandel matrix as B(θ) M Bᵀ(θ).
        km(x, ϑ) = plc_bond(ϑ) * plc_ortho_km(c11 * x, c22, c12, c66) * transpose(plc_bond(ϑ))
        x₀ = 1.3
        dx = ForwardDiff.derivative(x -> Matrix(mandel(inv(TensOrtho(c11 * x, c22, c12, c66, ℬ)))), x₀)
        @test dx ≈ ForwardDiff.derivative(x -> inv(km(x, θ)), x₀)
        # d(M⁻¹) = −M⁻¹ dM M⁻¹, with dM the canonical image of c₁₁ P₁⊗P₁.
        B = plc_bond(θ)
        Mi = inv(km(x₀, θ))
        @test dx ≈ -Mi * (B * plc_ortho_km(c11, 0.0, 0.0, 0.0) * transpose(B)) * Mi
        # With respect to the angle of the frame.
        dθ = ForwardDiff.derivative(ϑ -> Matrix(mandel(TensOrtho(c11, c22, c12, c66, Basis(ϑ)))), θ)
        @test dθ ≈ ForwardDiff.derivative(ϑ -> km(1.0, ϑ), θ)
        # Square class: the third row of B(θ) is (u, −u, cos 2θ), u = sin 2θ/√2, so
        # the entry (3, 3) is 2u²(C₁₁ − C₁₂) + 2C₄₄ cos² 2θ = γ + (β − γ) sin² 2θ,
        # and its second derivative 8(β − γ) cos 4θ.
        C11, C12, C44 = 10.0, 4.0, 2.0
        β, γ = C11 - C12, 2C44
        m33(ϑ) = mandel(cubic_from_constants(C11, C12, C44, Basis(ϑ)))[3, 3]
        @test m33(θ) ≈ γ + (β - γ) * sin(2θ)^2
        @test ForwardDiff.derivative(m33, θ) ≈ 2 * (β - γ) * sin(4θ)
        @test ForwardDiff.derivative(ϑ -> ForwardDiff.derivative(m33, ϑ), θ) ≈ 8 * (β - γ) * cos(4θ)
    end

    @testsection "symbols" begin
        C11, C22, C12, C66 = symbols("C11 C22 C12 C66", positive = true)
        θ = symbols("θ", real = true)
        O = TensOrtho(C11, C22, C12, C66, Basis(θ))
        @test all(iszero, simplify.(Matrix(mandel_material(O)) - [C11 C12 0; C12 C22 0; 0 0 2C66]))
        # 2C′₁₂₁₂ in the canonical frame: 2C₆₆ cos² 2θ + (C₁₁ + C₂₂ − 2C₁₂)/2 sin² 2θ.
        @test simplify(mandel(O)[3, 3] - (2C66 * cos(2θ)^2 + (C11 + C22 - 2C12) / 2 * sin(2θ)^2)) == 0
        Oc = TensOrtho(C11, C22, C12, C66, CanonicalBasis{2, typeof(C11)}())
        @test inv(Oc) isa TensOrtho{4, 2, typeof(C11), 4}
        @test all(iszero, simplify.(Matrix(mandel(inv(Oc))) * Matrix(mandel(Oc)) - I))
        # Symbolics: the inverse in closed form, evaluated against the numbers.
        Sy = TensND.Symbolics
        TensND.Symbolics.@variables s11 s22 s12 s66
        Os = inv(TensOrtho(s11, s22, s12, s66, CanonicalBasis{2, Sy.Num}()))
        @test Os isa TensOrtho{4, 2, Sy.Num, 4}
        vals = Dict(s11 => 10.0, s22 => 6.0, s12 => 3.0, s66 => 2.0)
        num = [Sy.value(Sy.substitute(x, vals; fold = Val(true))) for x in get_data(Os)]
        @test num ≈ collect(get_data(inv(TensOrtho(10.0, 6.0, 3.0, 2.0, ℬc))))
    end
end
