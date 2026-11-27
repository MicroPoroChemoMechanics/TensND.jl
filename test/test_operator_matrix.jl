# =============================================================================
#  test_operator_matrix.jl — a tensor as the matrix of the linear map it is
#  (`operator_matrix`, `from_operator_matrix`), its definiteness, its departure
#  from major symmetry, and the isotropic average of a Kelvin-Mandel matrix in
#  any dimension (`iso_average_mandel`).
#
#  References, independent of the code under test:
#  - the matrix maps the Kelvin-Mandel vector of a random symmetric ε onto that
#    of `A : ε`, the double contraction written as an explicit index sum
#    (`reference_contractions.jl`), the vectors written from the weights 1 and √2 of the
#    definition, the minor-symmetric part built by averaging index permutations;
#  - a tensor in a rotated basis is compared with its components rotated by
#    explicit index sums;
#  - composition: the matrix of `A ⊡ B` is the product of the matrices;
#  - the isotropic average is compared with a quadrature over the rotations
#    and reflections of the plane, with `isotropify` (a Frobenius projection on
#    the component arrays) in dimensions 2, 3 and 4, and with the closed form of
#    space written out;
#  - definiteness on tensors whose eigenvalues are known (isotropic ones,
#    `S + W` with S symmetric positive definite and W skew);
#  - the asymmetry of `S + x W` is `x ‖W‖/‖S‖`, whose derivative is checked by
#    ForwardDiff; the matrix of an isotropic tensor with symbolic moduli is
#    compared with the one written out.
# =============================================================================

using StaticArrays: SMatrix

const OPM_SEED = 20261126

"The Kelvin-Mandel vector of a symmetric matrix `ε`, from the weights 1 and √2 of the definition."
opm_vec(ε) = [(i == j ? 1.0 : sqrt(2.0)) * ε[i, j] for (i, j) in mandel_pairs(size(ε, 1))]

"The Kelvin-Mandel matrix of a minor-symmetric array, entry by entry from the definition."
function opm_km(A)
    P = mandel_pairs(size(A, 1))
    w(p) = p[1] == p[2] ? 1.0 : sqrt(2.0)
    return [w(p) * w(q) * A[p..., q...] for p in P, q in P]
end

"The minor-symmetric part of an array of order 4, by averaging its index permutations."
function opm_minor(A)
    A = (A + permutedims(A, (2, 1, 3, 4))) / 2
    return (A + permutedims(A, (1, 2, 4, 3))) / 2
end

opm_sym(rng, d) = (e = randn(rng, d, d); e + transpose(e))

"The canonical components of the tensor of components `A` in the basis of columns `R`."
opm_rotate(A, R) = ref_rotate(R, A)

opm_arr(t) = Array(canonical_components(t))

"Whether `M` maps the Kelvin-Mandel vector of a random symmetric ε onto that of `A : ε`, `A` minor-symmetric."
function opm_maps(M, A, rng)
    d = size(A, 1)
    return all(1:3) do _
        ε = opm_sym(rng, d)
        opm_vec(ref_dc(A, ε)) ≈ M * opm_vec(ε)
    end
end

"A symmetric positive definite matrix and a skew one, of size `n`."
function opm_spd_skew(rng, n)
    X, Y = randn(rng, n, n), randn(rng, n, n)
    return X * transpose(X) + n * I, Y - transpose(Y)
end

opm_inferred(t) = @inferred operator_matrix(t)
opm_inferred_basis(t, ℬ) = @inferred operator_matrix(t, ℬ)

@testsection "Operator matrix" begin
    rng = Xoshiro(OPM_SEED)
    ℬ2 = Basis(2π * rand(rng))
    ℬ3 = random_basis(rng)
    n̂ = random_axis(rng)

    @testsection "order 4: the map on symmetric tensors, any basis" begin
        for d in 2:4, _ in 1:2
            n = mandel_length(d)
            A = randn(rng, d, d, d, d)
            As = opm_minor(A)
            M = operator_matrix(Tens(A))
            @test opm_maps(M, As, rng)
            @test M ≈ operator_matrix(Tens(As); symmetrize = false)
            @test_throws ArgumentError operator_matrix(Tens(A); symmetrize = false)

            ℬ = random_basis(rng, Val(d))
            R = Matrix(vecbasis(ℬ))
            t = Tens(As, ℬ)
            Mr = operator_matrix(t)
            @test opm_maps(Mr, opm_rotate(As, R), rng)
            @test operator_matrix(Tens(A, ℬ)) ≈ Mr
            # In a basis: the matrix of the components in that basis.
            t0 = Tens(opm_rotate(As, R))
            @test operator_matrix(t, ℬ) ≈ opm_km(As)
            @test operator_matrix(t0, ℬ) ≈ opm_km(As)
            @test operator_matrix(t0, ℬ) ≈ opm_km(Array(get_array(express(t0, ℬ))))
            if d ≤ 3
                @test opm_inferred(t) isa SMatrix{n, n, Float64}
                @test opm_inferred(Tens(A)) isa SMatrix{n, n, Float64}
                @test opm_inferred_basis(t0, ℬ) isa SMatrix{n, n, Float64}
            else
                @test M isa Matrix{Float64}
                @test operator_matrix(t0, ℬ) isa Matrix{Float64}
            end
        end
    end

    @testsection "order 4: structured tensors" begin
        structured = (
            TensISO{3}(3.0, 2.0), TensISO{2}(2.0, 1.5), TensISO{4}(4.0, 1.0),
            TensOrtho(10.0, 8.0, 9.0, 3.0, 2.0, 4.0, 2.5, 3.0, 1.5, ℬ3), TensCubic(18.0, 6.0, 4.0, ℬ3),
            TensOrtho(10.0, 6.0, 3.0, 2.0, ℬ2), TensOrtho(10.0, 6.0, 3.0, 2.0, 2.5, ℬ2),
            cubic_from_constants(10.0, 4.0, 2.0, ℬ2),
            ti_from_constants(19.0, 5.0, 11.0, 12.0, 4.0, n̂; form = :hill),
            transverse_isotropify(Tens(opm_minor(randn(rng, 3, 3, 3, 3))), n̂),
        )
        for t in structured
            d = get_dim(t)
            n = mandel_length(d)
            A = opm_arr(t)
            M = operator_matrix(t)
            @test opm_maps(M, A, rng)
            @test M ≈ opm_km(A)
            @test M ≈ operator_matrix(t; symmetrize = false)
            ℬ = d == 2 ? ℬ2 : random_basis(rng, Val(d))
            @test operator_matrix(t, ℬ) ≈ opm_km(opm_rotate(A, transpose(Matrix(vecbasis(ℬ)))))
            @test opm_arr(from_operator_matrix(M, Val(4))) ≈ A
            if d ≤ 3
                @test opm_inferred(t) isa SMatrix{n, n, Float64}
            else
                @test M isa Matrix{Float64}
            end
        end
    end

    @testsection "order 4: composition and the reciprocal" begin
        for d in 2:4
            ℬa, ℬb = random_basis(rng, Val(d)), random_basis(rng, Val(d))
            A = Tens(opm_minor(randn(rng, d, d, d, d)), ℬa)
            B = Tens(opm_minor(randn(rng, d, d, d, d)), ℬb)
            @test operator_matrix(A ⊡ B) ≈ operator_matrix(A) * operator_matrix(B)
            for t in (A, B)
                @test opm_arr(from_operator_matrix(operator_matrix(t), Val(4))) ≈ opm_arr(t)
                @test opm_arr(from_operator_matrix(operator_matrix(t, ℬa), Val(4), ℬa)) ≈ opm_arr(t)
            end
        end
        @test_throws DimensionMismatch from_operator_matrix(zeros(4, 4), Val(4))
        @test_throws DimensionMismatch from_operator_matrix(zeros(3, 6), Val(4))
        @test_throws DimensionMismatch from_operator_matrix(zeros(2, 3), Val(2))
    end

    @testsection "order 2: the map on vectors" begin
        for d in 2:4
            a, v = randn(rng, d, d), randn(rng, d)
            @test operator_matrix(Tens(a)) * v ≈ a * v
            ℬ = random_basis(rng, Val(d))
            R = Matrix(vecbasis(ℬ))
            t = Tens(a, ℬ)
            @test operator_matrix(t) * v ≈ R * a * transpose(R) * v
            @test operator_matrix(t, ℬ) ≈ a
            @test opm_arr(from_operator_matrix(operator_matrix(t), Val(2))) ≈ R * a * transpose(R)
            @test opm_arr(from_operator_matrix(a, Val(2), ℬ)) ≈ R * a * transpose(R)
            if d ≤ 3
                @test opm_inferred(t) isa SMatrix{d, d, Float64}
            else
                @test operator_matrix(t) isa Matrix{Float64}
            end
        end
        @test operator_matrix(TensISO{3}(2.0)) ≈ 2.0 * Matrix(I, 3, 3)
    end

    @testsection "derivatives and symbols" begin
        A = Tens(opm_minor(randn(rng, 3, 3, 3, 3)), ℬ3)
        MA = operator_matrix(A)
        @test ForwardDiff.derivative(x -> operator_matrix(x * A), 0.7) ≈ MA
        @test ForwardDiff.derivative(y -> ForwardDiff.derivative(x -> operator_matrix(x * y * A), 0.7), 1.1) ≈ MA
        @test ForwardDiff.derivative(x -> operator_matrix(TensISO{3}(x, 2.0)), 3.0) ≈ operator_matrix(TensISO{3}(1.0, 0.0))
        @test ForwardDiff.derivative(x -> operator_matrix(TensOrtho(x, 6.0, 3.0, 2.0, ℬ2)), 10.0) ≈
            operator_matrix(TensOrtho(1.0, 0.0, 0.0, 0.0, ℬ2))
        a = Tens(randn(rng, 3, 3), ℬ3)
        @test ForwardDiff.derivative(x -> operator_matrix(x * a), 0.3) ≈ operator_matrix(a)

        k, μ = symbols("k μ", positive = true)
        C = iso_stiffness(Val(2), k, μ)
        expected = [k + μ k - μ 0; k - μ k + μ 0; 0 0 2μ]
        M = operator_matrix(C)
        @test M isa SMatrix{3, 3}
        @test all(iszero, simplify.(M - expected))
        @test all(iszero, simplify.(operator_matrix(C, ℬ2) - expected))
        @test all(iszero, simplify.(operator_matrix(Tens(Array(get_array(C)))) - expected))
    end

    @testsection "definiteness" begin
        C = TensISO{3}(3.0, 2.0)
        @test definiteness(C) === :positive
        @test definiteness(-C) === :negative
        @test definiteness(TensISO{3}(3.0, -2.0)) === :indefinite
        @test definiteness(TensISO{3}(0.0, 0.0)) === :null
        @test definiteness(TensISO{3}(3.0, 0.0)) === :positive_semidefinite
        @test definiteness(TensISO{3}(0.0, -2.0)) === :negative_semidefinite
        @test definiteness(TensISO{3}(3.0, 1.0e-14)) === :positive_semidefinite
        @test definiteness(TensISO{3}(3.0, 1.0e-14); rtol = 1.0e-16) === :positive
        @test definiteness(TensISO{4}(4.0, 1.0)) === :positive

        # Not major-symmetric: the sign is that of the symmetric part.
        S, W = opm_spd_skew(rng, 6)
        @test definiteness(from_operator_matrix(SMatrix{6, 6}(S + 10W), Val(4))) === :positive
        @test definiteness(from_operator_matrix(SMatrix{6, 6}(-S + 10W), Val(4), ℬ3)) === :negative
        @test definiteness(Tens([1.0 5.0; -5.0 1.0])) === :positive
        @test definiteness(TensISO{3}(-1.0)) === :negative

        # Duals are read on their values.
        D = ForwardDiff.Dual{Nothing}(1.0, 1.0)
        @test definiteness(D * C) === :positive
        @test definiteness(-D * C) === :negative
        DD = ForwardDiff.Dual{Missing}(ForwardDiff.Dual{Nothing}(1.0, 1.0), ForwardDiff.Dual{Nothing}(0.5, 0.0))
        @test definiteness(DD * TensISO{3}(3.0, -2.0)) === :indefinite

        k, μ = symbols("k μ", positive = true)
        @test_throws ArgumentError definiteness(iso_stiffness(k, μ))
    end

    @testsection "major asymmetry" begin
        @test major_asymmetry(TensISO{3}(3.0, 2.0)) == 0
        @test major_asymmetry(TensISO{3}(0.0, 0.0)) == 0
        A = opm_minor(randn(rng, 3, 3, 3, 3))
        @test major_asymmetry(Tens((A + permutedims(A, (3, 4, 1, 2))) / 2, ℬ3)) < 1.0e-14

        S, W = opm_spd_skew(rng, 6)
        @test major_asymmetry(from_operator_matrix(SMatrix{6, 6}(S + W), Val(4))) ≈ norm(W) / norm(S)
        asym(x) = major_asymmetry(from_operator_matrix(SMatrix{6, 6}(S + x * W), Val(4)))
        @test ForwardDiff.derivative(asym, 0.8) ≈ norm(W) / norm(S)
        @test ForwardDiff.derivative(asym, 2.5) ≈ norm(W) / norm(S)

        for t in (TensOrtho(10.0, 6.0, 3.0, 2.0, 2.5, ℬ2), Tens(opm_minor(randn(rng, 4, 4, 4, 4)), random_basis(rng, Val(4))))
            K = opm_km(opm_arr(t))
            @test major_asymmetry(t) ≈ norm(K - transpose(K)) / norm(K + transpose(K))
        end
    end

    @testsection "isotropic average of a Kelvin-Mandel matrix, any dimension" begin
        # Space: the closed form, its sums in the same order.
        M6 = randn(rng, 6, 6)
        α3 = (M6[1, 1] + M6[1, 2] + M6[1, 3] + M6[2, 1] + M6[2, 2] + M6[2, 3] + M6[3, 1] + M6[3, 2] + M6[3, 3]) / 3
        β3 = (M6[1, 1] + M6[2, 2] + M6[3, 3] + M6[4, 4] + M6[5, 5] + M6[6, 6] - α3) / 5
        @test iso_average_mandel(M6) === (α3, β3)
        @test iso_average_mandel(SMatrix{6, 6}(M6)) === (α3, β3)

        for d in 2:4
            n = mandel_length(d)
            A = opm_minor(randn(rng, d, d, d, d))
            M = opm_km(A)
            α, β = iso_average_mandel(M)
            @test [α, β] ≈ collect(get_data(isotropify(Tens(A))))
            @test collect(iso_average_mandel(SMatrix{n, n}(M))) ≈ [α, β]
        end

        # The plane: rotations and reflections, by a quadrature exact on the
        # trigonometric polynomials of the entries.
        𝐢 = [1.0, 1.0, 0.0]
        J = 𝐢 * transpose(𝐢) / 2
        M = opm_km(opm_minor(randn(rng, 2, 2, 2, 2)))
        P = mandel_rotation(SMatrix{2, 2}(1.0, 0.0, 0.0, -1.0))
        avg = sum(1:64) do m
            Q = mandel_rotation(2π * m / 64)
            (Q * M * transpose(Q) + Q * P * M * P * transpose(Q)) / 2
        end / 64
        α, β = iso_average_mandel(M)
        @test avg ≈ α * J + β * (I - J)
        # Rotations alone suffice for a major-symmetric matrix.
        Ms = M + transpose(M)
        avg = sum(1:64) do m
            Q = mandel_rotation(2π * m / 64)
            Q * Ms * transpose(Q)
        end / 64
        α, β = iso_average_mandel(Ms)
        @test avg ≈ α * J + β * (I - J)

        @test_throws DimensionMismatch iso_average_mandel(zeros(4, 4))
        @test_throws DimensionMismatch iso_average_mandel(zeros(2, 3))
        @test_throws DimensionMismatch iso_average_mandel(zeros(SMatrix{4, 4}))
        @test_throws ArgumentError iso_average_mandel(ones(1, 1))
    end
end
