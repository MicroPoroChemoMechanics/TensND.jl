# =============================================================================
#  test_rotations_any_dimension.jl — orthonormal bases, random rotations and
#  the Kelvin-Mandel rotation matrix in dimensions 2 to 5.
#
#  References, independent of the code under test: the components of a tensor
#  in a rotated frame computed by explicit sums over the indices (OMEinsum),
#  the closed form of the Kelvin-Mandel rotation of the plane, and the moments
#  of the Haar measure of SO(d).
# =============================================================================

using StaticArrays: SMatrix

const RAD_SEED = 20261126

"Components of an array of order 2 or 4 in the frame whose vectors are the columns of `R`."
rad_in_frame(A::AbstractMatrix, R) = ein"ia,jb,ij->ab"(R, R, A)
rad_in_frame(A::AbstractArray{<:Any, 4}, R) = ein"ia,jb,kc,ld,ijkl->abcd"(R, R, R, R, A)

"An array of order 4 with the two minor symmetries."
function rad_minor_symmetric(rng, d)
    A = randn(rng, d, d, d, d)
    A = (A + permutedims(A, (2, 1, 3, 4))) / 2
    return (A + permutedims(A, (1, 2, 4, 3))) / 2
end

"Kelvin-Mandel rotation of the plane, written out for the rotation of angle `θ`."
function rad_plane_bond(θ)
    c, s, r = cos(θ), sin(θ), sqrt(oftype(cos(θ), 2))
    return [c^2 s^2 -r * c * s; s^2 c^2 r * c * s; r * c * s -r * c * s c^2 - s^2]
end

rad_inferred_bond(R) = @inferred mandel_rotation(R)

@testsection "Rotations and orthonormal bases in any dimension" begin
    rng = Xoshiro(RAD_SEED)

    @testsection "random rotations: Haar measure of SO(d)" begin
        for d in 2:5
            R = @inferred random_rotation(rng, Val(d))
            @test R isa SMatrix{d, d, Float64}
            @test R' * R ≈ I
            @test det(R) ≈ 1
            @test random_rotation(rng, Val(d), Float32) isa SMatrix{d, d, Float32}
            @test random_rotation(Val(d)) isa SMatrix{d, d, Float64}
            @test random_rotation(Val(d), Float32) isa SMatrix{d, d, Float32}
            # First and second moments of an entry: 0 and 1/d. The entry squared
            # follows a Beta(1/2, (d-1)/2) law, of standard deviation below 0.36.
            n = 4000
            draws = [random_rotation(rng, Val(d)) for _ in 1:n]
            @test all(Q -> det(Q) > 0, draws)
            for (i, j) in ((1, 1), (d, d), (1, d))
                @test abs(sum(Q -> Q[i, j], draws) / n) < 5 / sqrt(n)
                @test abs(sum(Q -> Q[i, j]^2, draws) / n - 1 / d) < 5 * 0.36 / sqrt(n)
            end
        end
        # In dimension 3 the method with the dimension is the one without it.
        @test random_rotation(Xoshiro(7), Val(3)) == random_rotation(Xoshiro(7))
        @test random_rotation(Xoshiro(7), Val(3), Float32) == random_rotation(Xoshiro(7), Float32)
    end

    @testsection "random axes and bases" begin
        for d in 2:5
            n = @inferred random_axis(rng, Val(d))
            @test n isa NTuple{d, Float64}
            @test norm(n) ≈ 1
            @test random_axis(Val(d), Float32) isa NTuple{d, Float32}
            @test random_axis(rng, Val(d), Float32) isa NTuple{d, Float32}
            ℬ = @inferred random_basis(rng, Val(d))
            @test ℬ isa RotatedBasis{d, Float64}
            @test isorthonormal(ℬ)
            @test det(vecbasis(ℬ, :cov)) ≈ 1
            @test random_basis(Val(d)) isa RotatedBasis{d, Float64}
            @test random_basis(rng, Val(d), Float32) isa RotatedBasis{d, Float32}
        end
        @test random_axis(Xoshiro(7), Val(3)) == random_axis(Xoshiro(7))
        @test vecbasis(random_basis(Xoshiro(7), Val(3)), :cov) == vecbasis(random_basis(Xoshiro(7)), :cov)
        # Uniform on the sphere: the mean of a component is 0, that of its square 1/d.
        axes4 = [random_axis(rng, Val(4)) for _ in 1:4000]
        @test abs(sum(first, axes4) / 4000) < 5 * 0.5 / sqrt(4000)
        @test abs(sum(n -> n[4]^2, axes4) / 4000 - 1 / 4) < 5 * 0.36 / sqrt(4000)
    end

    @testsection "rotated bases beyond dimension 3" begin
        for d in (4, 5)
            R = Matrix(random_rotation(rng, Val(d)))
            ℬ = Basis(R)
            @test ℬ isa RotatedBasis{d, Float64}
            @test ℬ == RotatedBasis(R, Val(d))
            @test angles(ℬ) == (;)                      # no angle convention beyond dimension 3
            @test angles(R) == (;)
            @test vecbasis(ℬ, :cov) == R
            @test vecbasis(ℬ, :cont) == R
            @test isorthonormal(ℬ)
            @test Basis(Matrix{Float64}(I, d, d)) isa CanonicalBasis{d, Float64}
            @test_throws ArgumentError RotatedBasis(R + 0.1 * randn(rng, d, d), Val(d))
            @test_throws DimensionMismatch RotatedBasis(R, Val(d - 1))
            @test occursin("RotatedBasis{$d, Float64}", sprint(show, MIME"text/plain"(), ℬ))

            # Components in the rotated frame, order 2 and order 4.
            a, A = randn(rng, d, d), rad_minor_symmetric(rng, d)
            @test get_array(express(Tens(a), ℬ)) ≈ rad_in_frame(a, R)
            @test get_array(express(Tens(A), ℬ)) ≈ rad_in_frame(A, R)
            @test get_array(express(Tens(rad_in_frame(A, R), ℬ), CanonicalBasis{d, Float64}())) ≈ A
            @test Tens(rad_in_frame(a, R), ℬ) ≈ Tens(a)

            # A permuted basis keeps the orientation; a scaled one is orthogonal.
            σ = ntuple(k -> k == 1 ? 2 : k == 2 ? 1 : k, d)
            ℬσ = permute_basis(ℬ, σ)
            @test ℬσ isa RotatedBasis{d, Float64}
            @test det(vecbasis(ℬσ, :cov)) ≈ 1
            @test vecbasis(ℬσ, :cov)[:, 1] == R[:, 2]
            @test vecbasis(ℬσ, :cov)[:, d] == -R[:, d]
            λ = collect(1.0:d)
            ℬλ = OrthogonalBasis(ℬ, λ)
            @test isorthogonal(ℬλ) && !isorthonormal(ℬλ)
            @test vecbasis(ℬλ, :cov) ≈ R * Diagonal(λ)
        end
    end

    @testsection "Kelvin-Mandel rotation matrix" begin
        for d in 2:5
            R = random_rotation(rng, Val(d))
            n = d * (d + 1) ÷ 2
            Q = rad_inferred_bond(R)
            @test Q isa SMatrix{n, n, Float64}
            @test Q' * Q ≈ I
            @test mandel_rotation(Matrix(R), Val(d)) == Q
            # The rotation of the components and that of the Kelvin-Mandel matrix.
            A = rad_minor_symmetric(rng, d)
            @test Q' * mandel(Tens(A)) * Q ≈ mandel(Tens(rad_in_frame(A, Matrix(R))))
            ε = randn(rng, d, d); ε = (ε + ε') / 2
            @test Q' * mandel(Tens(ε)) ≈ mandel(Tens(rad_in_frame(ε, Matrix(R))))
            # The Kelvin-Mandel matrix of a tensor written in a rotated basis is
            # that of its components in this basis.
            ℬ = Basis(Matrix(R))
            @test mandel(express(Tens(A), ℬ)) ≈ Q' * mandel(Tens(A)) * Q
            # A composition of rotations is the product of their matrices.
            R₂ = random_rotation(rng, Val(d))
            @test mandel_rotation(R * R₂) ≈ Q * mandel_rotation(R₂)
        end
        # A plain matrix is three-dimensional unless the dimension is given.
        R₃ = random_rotation(rng)
        @test mandel_rotation(Matrix(R₃)) == mandel_rotation(R₃)
        @test (@inferred mandel_rotation(Matrix(R₃))) isa SMatrix{6, 6, Float64}
        @test_throws DimensionMismatch mandel_rotation(rand(2, 2))
        @test_throws DimensionMismatch mandel_rotation(rand(3, 3), Val(4))
        @test mandel_rotation(Tensors.Tensor{2, 3}(Tuple(R₃))) == mandel_rotation(R₃)
    end

    @testsection "Kelvin-Mandel rotation of the plane: closed form, derivatives, symbols" begin
        θ₀ = 0.37
        @test (@inferred mandel_rotation(θ₀)) isa SMatrix{3, 3, Float64}
        @test mandel_rotation(θ₀) ≈ rad_plane_bond(θ₀)
        @test mandel_rotation(θ₀) == mandel_rotation(TensND.rotation_matrix(θ₀))
        @test mandel_rotation(1) ≈ rad_plane_bond(1.0)
        # First and second derivatives with respect to the angle.
        D1 = ForwardDiff.derivative(mandel_rotation, θ₀)
        @test D1 ≈ ForwardDiff.derivative(rad_plane_bond, θ₀)
        D2 = ForwardDiff.derivative(θ -> ForwardDiff.derivative(mandel_rotation, θ), θ₀)
        @test D2 ≈ ForwardDiff.derivative(θ -> ForwardDiff.derivative(rad_plane_bond, θ), θ₀)
        # d/dθ of the entry c² - s² is -2 sin 2θ, d²/dθ² is -4 cos 2θ.
        @test D1[3, 3] ≈ -2sin(2θ₀)
        @test D2[3, 3] ≈ -4cos(2θ₀)
        # Derivative through a matrix of dual numbers, in dimension 4.
        R₄, G = random_rotation(rng, Val(4)), randn(rng, 4, 4)
        K = G - G'                                           # a generator of SO(4)
        curve(t) = R₄ * (I + t * K)                           # tangent to SO(4) at t = 0
        dQ = ForwardDiff.derivative(t -> mandel_rotation(curve(t), Val(4)), 0.0)
        A = rad_minor_symmetric(rng, 4)
        dref = ForwardDiff.derivative(t -> mandel(Tens(rad_in_frame(A, curve(t)))), 0.0)
        Q₄ = mandel_rotation(R₄)
        @test dQ' * mandel(Tens(A)) * Q₄ + Q₄' * mandel(Tens(A)) * dQ ≈ dref

        θ = symbols("θ", real = true)
        Qˢ = mandel_rotation(θ)
        @test Qˢ isa SMatrix{3, 3, <:Sym}
        @test all(iszero, simplify.(Qˢ - rad_plane_bond(θ)))
        @test all(iszero, simplify.(transpose(Qˢ) * Qˢ - I))
        # A symbolic stiffness of the plane, rotated: the isotropic one is invariant.
        k, μ = symbols("k μ", positive = true)
        ℂ = mandel(iso_stiffness(Val(2), k, μ))
        @test all(iszero, simplify.(transpose(Qˢ) * ℂ * Qˢ - ℂ))
    end
end
