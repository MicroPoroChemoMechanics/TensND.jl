# Distances between symmetric tensors and the closest tensor of a symmetry class
# in the sense of each. References, independent of the code under test:
#   1. the functions of a matrix of LinearAlgebra (`log`, `exp`, `sqrt` of a
#      `Symmetric`), and the closed forms on the coefficients of structured tensors;
#   2. the defining properties of a distance and the invariance under inversion;
#   3. the stiffness measured by François, Geymonat & Berthaud (1998) and the
#      moduli of its closest isotropic tensors published by Moakher & Norris
#      (2006) and recalled by Morin, Gilormini & Derrien (2020);
#   4. a cubic tensor, which commutes with every isotropic one, so that its
#      closest isotropic tensor has a closed form for the Riemannian distance too;
#   5. the stationarity of the squared distance at the tensor returned.

using TensND, Test, LinearAlgebra, ForwardDiff, Random
using JSON
import Symbolics
using SymPyPythonCall: symbols, simplify

# The stiffness of François et al. (1998), in the Kelvin-Mandel notation, GPa,
# and the published moduli of its closest isotropic tensors.
const _TD_DATA = JSON.parsefile(joinpath(pkgdir(TensND), "data", "literature", "francoisIJSS1998.json"))
const _TD_KM = Float64[_TD_DATA["tables"]["kelvin_mandel_GPa"][i][j] for i in 1:6, j in 1:6]
_td_stiffness() = from_mandel(_TD_KM)
_td_random_spd(rng, n = 6) = (M = randn(rng, n, n); from_mandel(M * M' + n * I))
_td_km(t) = Matrix(mandel(t, CanonicalBasis{3, Float64}()))

@testset "function of a symmetric tensor" begin
    rng = MersenneTwister(20261122)
    C = _td_random_spd(rng)
    M = _td_km(C)
    for (f, reference) in ((log, log(Symmetric(M))), (exp, exp(Symmetric(M / 10))), (sqrt, sqrt(Symmetric(M))), (inv, inv(M)))
        t = f === exp ? C / 10 : C
        @test _td_km(tensor_function(f, t)) ≈ reference rtol = 1.0e-11
    end
    @test get_array(tensor_function(exp, tensor_function(log, C))) ≈ get_array(C) rtol = 1.0e-11
    # order 2
    A = Tens([4.0 1.0 0.5; 1.0 3.0 0.2; 0.5 0.2 5.0])
    @test get_array(tensor_function(sqrt, A)) ≈ sqrt(Symmetric(get_array(A))) rtol = 1.0e-12
    # structured tensors: the function of their coefficients, the class kept
    Ciso = iso_stiffness(10.0, 4.0)
    @test tensor_function(log, Ciso) isa TensISO{4, 3}
    @test get_array(tensor_function(log, Ciso)) ≈ get_array(tensor_function(log, Tens(get_array(Ciso)))) rtol = 1.0e-12
    @test get_data(tensor_function(sqrt, iso_stiffness(Val(2), 10.0, 4.0))) == (sqrt(20.0), sqrt(8.0))
    frame = RotatedBasis(0.3, 0.9, -0.4)
    Ccub = TensCubic(30.0, 8.0, 12.0, frame)
    @test tensor_function(log, Ccub) isa TensCubic
    @test get_array(tensor_function(log, Ccub)) ≈ get_array(tensor_function(log, Tens(get_array(Ccub)))) rtol = 1.0e-11
    n = normalize([0.3, -0.5, 0.8])
    Cti = ti_from_constants(11.0, 4.0, 3.0, 20.0, 2.5, n)
    @test tensor_function(log, Cti) isa TensTI{4}
    @test get_array(tensor_function(log, Cti)) ≈ get_array(tensor_function(log, Tens(get_array(Cti)))) rtol = 1.0e-11
    # the two coupled eigenvalues equal: the limit of the divided difference
    Cdeg = TensTI{4}(7.0, 7.0, 0.0, 3.0, 2.0, n)
    @test get_array(tensor_function(log, Cdeg)) ≈ get_array(tensor_function(log, Tens(get_array(Cdeg)))) rtol = 1.0e-11
    @inferred tensor_function(log, Ciso)
    @inferred tensor_function(log, Cti)
    # derivatives: d/dt f(C + t D) by the divided differences, against the series of exp
    D = _td_random_spd(rng)
    g(t) = _td_km(tensor_function(exp, (C + t * D) / 50))
    h = 1.0e-4
    richardson = (8 * (g(h) - g(-h)) - (g(2h) - g(-2h))) / (12h)
    @test ForwardDiff.derivative(g, 0.0) ≈ richardson rtol = 1.0e-8
    # at a repeated eigenvalue, where the eigenvectors have no derivative: d/dt ln(α(t) 𝕁 + β 𝕂 + t D)
    iso(t) = _td_km(tensor_function(log, Tens(get_array(Ciso)) + t * D))
    richardson = (8 * (iso(h) - iso(-h)) - (iso(2h) - iso(-2h))) / (12h)
    @test ForwardDiff.derivative(iso, 0.0) ≈ richardson rtol = 1.0e-7
    # on the coefficients of a structured tensor, exactly
    @test ForwardDiff.derivative(k -> get_data(tensor_function(log, iso_stiffness(k, 4.0)))[1], 10.0) ≈ 1 / 10.0
    # symbols on structured tensors; a general symbolic tensor is refused
    a, b = symbols("α β", positive = true)
    @test all(simplify.(get_data(tensor_function(log, TensISO{3}(a, b))) .- (log(a), log(b))) .== 0)
    Symbolics.@variables x y
    @test isequal(get_data(tensor_function(inv, TensISO{3}(x, y))), (inv(x), inv(y)))
    @test_throws ArgumentError tensor_function(log, Tens(get_array(TensISO{3}(x, y))))
    # what is refused: an order, a tensor without the major symmetry
    @test_throws ArgumentError tensor_function(log, Tens([1.0, 2.0, 3.0]))
    @test_throws ArgumentError tensor_function(log, from_mandel(M + [i == 1 && j == 2 ? 5.0 : 0.0 for i in 1:6, j in 1:6]))
end

@testset "distances" begin
    rng = MersenneTwister(1998)
    A, B, C = (_td_random_spd(rng) for _ in 1:3)
    E₀ = 200.0
    distances = (
        EuclideanDistance(), LogEuclideanDistance(), RiemannianDistance(),
        ArctanEuclideanDistance(E₀), PowerEuclideanDistance(1, E₀), PowerEuclideanDistance(1 / 2, E₀),
    )
    for d in distances
        @test tensor_distance(A, A, d) ≤ 1.0e-10
        @test tensor_distance(A, B, d) > 0
        @test tensor_distance(A, B, d) ≈ tensor_distance(B, A, d)
        @test tensor_distance(A, B, d) ≤ tensor_distance(A, C, d) + tensor_distance(C, B, d)
    end
    # invariance under inversion: all but the Euclidean distance
    for d in distances[2:end]
        @test tensor_distance(inv(A), inv(B), d isa FunctionEuclideanDistance ? typeof(d)(d.f, d.inverse, inv(d.gauge)) : d) ≈ tensor_distance(A, B, d) rtol = 1.0e-9
    end
    @test !(tensor_distance(inv(A), inv(B)) ≈ tensor_distance(A, B))
    # the default, and its value
    @test tensor_distance(A, B) ≈ norm(_td_km(A) - _td_km(B))
    @test tensor_distance(A, B, LogEuclideanDistance()) ≈ norm(log(Symmetric(_td_km(A))) - log(Symmetric(_td_km(B)))) rtol = 1.0e-10
    # Riemannian: the eigenvalues of A⁻¹B, and the invariance under a congruence
    μ = real.(eigvals(_td_km(A) \ _td_km(B)))
    @test tensor_distance(A, B, RiemannianDistance()) ≈ sqrt(sum(abs2, log.(μ))) rtol = 1.0e-10
    G = randn(rng, 6, 6) + 3I
    congruent(t) = from_mandel(G * _td_km(t) * G')
    @test tensor_distance(congruent(A), congruent(B), RiemannianDistance()) ≈ tensor_distance(A, B, RiemannianDistance()) rtol = 1.0e-9
    # the log-Euclidean distance does not see the unit; the gauge of a function distance does
    @test tensor_distance(3 * A, 3 * B, LogEuclideanDistance()) ≈ tensor_distance(A, B, LogEuclideanDistance()) rtol = 1.0e-11
    @test tensor_distance(3 * A, 3 * B, ArctanEuclideanDistance(3 * E₀)) ≈ tensor_distance(A, B, ArctanEuclideanDistance(E₀)) rtol = 1.0e-11
    # isotropic tensors: closed forms, with numbers and with symbols
    C₁, C₂ = iso_stiffness(10.0, 4.0), iso_stiffness(14.0, 3.0)
    @test tensor_distance(C₁, C₂) ≈ sqrt((30 - 42)^2 + 5 * (8 - 6)^2)
    @test tensor_distance(C₁, C₂) ≈ tensor_distance(Tens(get_array(C₁)), Tens(get_array(C₂)))
    @test tensor_distance(C₁, C₂, LogEuclideanDistance()) ≈ sqrt(log(42 / 30)^2 + 5 * log(6 / 8)^2)
    @test tensor_distance(C₁, C₂, RiemannianDistance()) ≈ tensor_distance(C₁, C₂, LogEuclideanDistance())
    @test tensor_distance(C₁, C₂, RiemannianDistance()) ≈ tensor_distance(Tens(get_array(C₁)), Tens(get_array(C₂)), RiemannianDistance()) rtol = 1.0e-10
    # two moduli: the angle between the two lines of slopes E₁/E₀ and E₂/E₀
    @test tensor_distance(TensISO{3}(3.0, 3.0), TensISO{3}(5.0, 5.0), ArctanEuclideanDistance(2.0)) ≈ sqrt(6) * (atan(5 / 2) - atan(3 / 2))
    a, b, c, e = symbols("a b c e", positive = true)
    @test simplify(tensor_distance(TensISO{3}(a, b), TensISO{3}(c, e), LogEuclideanDistance())^2 - (log(a / c)^2 + 5 * log(b / e)^2)) == 0
    # derivative of a distance, through the eigenvalues
    D = _td_random_spd(rng)
    for d in distances
        φ(t) = tensor_distance(A, B + t * D, d)
        h = 1.0e-4
        @test ForwardDiff.derivative(φ, 0.0) ≈ (8 * (φ(h) - φ(-h)) - (φ(2h) - φ(-2h))) / (12h) rtol = 1.0e-7
    end
    @inferred tensor_distance(C₁, C₂, LogEuclideanDistance())
    @test_throws ArgumentError PowerEuclideanDistance(-1, 1.0)
end

@testset "closest tensor of a symmetry class" begin
    C = _td_stiffness()
    S = inv(C)
    # Euclidean: the projection, on the stiffness and on the compliance, against
    # the contractions with 𝕁 and 𝕂 written on the Kelvin-Mandel matrix
    J = [i ≤ 3 && j ≤ 3 ? 1 / 3 : 0.0 for i in 1:6, j in 1:6]
    K = I - J
    moduli(M) = (sum(J .* M) / 3, sum(K .* M) / 10)                 # 3k = 𝕁 ∷ ℂ, 2μ = 𝕂 ∷ ℂ / 5
    primal = moduli(_TD_KM)
    dual = (inv(9 * sum(J .* inv(_TD_KM)) / 3), inv(4 * sum(K .* inv(_TD_KM)) / 10))
    B, d = closest_tensor(C, :ISO)
    @test B isa TensISO{4, 3}
    @test get_data(B) == get_data(best_fit_iso(C))
    @test d ≈ tensor_distance(C, B)
    @test all(bulk_shear(B) .≈ primal)
    @test all(bulk_shear(inv(first(closest_tensor(S, :ISO)))) .≈ dual)
    # the published moduli (GPa): 170.11 and 96.87, 169.33 and 55.81; the last one is
    # 55.78 here, by both routes
    published = _TD_DATA["quantities"]
    @test primal[1] ≈ published["k_primal_GPa"] atol = 5.0e-3
    @test primal[2] ≈ published["mu_primal_GPa"] atol = 5.0e-3
    @test dual[1] ≈ published["k_dual_GPa"] atol = 5.0e-3
    @test dual[2] ≈ published["mu_dual_GPa"] atol = 3.0e-2
    # invariant under inversion: the tensor closest to the compliance is the inverse
    for distance in (LogEuclideanDistance(), RiemannianDistance(), ArctanEuclideanDistance(150.0), PowerEuclideanDistance(1, 150.0))
        Bc, dc = closest_tensor(C, :ISO; distance)
        inverted = distance isa FunctionEuclideanDistance ? typeof(distance)(distance.f, distance.inverse, inv(distance.gauge)) : distance
        Bs, ds = closest_tensor(S, :ISO; distance = inverted)
        @test all(isapprox.(get_data(inv(Bs)), get_data(Bc); rtol = 1.0e-8))
        @test ds ≈ dc rtol = 1.0e-8
        # between the two Euclidean ones
        k, μ = bulk_shear(Bc)
        @test dual[1] < k < primal[1]
        @test dual[2] < μ < primal[2]
        # a minimum over the isotropic tensors
        for (δk, δμ) in ((1.0e-3, 0.0), (-1.0e-3, 0.0), (0.0, 1.0e-3), (0.0, -1.0e-3))
            @test tensor_distance(C, iso_stiffness(k * (1 + δk), μ * (1 + δμ)), distance) > dc
        end
    end
    # the limits of the gauge: the two Euclidean projections
    small, large = 1.0e-4, 1.0e7
    @test all(isapprox.(bulk_shear(first(closest_tensor(C, :ISO; distance = PowerEuclideanDistance(1, small)))), primal; rtol = 1.0e-6))
    @test all(isapprox.(bulk_shear(first(closest_tensor(C, :ISO; distance = PowerEuclideanDistance(1, large)))), dual; rtol = 1.0e-6))
    @test all(isapprox.(bulk_shear(first(closest_tensor(C, :ISO; distance = ArctanEuclideanDistance(small)))), dual; rtol = 1.0e-6))
    @test all(isapprox.(bulk_shear(first(closest_tensor(C, :ISO; distance = ArctanEuclideanDistance(large)))), primal; rtol = 1.0e-6))
    # a cubic tensor commutes with the isotropic ones: ln(2μ) = (2 ln β + 3 ln γ)/5 for both
    frame = RotatedBasis(0.3, 0.9, -0.4)
    Ccub = TensCubic(30.0, 8.0, 12.0, frame)
    for distance in (LogEuclideanDistance(), RiemannianDistance())
        Bi, di = closest_tensor(Ccub, :ISO; distance)
        @test all(isapprox.(get_data(Bi), (30.0, exp((2 * log(8.0) + 3 * log(12.0)) / 5)); rtol = 1.0e-10))
        @test di ≈ sqrt(2 * log(8.0 / get_data(Bi)[2])^2 + 3 * log(12.0 / get_data(Bi)[2])^2) rtol = 1.0e-9
    end
    # the other classes, with their axis or frame: in the class, and closer than its neighbors
    n = normalize([0.3, -0.5, 0.8])
    for distance in (LogEuclideanDistance(), RiemannianDistance(), PowerEuclideanDistance(1 / 2, 150.0))
        Bti, dti = closest_tensor(C, :TI, n; distance)
        @test Bti isa TensTI{4}
        @test isparallel(axis(Bti), n)
        for k in 1:5
            data = collect(get_data(Bti))
            data[k] *= 1 + 1.0e-3
            @test tensor_distance(C, TensTI{4}(data..., n), distance) > dti
        end
        Bcub, dcub = closest_tensor(C, :CUBIC, frame; distance)
        @test Bcub isa TensCubic
        Bo, dortho = closest_tensor(C, :ORTHO, frame; distance)
        @test Bo isa TensOrtho
        # the classes are nested: iso ⊂ cubic ⊂ ortho (same frame), iso ⊂ TI
        diso = last(closest_tensor(C, :ISO; distance))
        @test dortho ≤ dcub ≤ diso
        @test dti ≤ diso
    end
    # a tensor of the class is its own closest tensor
    Cti = ti_from_constants(11.0, 4.0, 3.0, 20.0, 2.5, n)
    for distance in (EuclideanDistance(), LogEuclideanDistance(), RiemannianDistance(), ArctanEuclideanDistance(10.0))
        Bti, dti = closest_tensor(Tens(get_array(Cti)), :TI, n; distance)
        @test get_array(Bti) ≈ get_array(Cti) rtol = 1.0e-9
        @test dti ≤ 1.0e-8
    end
    # derivatives through the closed forms; symbols on a structured tensor
    kiso(x) = bulk_modulus(first(closest_tensor(TensCubic(x, 8.0, 12.0, frame), :ISO; distance = LogEuclideanDistance())))
    @test ForwardDiff.derivative(kiso, 30.0) ≈ 1 / 3
    Symbolics.@variables α β γ
    Bs, _ = closest_tensor(TensCubic(α, β, γ, CanonicalBasis{3, Float64}()), :ISO; distance = LogEuclideanDistance())
    at = Dict(α => 30.0, β => 8.0, γ => 12.0)
    number(x) = Symbolics.value(Symbolics.substitute(x, at; fold = Val(true)))
    @test number(get_data(Bs)[1]) ≈ 30.0
    @test number(get_data(Bs)[2]) ≈ exp((2 * log(8.0) + 3 * log(12.0)) / 5)
    # the Riemannian iteration takes numbers only
    @test_throws ArgumentError closest_tensor(TensCubic(ForwardDiff.Dual(30.0, 1.0), ForwardDiff.Dual(8.0, 0.0), ForwardDiff.Dual(12.0, 0.0), frame), :ISO; distance = RiemannianDistance())
end

@testset "in the plane" begin
    # the same functions on tensors of dimension 2: Kelvin-Mandel matrices of size 3
    M = [5.0 1.0 0.3; 1.0 4.0 -0.2; 0.3 -0.2 3.0]
    N = [6.0 0.5 0.1; 0.5 5.0 0.4; 0.1 0.4 2.5]
    A, B = from_mandel(M), from_mandel(N)
    km(t) = Matrix(mandel(t, CanonicalBasis{2, Float64}()))
    @test km(tensor_function(log, A)) ≈ log(Symmetric(M)) rtol = 1.0e-12
    @test tensor_distance(A, B) ≈ norm(M - N)
    @test tensor_distance(A, B, LogEuclideanDistance()) ≈ norm(log(Symmetric(M)) - log(Symmetric(N))) rtol = 1.0e-12
    @test tensor_distance(A, B, RiemannianDistance()) ≈ sqrt(sum(abs2, log.(real.(eigvals(M \ N))))) rtol = 1.0e-10
    for d in (LogEuclideanDistance(), RiemannianDistance(), ArctanEuclideanDistance(4.0))
        inverted = d isa FunctionEuclideanDistance ? typeof(d)(d.f, d.inverse, inv(d.gauge)) : d
        @test tensor_distance(inv(A), inv(B), inverted) ≈ tensor_distance(A, B, d) rtol = 1.0e-10
        Bi, di = closest_tensor(A, :ISO; distance = d)
        Si, _ = closest_tensor(inv(A), :ISO; distance = inverted)
        @test Bi isa TensISO{4, 2}
        @test all(isapprox.(get_data(inv(Si)), get_data(Bi); rtol = 1.0e-8))
        @test di ≤ tensor_distance(A, TensISO{2}((get_data(Bi) .* (1.001, 0.999))...), d)
    end
    # in the plane 𝕁 = 𝟏 ⊗ 𝟏 / 2 and 𝕂 has two dimensions: ln α = 𝕁 :: ln A, ln β = 𝕂 :: ln A / 2
    L = log(Symmetric(M))
    α, β = get_data(first(closest_tensor(A, :ISO; distance = LogEuclideanDistance())))
    @test log(α) ≈ (L[1, 1] + L[1, 2] + L[2, 1] + L[2, 2]) / 2
    @test log(β) ≈ (tr(L) - log(α)) / 2
    # second-order tensors of the plane
    a, b = [4.0 1.0; 1.0 3.0], [2.0 0.5; 0.5 5.0]
    @test tensor_distance(Tens(a), Tens(b), LogEuclideanDistance()) ≈ norm(log(Symmetric(a)) - log(Symmetric(b))) rtol = 1.0e-12
    Ki, _ = closest_tensor(Tens(a), :ISO; distance = LogEuclideanDistance())
    @test only(get_data(Ki)) ≈ sqrt(det(a))                              # the geometric mean of the eigenvalues
end
