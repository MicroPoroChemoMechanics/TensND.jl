# The battery of the bitwise lock (`test/test_golden.jl`): tensors of every
# structured class and generic ones, in canonical, rotated and orthogonal
# bases, and the outputs of the operations a refactoring may touch. Each case
# is a name and a function returning a vector of Float64 values; a case that
# throws records the type of the exception instead, so that a path which starts
# or stops failing is caught as well.

using TensND, LinearAlgebra

const _N = [1.0, 2.0, 2.0] / 3                      # a unit TI axis off every frame axis
const _R = RotatedBasis(0.3, 0.7, 1.1)              # a material frame
const _R2 = RotatedBasis(1.2, -0.4, 0.5)            # a second, to change bases into
const _CANON = CanonicalBasis{3, Float64}()

_gen_array() = [
    2.0 * (i == j) * (k == l) + 1.5 * ((i == k) * (j == l) + (i == l) * (j == k)) +
        0.1 * (i + j) * (k + l) + 0.01 * i * j * k * l
        for i in 1:3, j in 1:3, k in 1:3, l in 1:3
]

_tensors() = (
    iso = TensISO{3}(3 * 20.0, 2 * 8.0),
    ti = TensTI{4}(30.0, 9.0, 7.0, 6.0, 4.0, _N),
    ortho = TensOrtho(30.0, 25.0, 20.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, _R),
    cubic = TensCubic(60.0, 16.0, 10.0, _R),
    gen = Tens(_gen_array()),
    gen_rot = Tens(_gen_array(), _R),
)

_vals(x::Number) = Float64[x]
_vals(x::AbstractArray) = Float64.(vec(collect(x)))
_vals(t::AbstractTens) = _vals(get_array(t))
_canon(t) = _vals(express(t, _CANON))

function golden_cases()
    T = _tensors()
    S = Tens([2.0 0.5 0.1; 0.5 3.0 0.2; 0.1 0.2 4.0])
    cases = Pair{String, Function}[]
    for (name, t) in pairs(T)
        push!(cases, "$name/own" => () -> _vals(t))
        push!(cases, "$name/canonical" => () -> _canon(t))
        push!(cases, "$name/other_basis" => () -> _vals(express(t, _R2)))
        push!(cases, "$name/KM" => () -> _vals(mandel(t)))
        push!(cases, "$name/inv" => () -> _canon(inv(t)))
        push!(cases, "$name/dcontract_order2" => () -> _canon(t ⊡ S))
        push!(cases, "$name/acoustic" => () -> _canon(dotdot(Tens(_N), t, Tens(_N))))
        push!(cases, "$name/getindex" => () -> Float64[t[1, 1, 1, 1], t[1, 2, 1, 2], t[2, 3, 1, 1], t[3, 1, 2, 3]])
        push!(cases, "$name/best_fit_iso" => () -> _canon(best_fit_iso(t)))
        push!(cases, "$name/best_fit_ti" => () -> _canon(best_fit_ti(t, _N)))
    end
    for (a, b) in (
            (:iso, :ti), (:ti, :ti), (:ortho, :ortho), (:iso, :ortho), (:cubic, :ortho),
            (:cubic, :cubic), (:gen, :gen), (:gen_rot, :iso), (:ti, :gen_rot),
        )
        push!(cases, "dcontract/$a⊡$b" => () -> _canon(T[a] ⊡ T[b]))
        push!(cases, "sum/$a+$b" => () -> _canon(T[a] + T[b]))
    end
    push!(cases, "order2/otimes" => () -> _canon(S ⊗ S))
    push!(cases, "order2/otimesu" => () -> _canon(otimesu(S, S)))
    push!(cases, "order2/sotimes" => () -> _canon(sotimes(S, S)))
    push!(cases, "order2/inv" => () -> _canon(inv(S)))
    push!(cases, "order2/rotated" => () -> _vals(express(S, _R)))
    # Elastic waves: the velocities, polarizations and group velocities of each
    # class along an oblique direction, the first-order and generalized
    # parameters, an interface and a prestressed solid.
    n = [0.3, -0.5, 0.8]
    for (name, t) in pairs(T)
        name === :gen_rot && continue
        push!(cases, "waves/$name/plane_waves" => () -> (w = plane_waves(t, 2.5, n); vcat(_vals(w.velocity), _vals(abs.(w.polarization)), _vals(w.group))))
        push!(cases, "waves/$name/weak" => () -> (w = weak_plane_waves(t, 2.5, n); Float64[w.qP[1], w.qS1[1], w.qS2[1]]))
        push!(cases, "waves/$name/generalized" => () -> _vals(collect(values(generalized_thomsen_parameters(t)))))
    end
    push!(cases, "waves/ti/thomsen" => () -> _vals(collect(values(thomsen_parameters(T.ti, 2.5)))))
    push!(
        cases, "waves/interface" => () -> begin
            f, _ = interface_waves(T.ortho, 2.5, _N, [0.05, 0.02, 0.0])
            r = reflection_transmission(T.ortho, 2.5, T.iso, 2.0, _N, f[1])
            vcat(_vals(abs.(r.R)), _vals(abs.(r.T)), _vals(collect(r.energy[1])), _vals(collect(r.energy[2])))
        end
    )
    push!(
        cases, "waves/third_order" => () -> begin
            s3 = ThirdOrderIsotropic(8.3, 11.9, -3.4e3, -6.8e3, -6.7e3)
            F = [1.001 0.0002 0.0; 0.0 0.9995 0.0003; 0.0001 0.0 1.0004]
            vcat(_vals(prestressed_waves(s3, 2.3, F, n).velocity), _vals(collect(values(acoustoelastic_coefficients(s3)))))
        end
    )
    push!(
        cases, "plane/characteristic_pair" => () -> begin
            s2, p2 = characteristic_pair(from_mandel([30.0 4.0 3.0; 4.0 9.0 -2.0; 3.0 -2.0 7.0]))
            Float64[real(s2), imag(s2), real(p2), imag(p2)]
        end
    )
    return cases
end

# The values of a case as hex strings, or the name of the exception it throws.
function golden_record(f)
    v = try
        f()
    catch err
        return "ERROR:" * string(nameof(typeof(err)))
    end
    return join((string(reinterpret(UInt64, x); base = 16, pad = 16) for x in v), ",")
end
