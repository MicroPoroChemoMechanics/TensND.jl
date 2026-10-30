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
