# ─────────────────────────────────────────────────────────────────────────────
#  waves_splitting.jl — shear-wave splitting.
#
#  In an anisotropic medium the two quasi-transverse waves of a direction
#  travel at different velocities, with orthogonal polarizations: a shear
#  pulse splits into a fast and a slow one. The measure is the relative
#  difference of the two velocities, read on the waves of
#  `waves_christoffel.jl`, and the delay that a path of unit length builds up.
# ─────────────────────────────────────────────────────────────────────────────

"""
    shear_splitting(ℂ::AbstractTens{4,3}, ρ, n) → (splitting, fast, delay)

The splitting of the shear waves that travel along `n`:

- `splitting`: ``(V_{S1} - V_{S2})/V_{S1}``, between 0 and 1, the fast velocity
  as the reference;
- `fast`: the unit polarization of the fast wave;
- `delay`: ``1/V_{S2} - 1/V_{S1}``, the time by which the slow wave lags per
  unit length of path.

Read on [`plane_waves`](@ref); numbers and dual numbers
([Shear-wave splitting](@ref th-splitting)).
"""
function shear_splitting(ℂ::AbstractTens{4, 3}, ρ, n)
    w = plane_waves(ℂ, ρ, n)
    V₁, V₂ = w.velocity[2], w.velocity[3]
    return (splitting = (V₁ - V₂) / V₁, fast = w.polarization[:, 2], delay = inv(V₂) - inv(V₁))
end

export shear_splitting
