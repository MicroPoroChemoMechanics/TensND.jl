# Tensor projection onto symmetry subspaces

This tutorial demonstrates how to use the projection tools in TensND.jl to
approximate tensors by their closest isotropic, transversely isotropic (TI),
or orthotropic (ORTHO) counterpart.  These capabilities are useful in
micromechanics and computational homogenization, where effective properties
often need to be identified as belonging to a specific material symmetry class.

## Mathematical background

Given a tensor ``A`` (2nd or 4th order), the projection onto a symmetry class
``\mathcal{S}`` solves

```math
B^* = \arg\min_{B \in \mathcal{S}} \lVert B - A \rVert_F
```

where ``\lVert \cdot \rVert_F`` is the Frobenius norm.  The result is the
orthogonal projection of ``A`` onto the linear subspace ``\mathcal{S}``.

For **transverse isotropy** (TI), the subspace is parametrized by a symmetry
axis ``\mathbf{n}`` and either 2 scalars (order 2) or 5 Walpole coefficients
(order 4).

For **orthotropy** (ORTHO), the subspace is parametrized by a material frame
``(\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3)`` and 9 elastic constants (order 4)
or 3 diagonal entries (order 2).

## Setup

```@example proj
using TensND, LinearAlgebra
```

## Part 1 — 2nd-order TI tensors

### Constructing a TI tensor of order 2

A 2nd-order TI tensor ``\mathbf{A} = a\,\mathbf{n}_T + b\,\mathbf{n}_n`` is
built with the `TensTI{2}` constructor:

```@example proj
n = [0., 0., 1.]   # symmetry axis = e₃

A = TensTI{2}(5.0, 8.0, n)
println("data = ", A.data)
println("trace = ", tr(A))
println("isISO = ", isISO(A), "  isTI = ", isTI(A))
```

The 3×3 matrix is ``\mathrm{diag}(a, a, b)`` when ``\mathbf{n} = \mathbf{e}_3``:

```@example proj
getarray(A)
```

### Projecting onto the TI subspace

Starting from a matrix that is *almost* TI, project it:

```@example proj
M = Float64[5.1 0.3 0; 0.3 4.9 0; 0 0 8]
B, d, drel = proj_tens(:TI, M, n)
println("Projected: a = ", B.data[1], ", b = ", B.data[2])
println("Distance: d = ", round(d, sigdigits=4), ", drel = ", round(drel, sigdigits=4))
```

The projected tensor is exactly TI by construction:

```@example proj
getarray(B)
```

### TI projection with a tilted axis

The axis does not have to be aligned with a canonical direction:

```@example proj
n45 = [1/√2, 0., 1/√2]
A45 = TensTI{2}(3.0, 7.0, n45)
M45 = getarray(A45)
println("Original matrix:")
display(round.(M45, digits=4))
```

Projecting back onto TI with the same axis recovers the original parameters:

```@example proj
B45, d45, _ = proj_tens(:TI, M45, n45)
println("Recovered: a = ", round(B45.data[1], digits=10), ", b = ", round(B45.data[2], digits=10))
println("Distance = ", d45)
```

## Part 2 — 4th-order TI projection

### Round-trip: TI tensor → array → projection

```@example proj
C_ti = tensTI(10., 3., 2.5, 12., 2., n)
A4 = getarray(C_ti)

B4, d4, drel4 = proj_tens(:TI, A4, n)
println("Relative distance = ", drel4)
println("Walpole coefficients: ", collect(argTI(B4)))
```

### Projection of an isotropic tensor onto TI

An isotropic tensor is a special case of TI — the projection should give zero distance:

```@example proj
𝕀, 𝕁, 𝕂 = ISO(Val(3), Val(Float64))
k, μ = 10.0, 5.0
C_iso = 3k * 𝕁 + 2μ * 𝕂
B_iso, d_iso, drel_iso = proj_tens(:TI, getarray(C_iso), n)
println("Relative distance (ISO → TI) = ", drel_iso)
```

### Measuring the departure from TI symmetry

Start with a TI tensor and add a perturbation that breaks the symmetry:

```@example proj
A_pert = copy(A4)
A_pert[1,1,1,1] += 2.0   # break TI symmetry
A_pert[2,2,2,2] -= 1.0

B_pert, d_pert, drel_pert = proj_tens(:TI, A_pert, n)
println("Perturbation distance = ", round(d_pert, sigdigits=4))
println("Relative error = ", round(drel_pert, sigdigits=4))
```

## Part 3 — Orthotropic projection

### Round-trip with a canonical frame

```@example proj
frame = CanonicalBasis{3,Float64}()
t_ort = TensOrtho(10., 8., 12., 3., 2.5, 1.5, 2., 3., 3.5, frame)

Bo, do_, drelo = proj_tens(:ORTHO, getarray(t_ort), frame)
println("ORTHO round-trip relative distance = ", drelo)
```

### ORTHO projection of a TI tensor

A TI tensor is a special case of orthotropy:

```@example proj
B_ti_ortho, d_ti_ortho, _ = proj_tens(:ORTHO, A4, frame)
println("TI → ORTHO distance = ", round(d_ti_ortho, sigdigits=4))
```

### Orthotropic projection with a rotated frame

```@example proj
using Rotations
frame_rot = RotatedBasis(0.3, 0.5, 0.7)
t_rot = TensOrtho(10., 8., 12., 3., 2.5, 1.5, 2., 3., 3.5, frame_rot)

Br, dr, drelr = proj_tens(:ORTHO, getarray(t_rot), frame_rot)
println("Rotated frame round-trip: drel = ", drelr)
```

### 2nd-order ORTHO projection

For a 2nd-order tensor, the orthotropic projection in a given frame keeps only the
diagonal entries:

```@example proj
M_full = Float64[5 1 2; 1 8 3; 2 3 12]
Bm, dm, drelm = proj_tens(:ORTHO, M_full, frame)
println("Projected matrix:")
display(round.(Bm, digits=4))
println("Off-diagonal terms removed → d = ", round(dm, sigdigits=4))
```

## Part 4 — Best symmetry detection

The function `best_sym_tens` tries symmetries from the most restrictive to the
least and returns the first whose relative error falls below a threshold:

```@example proj
# A TI tensor should be detected as TI
_, _, drel_ti, sym_ti = best_sym_tens(C_ti, n; proj = (:TI, :ORTHO))
println("TI tensor → detected symmetry: ", sym_ti, " (drel = ", drel_ti, ")")

# An ORTHO tensor should be detected as ORTHO
_, _, drel_ort, sym_ort = best_sym_tens(t_ort, frame; proj = (:TI, :ORTHO))
println("ORTHO tensor → detected symmetry: ", sym_ort, " (drel = ", round(drel_ort, sigdigits=4), ")")
```

## Part 5 — Rotation-optimized projection (NLopt extension)

When the optimal axis (TI) or frame (ORTHO) is **unknown**, calling `proj_tens`
without providing it triggers a global optimization over the rotation angles.
This requires loading the `NLopt` package:

```julia
using NLopt   # activates the TensNDNLoptExt extension

# Build a TI tensor with a non-trivial axis
n_tilt = [1/√3, 1/√3, 1/√3]
C_tilt = tensTI(10., 3., 2.5, 12., 2., n_tilt)

# Optimise the axis — should recover n_tilt
B_opt, d_opt, drel_opt = proj_tens(:TI, getarray(C_tilt))
println("Optimized relative distance = ", drel_opt)

# Optimise the ORTHO frame
B_ort_opt, d_ort_opt, drel_ort_opt = proj_tens(:ORTHO, getarray(C_tilt))
println("ORTHO optimized relative distance = ", drel_ort_opt)
```

The optimizer uses a two-pass strategy matching the ECHOES C++ library:

1. **Pass 1**: global search (`GD_MLSL`) with a local sub-optimizer (`LD_TNEWTON`),
   coarse tolerances (`xtol = 1e-2`, `ftol = 1e-3`), up to 1000 evaluations.
2. **Pass 2**: local refinement (`LD_TNEWTON`), fine tolerances (`xtol = ftol = 1e-6`),
   up to 100 evaluations.

Gradients are computed automatically via `ForwardDiff.jl`.

## Summary of the projection API

| Function | Description |
| -------- | ----------- |
| `proj_tens(:TI, A, n)` | TI projection with fixed axis `n` |
| `proj_tens(:TI, A)` | TI projection, axis optimized (requires NLopt) |
| `proj_tens(:ORTHO, A, frame)` | ORTHO projection with fixed frame |
| `proj_tens(:ORTHO, A, frame)` | ORTHO projection, frame optimized (requires NLopt) |
| `best_sym_tens(t, n_or_frame)` | Best symmetry with fixed axis/frame |
| `best_sym_tens(t)` | Best symmetry, optimized (requires NLopt) |

All projection functions return `(B, d, drel)` (or `(B, d, drel, sym)` for `best_sym_tens`).
