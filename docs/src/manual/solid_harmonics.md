# [Solid harmonics and the Lamé solutions](@id man-solid-harmonics)

!!! info "See also"
    **Theory**: [Solid harmonics, vector harmonics and the Lamé solutions](@ref th-solid-harmonics) · **API**: [Solid harmonics](@ref api-solid-harmonics)

The complex solid harmonics expand a harmonic field about a point; the vector
harmonics and the solutions of the Lamé equation built on them expand the
displacement of an isotropic elastic solid. This page computes them, translates
them from one center to another and sums them over a periodic array.

## Tables of harmonics at a point

[`regular_solid_harmonics`](@ref) and [`irregular_solid_harmonics`](@ref)
return a [`SolidHarmonicTable`](@ref) of all the harmonics of degree
``n \le N`` at a point, indexed by `[n, m]` with ``|m| \le n``:

```@example solid
using TensND, LinearAlgebra
x = [0.3, -0.7, 0.45]
y = regular_solid_harmonics(x, 2)
Y = irregular_solid_harmonics(x, 2)
(y[1, 0], y[1, 1], Y[1, 0] * norm(x)^3)       # x₃, (x₁ + i x₂)/2, x₃
```

An order beyond the degree reads as zero, and the negative orders follow from
the conjugates:

```@example solid
(y[1, 2], y[1, -1] ≈ -conj(y[1, 1]))
```

## From one center to another

A decaying harmonic about one center is the series of the regular ones about
another, whose coefficients are decaying harmonics of the vector that joins the
two centers:

```@example solid
R, r = [1.1, -0.6, 0.9], [0.15, 0.2, -0.1]
YR, yr = irregular_solid_harmonics(R, 22), regular_solid_harmonics(r, 20)
series = sum((-1)^(n + m) * YR[2 + n, 1 - m] * yr[n, m] for n in 0:20 for m in -n:n)
series ≈ irregular_solid_harmonics(R + r, 2)[2, 1]
```

## Vector harmonics and Lamé solutions

[`regular_vector_harmonic`](@ref) and [`irregular_vector_harmonic`](@ref) read
the vector harmonics of a family from a table that reaches the degree
``n + 1``; [`regular_lame_solution`](@ref) and
[`irregular_lame_solution`](@ref) give the displacements of an isotropic solid
of Poisson ratio ``\nu``. The lowest ones are the dilatation and the uniform
strains:

```@example solid
ν = 0.3
(regular_lame_solution(3, 0, 0, x, ν) ≈ (4ν - 2) / 3 * x,
 regular_lame_solution(1, 2, 0, x, ν) ≈ [-x[1] / 2, -x[2] / 2, x[3]])
```

[`lame_solution_translation`](@ref) gives the coefficients that carry a
decaying solution to the regular ones about another center, as the two parts
of an affine function of ``4(1-\nu)``:

```@example solid
nmax = 24
YR = irregular_solid_harmonics(R, nmax + 3)
series = zeros(ComplexF64, 3)
for n in 0:nmax, m in -n:n, j in 1:3
    a, b = lame_solution_translation(3, j, 2, n, 1, m, YR, R)
    global series += (-1)^(n + m) * (a + 4 * (1 - ν) * b) * regular_lame_solution(j, n, m, r, ν)
end
series ≈ irregular_lame_solution(3, 2, 1, R + r, ν)
```

## Periodic arrays

[`lattice_solid_harmonics`](@ref) sums the decaying harmonics over the images
of a source repeated on an orthogonal lattice. At the source itself, the
degree 0 of the simple cubic lattice of unit period is the constant of its
periodic potential:

```@example solid
S = lattice_solid_harmonics(zeros(3), (1.0, 1.0, 1.0), 4)
real(S[0, 0])
```

and the cubic symmetry leaves no harmonic of degree 2:

```@example solid
abs(S[2, 0]) < 1e-13
```

[`lattice_solid_biharmonics`](@ref) returns, with the same table, the sums
that the translation of the Lamé solutions needs over the images; with
`lame_solution_translation(i, j, n′, n, m′, m, Y, B)` the coefficients are
those of the sum of a solution over the array.
