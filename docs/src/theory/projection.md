# [Projection onto a symmetry class](@id th-projection)

Given an arbitrary tensor, find the closest one having a prescribed material
symmetry — and, if the orientation is unknown, find that too. This page states
the problem, derives the closed forms `TensND` uses at fixed orientation, and
describes the optimizer used when the orientation is free.

Implementation: `src/tens_projection.jl` and, for the orientation search, the
weak-dependency extension `ext/TensNDNLoptExt.jl`. The formulation follows the
one of the Echoes C++ library; the general theory of closest tensors of
prescribed symmetry is [moakher2006](@cite), and the harmonic-decomposition
alternative is [browaeys2004](@cite).

## The problem

Let ``\mathcal{S}`` be a symmetry class — a **linear subspace** once the
orientation is fixed. With the Frobenius scalar product of
[Tensor algebra](@ref th-tensor-algebra),

```math
\mathbb{B}^\star=\arg\min_{\mathbb{B}\in\mathcal{S}}\ \|\mathbb{B}-\mathbb{A}\| ,
```

so ``\mathbb{B}^\star`` is the **orthogonal projection** of ``\mathbb{A}`` onto
``\mathcal{S}``. [`proj_tens`](@ref) returns the triple

| Returned | Meaning |
| :------- | :------ |
| ``\mathbb{B}`` | the projected tensor, in the storage type of the class |
| ``d`` | absolute distance ``\|\mathbb{B}-\mathbb{A}\|`` |
| ``d_{\text{rel}}`` | relative distance ``d/\|\mathbb{A}\|`` |

Because the Kelvin–Mandel map is an isometry
([Kelvin–Mandel representation](@ref th-kelvin-mandel)), the whole computation is
done on ``6\times6`` matrices.

## Fixed orientation: the normal equations

Write a basis ``(\mathbb{T}_i)_{i=1..N}`` of ``\mathcal{S}`` and
``\mathbb{B}=\sum_i\alpha_i\mathbb{T}_i``. Stationarity of
``\|\mathbb{B}-\mathbb{A}\|^2`` gives the **normal equations**

```math
\sum_{j}\langle\mathbb{T}_i,\mathbb{T}_j\rangle\,\alpha_j^\star
=\langle\mathbb{T}_i,\mathbb{A}\rangle
\qquad\forall i .
```

In general this is an ``N\times N`` linear solve. The reason `TensND` never
performs one is that **in the bases it uses, the Gram matrix is diagonal**.

### Transverse isotropy

The major-symmetric TI subspace is spanned by the symmetrized Walpole basis
``\mathbb{W}^s=(\mathbb{W}_1,\ \mathbb{W}_2,\ \mathbb{W}_3+\mathbb{W}_4,\
\mathbb{W}_5,\ \mathbb{W}_6)``, whose Gram matrix is
([The Walpole basis](@ref th-walpole))

```math
\langle\mathbb{W}^s_i,\mathbb{W}^s_j\rangle=g_i\,\delta_{ij},
\qquad (g_i)=(1,\,1,\,2,\,2,\,2).
```

Each coefficient is therefore an independent quotient,
``\ell_i=\langle\mathbb{A},\mathbb{W}^s_i\rangle/g_i``. Expanding the scalar
products on the ``6\times6`` matrix ``C=\mathrm{Mat}(\mathbb{A})`` in the frame
where ``\underline{n}=\underline{e}_3``, and writing

```math
c=\tfrac{C_{11}+C_{22}}{2},\qquad d=\tfrac{C_{12}+C_{21}}{2},
```

gives the closed forms implemented in `_project_TI_KM`:

```math
\ell_1=C_{33},\qquad
\ell_2=c+d,\qquad
\ell_3=\frac{C_{13}+C_{23}+C_{31}+C_{32}}{2\sqrt2},
```

```math
\ell_5=\frac{c-d+C_{66}}{2},\qquad
\ell_6=\frac{C_{44}+C_{55}}{2}.
```

Every one of them is a plain average of the entries the corresponding
``\mathbb{W}^s_i`` selects — which is what an orthogonal projection onto an
orthogonal basis must look like.

### Orthotropy

Same argument with the nine orthotropic generators, which are likewise mutually
orthogonal. The projection at fixed frame simply **keeps the block-diagonal part
and symmetrizes** (`_project_ORTHO_KM`):

```math
C_{ii}\ \text{unchanged},\qquad
C_{ij}=\tfrac{C_{ij}+C_{ji}}{2}\ (i<j\le3),\qquad
C_{44},C_{55},C_{66}\ \text{halved from the KM diagonal},
```

everything outside the two blocks being discarded. For an order-2 tensor the
same reasoning reduces to *keep the diagonal in the material frame*.

### Isotropy

The two-dimensional case, already given in [Isotropic tensors](@ref th-isotropic):
``\mathrm{ISO}(\mathbb{T})=(\mathbb{T}::\mathbb{J})\,\mathbb{J}
+\tfrac15(\mathbb{T}::\mathbb{K})\,\mathbb{K}``. No orientation is involved, so
this case is complete on its own.

## Free orientation

When the axis or frame is unknown, ``\mathcal{S}`` is no longer a subspace: the
admissible set is the union of the rotated subspaces

```math
\mathcal{S}_{\text{rot}}=\bigl\{\,{}^{t}Q_\Theta\,
\bigl(\textstyle\sum_i\alpha_i\mathbb{T}_i\bigr)\,Q_\Theta
\;:\;\alpha\in\mathbb{R}^N,\ \Theta\in SO(3)\,\bigr\},
```

with ``Q_\Theta`` the orthogonal Kelvin–Mandel rotation of
[Rotations](@ref th-rotations). The problem becomes

```math
\min_{\alpha,\,\Theta}\ J(\alpha,\Theta),
\qquad
J=\frac{\bigl\|\,{}^{t}Q_\Theta\,C\,Q_\Theta-\sum_i\alpha_i\mathbb{T}_i\,\bigr\|^2}
       {\|C\|^2}.
```

### Condensing out the linear part

For a *given* ``\Theta`` the inner minimization over ``\alpha`` is the linear
problem already solved. Substituting its solution, and using that
``\mathbb{B}`` is the orthogonal projection of the rotated tensor — so that
``\|C\|^2=\|\mathbb{B}\|^2+\|C-\mathbb{B}\|^2`` by Pythagoras — the objective
collapses to a function of the orientation alone:

```math
\boxed{\;
j(\Theta)=J\bigl(\alpha^\star(\Theta),\Theta\bigr)
=1-\frac{\bigl\|\mathbb{B}(\Theta)\bigr\|^2}{\|C\|^2}\;}
```

This is literally what `_obj_TI4` and `_obj_ORTHO4` compute: rotate, project,
compare norms. Maximizing the norm of the projection *is* minimizing the
distance, and only two or three angles remain — ``(\theta,\varphi)`` for
transverse isotropy, ``(\theta,\varphi,\psi)`` for orthotropy.

Differentiating the normal equations shows that the gradient needs no derivative
of ``\alpha^\star``:

```math
\frac{\partial j}{\partial\Theta}
=-\frac{4}{\|C\|^2}\,
\Bigl\langle\,\mathbb{B}(\Theta),\;
{}^{t}Q_\Theta\,C\,\frac{\partial Q_\Theta}{\partial\Theta}\Bigr\rangle .
```

In practice `TensND` obtains this gradient by **automatic differentiation**
(`ForwardDiff` through the objective), which is exact to machine precision and
removes any risk of the analytic expression drifting from the code.

### The optimizer

Loading `NLopt` activates `TensNDNLoptExt` and enables the no-orientation
methods of [`proj_tens`](@ref). The strategy is a **deterministic multi-start**:

1. a candidate from the eigenstructure of ``C`` (`_candidate_TI_axis`,
   `_candidate_ORTHO_frame`) — exact whenever the tensor genuinely has the
   symmetry sought;
2. a fixed angular grid containing the canonical axes;
3. `LD_TNEWTON` local refinement from every start, with ForwardDiff gradients;
4. the best objective over all starts *and* all refined starts.

Two properties follow, and both matter:

- **reproducibility** — repeated calls return bit-identical angles;
- **no regression** — because the grid contains the canonical frame, the
  optimized projection is never worse than the fixed-frame projection along any
  grid point.

!!! note "Why the stochastic global pass was removed"
    Earlier versions followed the Echoes C++ strategy: a `GD_MLSL` global stage
    followed by a local one. `GD_MLSL` is stochastic and NLopt seeds it from the
    clock, so results were irreproducible from call to call. On a tensor exactly
    orthotropic about a tilted frame, roughly 0.5 % of runs returned a spurious
    local minimum (``d_{\text{rel}}\approx0.44`` instead of ``\approx10^{-13}``),
    which surfaced as an intermittent CI failure. The multi-start above matched
    or beat the best of 30 `GD_MLSL` runs on every tensor tested — 150 random
    anisotropic ones included — and is 3 to 8 times faster.

## Detecting the symmetry: the cascade

[`best_sym_tens`](@ref) tries the classes from the most restrictive to the least
and returns the first whose *relative* error falls below a tolerance ``\varepsilon``:

```math
\text{ISO}\ \longrightarrow\ \text{TI}\ \longrightarrow\ \text{ORTHO}
\ \longrightarrow\ \text{ANISO} .
```

The relative criterion is what makes the tolerance dimensionless and independent
of the units of the moduli. Two modes:

| `optimize_angles` | TI axis / ORTHO frame | Needs NLopt |
| :---------------- | :-------------------- | :---------- |
| `false` (default) | taken from the tensor if it is a structured container, otherwise from the Kelvin–Mandel eigenstructure | no |
| `true` | found by the multi-start above | yes |

The value-level predicates [`is_ISO`](@ref), [`is_TI`](@ref), [`is_ORTHO`](@ref)
are the same computation with a boolean answer.

## Projection is not averaging

Two different operations are easily confused, and only one of them is on this
page.

| | **Projection** (this page) | **Group average** |
| :--- | :--- | :--- |
| definition | closest element of the subspace | ``\bar{\mathbb{A}}=\int_G{}^{t}Q_g\,\mathbb{A}\,Q_g\,\mathrm{d}g`` |
| result | best *fit* in the class | exact invariant part |
| target space | major-symmetric subspace (5 or 9 parameters) | full commutant (8 parameters for TI) |
| discards | everything orthogonal to the subspace | nothing invariant |

For a major-symmetric input the two coincide. For an input that is **not**
major-symmetric — a strain-concentration tensor, typically — they differ: the
average preserves the ``\ell_3\neq\ell_4`` split and the couplings ``\ell_7``,
``\ell_8`` of [The extended Walpole algebra](@ref th-walpole-extended), whereas
the best fit forces ``\ell_3=\ell_4`` and drops ``\ell_7,\ell_8``. Reporting one
where the other is meant silently deletes physical content.

## The choice of distance

Everything above minimizes the **Euclidean** (Frobenius) distance. That choice is
not neutral: it is not invariant under inversion, so projecting a stiffness and
projecting its compliance give different materials
([Isotropic tensors](@ref th-isotropic)). Distances that repair this —
log-Euclidean, power-Euclidean, arctan-Euclidean — are constructed in
[morin2020](@cite). `TensND` implements the Euclidean one only, so a reported
projection should always state which of ``\mathbb{C}`` or ``\mathbb{S}`` was
projected.
