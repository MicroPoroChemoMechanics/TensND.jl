# [Projection onto a symmetry class](@id th-projection)

Given an arbitrary tensor, find the closest one having a prescribed material
symmetry — and, if the orientation is unknown, find that too. This page states
the problem, derives the closed forms `TensND` uses at fixed orientation, and
describes the optimizer used when the orientation is free.

Implementation: `src/tens_projection.jl` and, for the orientation search,
`src/orientation_search.jl`. The formulation follows the
one of the Echoes C++ library; the general theory of closest tensors of
prescribed symmetry is due to [moakher2006](@citet), and the
harmonic-decomposition alternative to [browaeys2004](@citet).

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
\begin{aligned}
&C_{ii}\ \text{unchanged},\\
&C_{ij}=\tfrac{1}{2}\bigl(C_{ij}+C_{ji}\bigr)\quad(i<j\le3),\\
&C_{44},C_{55},C_{66}\ \text{halved from the KM diagonal},
\end{aligned}
```

everything outside the two blocks being discarded. For an order-2 tensor the
same reasoning reduces to *keep the diagonal in the material frame*.

### Cubic symmetry

The cleanest case of all, because the three generators are not merely mutually
orthogonal but are **projectors** summing to the identity. The normal equations
degenerate into three traces (`_project_CUBIC_KM`):

```math
\alpha = \operatorname{tr}(\mathbb{C}\,\mathbb{J}),\qquad
\beta  = \tfrac12\operatorname{tr}(\mathbb{C}\,\mathbb{E}),\qquad
\gamma = \tfrac13\operatorname{tr}(\mathbb{C}\,\mathbb{T}),
```

the divisors being the traces of the projectors themselves. Concretely, ``\alpha`` is
the mean of the whole upper ``3\times3`` block, ``\beta`` the mean of its diagonal
minus its off-diagonal, ``\gamma`` the mean of the shear diagonal — so the projection
*averages* rather than discards, and a non-major-symmetric input is symmetrized
implicitly by the same traces.

At order 2 there is nothing to do: the octahedral group leaves no second-order
tensor invariant but a multiple of the identity, so the cubic projection **is**
the isotropic one. See [Cubic symmetry](@ref th-cubic).

With a **free** orientation the treatment is the orthotropic one unchanged: the
cube orientation is an ordinary rotation, so the objective
``j(\theta,\phi,\psi)`` is smooth and the same optimizer applies, started from
the same eigenstructure candidate. The discreteness of ``O_h`` shows up only as
a 24-fold degenerate minimum.

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

### The search

When the orientation is free, [`proj_tens`](@ref) searches it, with no
dependency beyond `ForwardDiff`:

1. the objective is **screened** on a fixed, quasi-uniform set of orientations:
   200 axes spread over a hemisphere along a golden-angle spiral for the TI
   axis; 100 such axes, each with four turns about it, for an orthotropic or a
   cubic frame — a quarter turn about an axis of the frame maps it onto an
   equivalent one. An evaluation costs about a third of a microsecond;
2. **Newton's method** on the angles refines the three (TI) or four (frames)
   best screened orientations, the eigenstructure candidate
   (`_candidate_TI_axis`, `_candidate_ORTHO_frame` — exact whenever the tensor
   genuinely has the symmetry sought; there is none for the cubic class, every
   second-order contraction of a cubic tensor being isotropic) and the
   canonical orientation. The Hessian is exact (`ForwardDiff`, a ``2\times2``
   or ``3\times3`` matrix), so the refinement converges quadratically near a
   minimum; away from one, each eigen-direction is scaled by the absolute value
   of its curvature and a backtracking line search keeps the objective
   decreasing;
3. the best objective over every refined orientation is kept.

Three properties follow:

- **reproducibility** — no random step, so repeated calls return bit-identical
  angles;
- **no regression** — the result is never worse than the projection about any
  refined orientation, the eigenstructure candidate and the canonical one
  included;
- **globality, as measured** — on 150 random anisotropic stiffnesses, the
  minimum found matched, to rounding, that of a dense search (3,000 axes; 7,200
  frames, the best of each refined), which a multi-start from a dozen fixed
  points missed for about one TI tensor in twenty. The test suite repeats the
  comparison on random tensors.

At **order 2** no search is needed. In the frame of the eigenvectors of the
symmetric part of ``\boldsymbol{A}``, the orthotropic projection keeps all of
that part, the most any frame can keep, so that frame is the optimum. About an
axis ``\underline{n}``, the transversely isotropic projection is
``a\,(\boldsymbol{1}-\underline{n}\otimes\underline{n})+b\,\underline{n}\otimes\underline{n}``
with ``b=q=\underline{n}\cdot\boldsymbol{A}\cdot\underline{n}`` and
``a=(\operatorname{tr}\boldsymbol{A}-q)/2``; its squared norm
``(\operatorname{tr}\boldsymbol{A}-q)^2/2+q^2`` is convex in ``q``, which ranges
over the eigenvalues of the symmetric part, so the optimum is at the smallest or
the largest one, and the axis is the corresponding eigenvector.

## Detecting the symmetry: the cascade

[`best_sym_tens`](@ref) tries the classes from the most restrictive to the least
and returns the first whose *relative* error falls below a tolerance ``\varepsilon``:

```math
\text{ISO}\ \longrightarrow\ \text{CUBIC}\ \longrightarrow\
\text{TI}\ \longrightarrow\ \text{ORTHO}
\ \longrightarrow\ \text{ANISO} .
```

The chain is ordered by **number of constants** — 2, 3, 5, 9 — and not by
inclusion, because it cannot be: ``\text{CUBIC}`` and ``\text{TI}`` are
incomparable, their intersection being ``\text{ISO}`` and neither containing
the other. A tensor satisfying both is reported cubic, having the fewer
constants; one satisfying neither falls through to ``\text{ORTHO}``, which
contains them both.

The relative criterion is what makes the tolerance dimensionless and independent
of the units of the moduli. Two modes:

| `optimize_angles` | TI axis / ORTHO frame / CUBIC frame |
| :---------------- | :-------------------- |
| `false` (default) | taken from the tensor if it is a structured container, otherwise from the Kelvin–Mandel eigenstructure |
| `true` | found by the search above |

The value-level predicates [`isisotropic`](@ref), [`istransverselyisotropic`](@ref), [`isorthotropic`](@ref),
[`iscubic`](@ref)
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
log-Euclidean, power-Euclidean, arctan-Euclidean — are constructed by
[morin2020](@citet). `TensND` implements the Euclidean one only, so a reported
projection should always state which of ``\mathbb{C}`` or ``\mathbb{S}`` was
projected.
