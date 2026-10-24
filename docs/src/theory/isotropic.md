# [Isotropic tensors](@id th-isotropic)

The most constrained symmetry class, and the template for every other one: a
small closed algebra in which products and inverses are computed on a couple of
scalars rather than on ``d^4`` components. Implemented by
[`TensISO`](@ref) in `src/tens_isotropic.jl`, in **arbitrary dimension** ``d``.

## The two projectors

With ``\boldsymbol{1}`` the order-2 identity and
``\mathbb{I}=\boldsymbol{1}\stackrel{s}{\boxtimes}\boldsymbol{1}`` the order-4
identity on symmetric tensors,

```math
\mathbb{J}=\frac{1}{d}\,\boldsymbol{1}\otimes\boldsymbol{1},
\qquad
\mathbb{K}=\mathbb{I}-\mathbb{J} .
```

``\mathbb{J}`` extracts the spherical part of a symmetric tensor and
``\mathbb{K}`` the deviatoric one. They are **complementary orthogonal
projectors**:

```math
\mathbb{J}:\mathbb{J}=\mathbb{J},
\qquad
\mathbb{K}:\mathbb{K}=\mathbb{K},
\qquad
\mathbb{J}:\mathbb{K}=\mathbb{K}:\mathbb{J}=0,
\qquad
\mathbb{J}+\mathbb{K}=\mathbb{I}.
```

Their Frobenius norms are the dimensions of the subspaces they project onto:

```math
\mathbb{J}::\mathbb{J}=1,
\qquad
\mathbb{K}::\mathbb{K}=\frac{d(d+1)}{2}-1
\quad(=5\text{ in 3-D}),
\qquad
\mathbb{J}::\mathbb{K}=0 .
```

In `TensND` the triple is [`iso_projectors`](@ref) / [`iso_projectors`](@ref), and
individually [`tens_J4`](@ref), [`tens_K4`](@ref), [`tens_Id4`](@ref),
[`tens_Id2`](@ref).

## The algebra

Every isotropic minor-symmetric order-4 tensor is

```math
\mathbb{A}=\alpha\,\mathbb{J}+\beta\,\mathbb{K},
```

and since ``\mathbb{J}`` and ``\mathbb{K}`` are orthogonal idempotents the
algebra is that of a pair of independent scalars:

```math
\begin{aligned}
\mathbb{A}:\mathbb{B}&=\alpha\alpha'\,\mathbb{J}+\beta\beta'\,\mathbb{K},\\[4pt]
\mathbb{A}^{-1}&=\frac{1}{\alpha}\,\mathbb{J}+\frac{1}{\beta}\,\mathbb{K},\\[4pt]
\mathbb{A}::\mathbb{B}&=\alpha\alpha'+\Bigl(\frac{d(d+1)}{2}-1\Bigr)\beta\beta' .
\end{aligned}
```

In elasticity ``\alpha=3k`` and ``\beta=2\mu`` in 3-D. `TensISO{4}` stores
exactly the pair ``(\alpha,\beta)``; `TensISO{2}` stores the single scalar of
``\lambda\boldsymbol{1}``.

## Isotropization

The **closest isotropic tensor** to an arbitrary ``\mathbb{T}``, for the
Frobenius distance, is its orthogonal projection onto
``\mathrm{span}(\mathbb{J},\mathbb{K})``. Because the two projectors are
orthogonal with known norms, the normal equations are diagonal and the answer is
closed-form [bornert2001](@cite):

```math
\mathrm{ISO}(\mathbb{T})
=\bigl(\mathbb{T}::\mathbb{J}\bigr)\,\mathbb{J}
+\frac{\mathbb{T}::\mathbb{K}}{5}\,\mathbb{K}
\qquad(d=3),
```

the general-``d`` denominator being ``\mathbb{K}::\mathbb{K}``. This is
[`isotropify`](@ref), and it is the simplest instance of the general machinery
of [Projection onto a symmetry class](@ref th-projection): no orientation to
optimize, because an isotropic tensor has none.

!!! warning "Isotropization does not commute with inversion"
    ```math
    \mathrm{ISO}\bigl(\mathbb{T}^{-1}\bigr)\neq\bigl(\mathrm{ISO}(\mathbb{T})\bigr)^{-1}
    ```

    Projecting a stiffness and projecting the corresponding compliance give
    *different* isotropic materials. This is not a defect of the
    implementation but of the Euclidean distance itself, which is not invariant
    under inversion. Distances that are — log-Euclidean, power-Euclidean,
    arctan-Euclidean — are constructed and compared by [morin2020](@citet);
    `TensND` implements the Euclidean one, so the choice of which of
    ``\mathbb{C}`` or ``\mathbb{S}`` to project is the user's and must be
    stated when a result is reported.

## Storage

| Storage | Coefficients | Case |
| :------ | :----------- | :--- |
| `TensISO{2,dim,T,1}` | ``\lambda`` | ``\lambda\boldsymbol{1}`` |
| `TensISO{4,dim,T,2}` | ``(\alpha,\beta)`` | ``\alpha\mathbb{J}+\beta\mathbb{K}`` |

One shape per order, invariant under the whole orthogonal group: an isotropic
tensor is its own transpose, and the products and inverses above stay in the
shape, so there is no wider one as there is for the transversely isotropic and
orthotropic classes.

## Relation to the wider classes

Isotropy is the intersection of every other class, so an isotropic tensor
satisfies all three predicates:

```math
\text{ISO}\subset\text{TI}\subset\text{ORTHO}\subset\text{ANISO} .
```

[`isisotropic`](@ref), [`istransverselyisotropic`](@ref) and [`isorthotropic`](@ref) reflect this — all three
return `true` on ``\mathbb{I}``. Conversion up the chain is explicit:
[`iso_to_ti`](@ref) re-expresses an isotropic tensor on the Walpole basis about a
chosen axis, and [`iso_to_ortho`](@ref) on an orthotropic frame; the coefficients
that result are computed on [Walpole basis](@ref th-walpole) and
[Orthotropy](@ref th-orthotropy).
