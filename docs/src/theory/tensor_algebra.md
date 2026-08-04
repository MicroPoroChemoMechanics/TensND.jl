# [Tensor algebra](@id th-tensor-algebra)

Definitions and identities of the products `TensND` implements. Index formulas
are those of the OMEinsum contraction codes in `src/array_utils.jl`; the
algebraic identities follow the
[Echoes manual](https://jfbarthelemy.github.io/echoes/) appendix on tensor
algebra.

Components below are written in an orthonormal frame ``(\underline{e}_i)`` so
that index position carries no extra information; the general, variance-aware
statement is on [Bases and variance](@ref th-bases-variance).

## Tensor products

For ``\mathcal{T}`` of order ``p`` and ``\mathcal{T}'`` of order ``q``, the
tensor product is the order-``(p+q)`` tensor

```math
\mathcal{T}\otimes\mathcal{T}'
= \mathcal{T}_{i_1\ldots i_p}\,\mathcal{T}'_{i_{p+1}\ldots i_{p+q}}\;
\underline{e}_{i_1}\otimes\ldots\otimes\underline{e}_{i_{p+q}} .
```

The **symmetrized** tensor product ``\stackrel{s}{\otimes}`` symmetrizes over
the last index of ``\mathcal{T}`` and the first of ``\mathcal{T}'``:

```math
\mathcal{T}\stackrel{s}{\otimes}\mathcal{T}'
= \frac{\mathcal{T}_{i_1\ldots i_p}\mathcal{T}'_{i_{p+1}\ldots i_{p+q}}
      + \mathcal{T}_{i_1\ldots i_{p-1}i_{p+1}}\mathcal{T}'_{i_p i_{p+2}\ldots i_{p+q}}}{2}\;
\underline{e}_{i_1}\otimes\ldots\otimes\underline{e}_{i_{p+q}} ,
```

so that for two vectors

```math
\underline{u}\stackrel{s}{\otimes}\underline{v}
= \frac{\underline{u}\otimes\underline{v}+\underline{v}\otimes\underline{u}}{2} .
```

## Box products

Two order-2 tensors combine into an order-4 tensor in two further ways, which
are what make the order-4 identity and the deviatoric projector expressible in
closed form:

```math
(\boldsymbol{a}\boxtimes\boldsymbol{b})_{ijkl}=a_{ik}\,b_{jl},
\qquad
(\boldsymbol{a}\stackrel{s}{\boxtimes}\boldsymbol{b})_{ijkl}
=\frac{a_{ik}b_{jl}+a_{il}b_{jk}}{2} .
```

`TensND` also exposes the intermediate `otimesl`, with
``(\boldsymbol{a}\boxtimes_{\!l}\boldsymbol{b})_{ijkl}=a_{il}b_{jk}``, from
which ``\stackrel{s}{\boxtimes}`` is the half-sum.

With ``\boldsymbol{1}`` the order-2 identity, the two order-4 identities are

```math
\mathbb{1}=\boldsymbol{1}\boxtimes\boldsymbol{1}
\quad(\mathbb{1}_{ijkl}=\delta_{ik}\delta_{jl}),
\qquad
\mathbb{I}=\boldsymbol{1}\stackrel{s}{\boxtimes}\boldsymbol{1}
\quad\Bigl(\mathbb{I}_{ijkl}=\tfrac{\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}}{2}\Bigr).
```

``\mathbb{1}`` is the identity of the space of *all* order-2 tensors,
``\mathbb{I}`` that of the **symmetric** ones — which is the one
[`tens_Id4`](@ref) returns, since every symmetry class in this library lives in
the minor-symmetric subspace.

## Contractions

| Operation | Definition | Resulting order |
| :-------- | :--------- | :-------------- |
| ``\mathcal{T}\cdot\mathcal{T}'`` | ``\mathcal{T}_{i_1\ldots i_{p-1}k}\,\mathcal{T}'_{k\,i_{p}\ldots}`` | ``p+q-2`` |
| ``\mathcal{T}:\mathcal{T}'`` | ``\mathcal{T}_{i_1\ldots i_{p-2}kl}\,\mathcal{T}'_{kl\,i_{p-1}\ldots}`` | ``p+q-4`` |
| ``\mathcal{T}::\mathcal{T}'`` | ``\mathcal{T}_{i_1\ldots i_{p-4}klmn}\,\mathcal{T}'_{klmn\ldots}`` | ``p+q-8`` |

Each consumes the *last* indices of the left operand against the *first* of the
right one — see the warning on [Notation](@ref th-notation) about the competing
double-contraction convention.

Between two order-4 tensors, ``::`` is the **Frobenius scalar product**

```math
\mathbb{T}::\mathbb{T}'=T_{ijkl}T'_{ijkl}
=\langle\mathbb{T},\mathbb{T}'\rangle ,
\qquad
\|\mathbb{T}\|^2=\mathbb{T}::\mathbb{T} ,
```

and it is this scalar product that every projection in
[Projection onto a symmetry class](@ref th-projection) minimizes against.

## Transpose of an order-4 tensor

Consistently with the pair-wise double contraction, the transpose exchanges the
two index pairs:

```math
{}^{t}\mathbb{T}:\boldsymbol{a}=\boldsymbol{a}:\mathbb{T}
\quad\Longleftrightarrow\quad
({}^{t}\mathbb{T})_{ijkl}=\mathbb{T}_{klij} .
```

A tensor with ``\mathbb{T}={}^{t}\mathbb{T}`` is **major-symmetric**; one with
``T_{ijkl}=T_{jikl}=T_{ijlk}`` is **minor-symmetric**. The two are independent,
and the distinction decides how many coefficients a symmetry class needs — five
against six for transverse isotropy ([Walpole basis](@ref th-walpole)), nine
against twelve for orthotropy ([Orthotropy](@ref th-orthotropy)).

## Identities

The following hold for order-2 ``\boldsymbol{a},\boldsymbol{b},\boldsymbol{c},\boldsymbol{d}``
and are what allow the structured types to compute products and inverses in
closed form rather than by expanding 81 components:

```math
\begin{aligned}
&(\boldsymbol{a}\boxtimes\boldsymbol{b}):(\boldsymbol{c}\boxtimes\boldsymbol{d})
 =(\boldsymbol{a}\cdot\boldsymbol{c})\boxtimes(\boldsymbol{b}\cdot\boldsymbol{d})
 &&\text{(a)}\\[2pt]
&(\boldsymbol{a}\boxtimes\boldsymbol{b}):\boldsymbol{c}
 =\boldsymbol{a}\cdot\boldsymbol{c}\cdot{}^{t}\boldsymbol{b}
 &&\text{(b)}\\[2pt]
&(\boldsymbol{a}\stackrel{s}{\boxtimes}\boldsymbol{b}):\boldsymbol{c}
 =\frac{\boldsymbol{a}\cdot\boldsymbol{c}\cdot{}^{t}\boldsymbol{b}
       +\boldsymbol{a}\cdot{}^{t}\boldsymbol{c}\cdot{}^{t}\boldsymbol{b}}{2}
 &&\text{(c)}\\[2pt]
&(\boldsymbol{a}\otimes\boldsymbol{b}):(\boldsymbol{c}\otimes\boldsymbol{d})
 =(\boldsymbol{b}:\boldsymbol{c})\;\boldsymbol{a}\otimes\boldsymbol{d}
 &&\text{(d)}\\[2pt]
&(\boldsymbol{a}\boxtimes\boldsymbol{b})^{-1}
 =\boldsymbol{a}^{-1}\boxtimes\boldsymbol{b}^{-1}
 &&\text{(e)}
\end{aligned}
```

!!! warning "The symmetrized box product does not invert termwise"
    Identity (e) has **no counterpart** for ``\stackrel{s}{\boxtimes}``:

    ```math
    (\boldsymbol{a}\stackrel{s}{\boxtimes}\boldsymbol{b})^{-1}
    \neq\boldsymbol{a}^{-1}\stackrel{s}{\boxtimes}\boldsymbol{b}^{-1}
    \quad\text{in general.}
    ```

    Equality requires ``\boldsymbol{a}`` and ``\boldsymbol{b}`` to be
    **proportional**, ``\boldsymbol{b}=\alpha\,\boldsymbol{a}``, in which case
    the scalar cancels between the two factors. Commuting is *not* enough: two
    distinct diagonal tensors commute, yet their symmetrized box product still
    fails to invert termwise. Nor does taking ``\boldsymbol{b}=\boldsymbol{1}``
    help.

    This is why inversion is implemented per symmetry class — closed forms on
    the ``(\alpha,\beta)`` pair, on the Walpole triplet, or on the orthotropic
    block — rather than through one generic formula.

Note also that ``\mathbb{I}`` inverts and acts as an identity only on
**symmetric** order-2 tensors: ``\mathbb{I}:\boldsymbol{a}`` is the symmetric
part of ``\boldsymbol{a}``, whereas ``\mathbb{1}:\boldsymbol{a}=\boldsymbol{a}``
for every ``\boldsymbol{a}``.

## Symbolic and numerical evaluation

All of the above is implemented once, generically, on `AbstractArray`s through
OMEinsum contraction codes, and therefore applies unchanged to `Float64`,
`ForwardDiff.Dual`, `SymPy.Sym` and `Symbolics.Num` element types. The
structured types ([`TensISO`](@ref), [`TensTI`](@ref), [`TensOrtho`](@ref))
override the generic route with closed forms whenever the result stays in the
class — see [Structured tensors](@ref man-structured).
