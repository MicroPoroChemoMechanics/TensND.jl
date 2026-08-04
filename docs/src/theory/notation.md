# [Notation and conventions](@id th-notation)

This page fixes the notation used throughout the documentation. Every symbol
that is *not* listed here is redefined on the page where it appears.

## Tensor order is carried by the typeface

The convention is that of the
[Echoes manual](https://jfbarthelemy.github.io/echoes/), so that formulas can be
compared side by side with it.

| Object | Typeset as | Example |
| :----- | :--------- | :------ |
| scalar | italic | ``\lambda``, ``\mu``, ``\chi_i``, ``\ell_i`` |
| vector (order 1) | underlined | ``\underline{u}``, ``\underline{n}``, ``\underline{e}_i`` |
| tensor of order 2 | bold | ``\boldsymbol{a}``, ``\boldsymbol{1}``, ``\boldsymbol{g}`` |
| tensor of order 4 | blackboard bold | ``\mathbb{C}``, ``\mathbb{I}``, ``\mathbb{J}``, ``\mathbb{W}_i`` |
| tensor of arbitrary order | calligraphic | ``\mathcal{T}``, ``\mathcal{T}'`` |
| basis, coordinate system | script | ``\mathcal{B}``, ``\mathcal{S}`` |

## Index conventions

Einstein summation over repeated indices is implicit. Because `TensND` supports
**non-orthonormal bases**, the position of an index is meaningful: a subscript
is a covariant component, a superscript a contravariant one. See
[Bases and variance](@ref th-bases-variance).

| Symbol | Meaning |
| :----- | :------ |
| ``\delta_{ij}``, ``\delta^i_j`` | Kronecker symbol |
| ``\underline{e}_i`` / ``\underline{e}^i`` | basis vector / dual (reciprocal) basis vector |
| ``\underline{a}_i`` / ``\underline{a}^i`` | natural / dual natural vector of a coordinate system |
| ``g_{ij}`` / ``g^{ij}`` | covariant / contravariant metric components |
| ``\varepsilon_{ijk}`` | Levi-Civita symbol |

## Products

The full definitions, with their index formulas and algebraic identities, are on
[Tensor algebra](@ref th-tensor-algebra). This table is the dictionary between
mathematical symbol, Julia operator and function name.

| Product | Julia | Function | Acts on |
| :------ | :---- | :------- | :------ |
| ``\mathcal{T}\otimes\mathcal{T}'`` | `⊗` | `otimes` | any orders |
| ``\mathcal{T}\stackrel{s}{\otimes}\mathcal{T}'`` | `⊗ˢ` | `sotimes` | any orders |
| ``\boldsymbol{a}\boxtimes\boldsymbol{b}`` | `⊠` | `otimesu` | any orders |
| ``\boldsymbol{a}\stackrel{s}{\boxtimes}\boldsymbol{b}`` | `⊠ˢ` | `otimesul` (alias `sboxtimes`) | any orders |
| ``\mathcal{T}\cdot\mathcal{T}'`` | `⋅` | `dot` | 1 contracted index |
| ``\mathcal{T}:\mathcal{T}'`` | `⊡` | `dcontract` | 2 contracted indices |
| ``\mathcal{T}::\mathcal{T}'`` | `⊙` | `qcontract` | 4 contracted indices |

!!! warning "Double contraction: the pair-wise convention"
    Two conventions coexist for the double dot product. The classical one
    consumes indices going *inwards from the extremities*. The one adopted here
    — and by the Echoes manual — treats the **two last indices of the left
    operand as a pair** and the **two first indices of the right operand** as
    the corresponding pair:

    ```math
    \boldsymbol{a}:\boldsymbol{b}=a_{ij}b_{ij},
    \qquad
    \mathbb{T}:\boldsymbol{a}=T_{ijkl}\,a_{kl}\;
    \underline{e}_i\otimes\underline{e}_j .
    ```

    The two conventions differ by a transpose on the right operand and agree
    only when it is symmetric. `TensND` implements the pair-wise one
    (`dcontract` in `src/array_utils.jl`).

## Two rules this documentation follows

**Every symbol is defined where it is used.** A page never relies on a symbol
introduced only on another page, even at the cost of repeating a definition.

!!! warning "No formula without a traceable source"
    Every expression in this documentation is either (i) accompanied by a
    citation to published work, (ii) derived explicitly on the page from
    expressions that are, or (iii) read off the implementation, with the source
    file named. Where a convention differs between references, the competing
    conventions are named and the one implemented by `TensND` is stated.

    Formulas that a reader is likely to want to check — the Walpole
    multiplication table, the Kelvin–Mandel congruence, the operator index
    placement — are additionally **pinned by tests**, so that the documentation
    cannot silently drift from the code. See
    [Testing and conventions](@ref dev-testing).
