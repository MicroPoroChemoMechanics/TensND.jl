# Tensors

A tensor, parametrized by an order and a dimension, is in general defined by

- an array or a set of condensed parameters (e.g. isotropic tensors),
- a basis,
- a set of variances (covariant `:cov` or contravariant `:cont`) useful if the basis is not orthonormal.

In practice, the type of basis conditions the type of tensor (`TensCanonical`, `TensRotated`, `TensOrthogonal`, `Tens` or even `TensISO` in case of isotropic tensor).

```@repl tensors
using TensND, SymPy, Tensors
ℬ = Basis(Sym[0 1 1; 1 0 1; 1 1 0])
V = Tens(Tensor{1,3}(i -> symbols("v$i", real = true)))
components(V, ℬ, (:cont,))
components(V, ℬ, (:cov,))
ℬ̄ = normalize(ℬ)
components(V, ℬ̄, (:cov,))
T = Tens(Tensor{2,3}((i, j) -> symbols("t$i$j", real = true)))
components(T, ℬ, (:cov, :cov))
factor(simplify(components(T, ℬ, (:cont, :cov))))
```

Special tensors are available

- `tensId2(::Val{dim} = Val(3), ::Val{T} = Val(Sym)) where {dim,T<:Number}`: second-order identity (`𝟏ᵢⱼ = δᵢⱼ = 1 if i=j otherwise 0`)
- `tensId4(::Val{dim} = Val(3), ::Val{T} = Val(Sym)) where {dim,T<:Number}`: fourth-order identity with minor symmetries (`𝕀 = 𝟏 ⊠ˢ 𝟏` i.e. `(𝕀)ᵢⱼₖₗ = (δᵢₖδⱼₗ+δᵢₗδⱼₖ)/2`)
- `tensJ4(::Val{dim} = Val(3), ::Val{T} = Val(Sym)) where {dim,T<:Number}`: fourth-order spherical projector (`𝕁 = (𝟏 ⊗ 𝟏) / dim` i.e. `(𝕁)ᵢⱼₖₗ = δᵢⱼδₖₗ/dim`)
- `tensK4(::Val{dim} = Val(3), ::Val{T} = Val(Sym)) where {dim,T<:Number}`: fourth-order deviatoric projector (`𝕂 = 𝕀 - 𝕁` i.e. `(𝕂)ᵢⱼₖₗ = (δᵢₖδⱼₗ+δᵢₗδⱼₖ)/2 - δᵢⱼδₖₗ/dim`)
- `ISO(::Val{dim} = Val(3), ::Val{T} = Val(Sym)) where {dim,T<:Number}`: returns `𝕀, 𝕁, 𝕂`

The useful tensor products are the following:

- `⊗` tensor product
- `⊗ˢ` symmetrized tensor product
- `⊠` modified tensor product
- `⊠ˢ` symmetrized modified tensor product
- `⋅` contracted product
- `⊡` double contracted product
- `⊙` quadruple contracted product

```@repl tensors
𝟏 = tensId2(3, Sym)
𝕀, 𝕁, 𝕂 = ISO(3, Sym) ;
𝕀 == 𝟏 ⊠ˢ 𝟏
𝕁 == (𝟏 ⊗ 𝟏)/3
a = Tens(Vec{3}((i,) -> symbols("a$i", real = true))) ;
b = Tens(Vec{3}((i,) -> symbols("b$i", real = true))) ;
a ⊗ b
a ⊗ˢ b
```

The predefined spherical coordinate system `init_spherical()` provides the local orthonormal basis
``(\mathbf{e}_\theta, \mathbf{e}_\varphi, \mathbf{e}_r)`` in terms of polar angle ``\theta`` (from the ``z``-axis) and azimuthal angle ``\varphi``:

```math
\mathbf{e}_\theta = \cos\theta\cos\varphi\,\mathbf{e}_1 + \cos\theta\sin\varphi\,\mathbf{e}_2 - \sin\theta\,\mathbf{e}_3
```

```math
\mathbf{e}_\varphi = -\sin\varphi\,\mathbf{e}_1 + \cos\varphi\,\mathbf{e}_2
```

```math
\mathbf{e}_r = \sin\theta\cos\varphi\,\mathbf{e}_1 + \sin\theta\sin\varphi\,\mathbf{e}_2 + \cos\theta\,\mathbf{e}_3
```

The rotation matrix ``R = [\mathbf{e}_\theta \mid \mathbf{e}_\varphi \mid \mathbf{e}_r]`` encodes this change of basis.
For any vector ``\mathbf{A}`` in the canonical frame, `change_tens(A, ℬˢ)` returns its components in the spherical basis.
The example below verifies that if ``\mathbf{A} = R\,\mathbf{a}``, then expressing ``\mathbf{A}`` in ``\mathcal{B}^s`` recovers the original components ``(a_1, a_2, a_3)``:

```@repl tensors
(θ, ϕ, r), (𝐞ᶿ, 𝐞ᵠ, 𝐞ʳ), ℬˢ = init_spherical()
R = rot3(θ, ϕ)
A = Tens(R * a)
simplify(change_tens(A, ℬˢ))
```

## Isotropic tensors (TensISO)

Isotropic tensors are stored compactly: a second-order isotropic tensor ``\lambda\mathbf{1}`` is
parametrized by one scalar, while a fourth-order isotropic tensor ``\alpha\mathbb{J} + \beta\mathbb{K}``
is parametrized by two scalars.  All arithmetic operations (``+``, ``-``, ``\times``, ``\mathbb{A}:\mathbb{B}``,
``\mathbb{A}^{-1}``) exploit this compact form and remain in the `TensISO` type whenever possible.

The type predicates `isISO`, `isTI`, `isOrtho` allow querying the symmetry class of any tensor:

```@repl tensors_iso
using TensND, Tensors
𝟏 = tensId2(Val(3), Val(Float64))
𝕀, 𝕁, 𝕂 = ISO(Val(3), Val(Float64)) ;
isISO(𝕀)
isTI(𝕀)
isOrtho(𝕀)
```

The compact display reflects the algebraic form directly:

```@repl tensors_iso
show(stdout, 𝕁 + 𝕂)   # prints "(1.0) 𝕁 + (1.0) 𝕂"
show(stdout, 2.0 * 𝟏)  # prints "(2.0) 𝟏"
```

## Transverse isotropy and orthotropy

### TensWalpole

A transversely isotropic 4th-order tensor with symmetry axis ``\mathbf{n}`` is decomposed in the Walpole basis:

```math
L = \ell_1 W_1 + \ell_2 W_2 + \ell_3 W_3 + \ell_4 W_4 + \ell_5 W_5 + \ell_6 W_6
```

where ``\mathbf{n}_n = \mathbf{n}\otimes\mathbf{n}``, ``\mathbf{n}_T = \mathbf{1} - \mathbf{n}_n`` and

| Tensor | Expression |
| ------ | ----------- |
| ``W_1`` | ``\mathbf{n}_n\otimes\mathbf{n}_n`` |
| ``W_2`` | ``(\mathbf{n}_T\otimes\mathbf{n}_T)/2`` |
| ``W_3`` | ``(\mathbf{n}_n\otimes\mathbf{n}_T)/\sqrt{2}`` |
| ``W_4`` | ``(\mathbf{n}_T\otimes\mathbf{n}_n)/\sqrt{2}`` |
| ``W_5`` | ``\mathbf{n}_T\,\overline{\boxtimes}^s\,\mathbf{n}_T - (\mathbf{n}_T\otimes\mathbf{n}_T)/2`` |
| ``W_6`` | ``\mathbf{n}_T\,\overline{\boxtimes}^s\,\mathbf{n}_n + \mathbf{n}_n\,\overline{\boxtimes}^s\,\mathbf{n}_T`` |

The double contraction follows the **synthetic Walpole rule**:

```math
L\colon M \equiv \left(\begin{bmatrix}\ell_1 & \ell_3\\\ell_4 & \ell_2\end{bmatrix}\begin{bmatrix}m_1 & m_3\\m_4 & m_2\end{bmatrix},\; \ell_5 m_5,\; \ell_6 m_6\right)
```

For major-symmetric tensors (``\ell_3=\ell_4``), use `N=5`; for general tensors, `N=6`.

The `show` method displays the tensor in its compact Walpole form, including the symmetry axis:

```@repl tensors_walpole
using TensND, Tensors
n = 𝐞(3) ;
W1, W2, W3, W4, W5, W6 = Walpole(n) ;
L = TensWalpole(2., 1., 0.5, 0.3, 0.8, n)
show(stdout, L)
maximum(abs.(getarray(L ⊡ inv(L)) - getarray(tensId4(Val(3), Val(Float64)))))
𝕀, 𝕁, 𝕂 = ISO() ; L2 = fromISO(3𝕁 + 2𝕂, n)
isTI(L)
isISO(L)
isOrtho(L)
```

An isotropic tensor converted to `TensWalpole` via `fromISO` retains the `isTI` predicate, and
symbolic manipulations via `tsimplify`, `tsubs`, `tdiff`, etc. preserve the `TensWalpole` type:

```@repl tensors_walpole
using SymPy
ℓ₁, ℓ₂, ℓ₃ = symbols("ℓ₁ ℓ₂ ℓ₃", real = true) ;
ns = 𝐞(Val(3), Val(3), Val(Sym)) ;
Ls = TensWalpole(ℓ₁, ℓ₂, ℓ₃, ℓ₁ + ℓ₂, ℓ₂ + ℓ₃, ns) ;
Ls_simp = tsimplify(Ls) ;
Ls_simp isa TensWalpole
```

### TensOrtho

An orthotropic 4th-order tensor in material frame ``(\mathbf{e}_1,\mathbf{e}_2,\mathbf{e}_3)``
with ``P_m = \mathbf{e}_m\otimes\mathbf{e}_m`` has 9 independent elastic constants:

```math
\mathbb{C} = C_{11}P_1{\otimes}P_1 + C_{22}P_2{\otimes}P_2 + C_{33}P_3{\otimes}P_3
+ C_{12}(P_1{\otimes}P_2+P_2{\otimes}P_1) + C_{13}(P_1{\otimes}P_3+P_3{\otimes}P_1) + C_{23}(P_2{\otimes}P_3+P_3{\otimes}P_2)
+ 2C_{44}(P_2\,\overline{\boxtimes}^s P_3) + 2C_{55}(P_1\,\overline{\boxtimes}^s P_3) + 2C_{66}(P_1\,\overline{\boxtimes}^s P_2)
```

The Kelvin-Mandel matrix in the material frame (ordering ``11,22,33,23,13,12``) is block-diagonal.
Use `KM_material(t)` to retrieve it; `KM(t)` gives the matrix in the canonical frame.

The `show` method displays all 9 constants and the material frame, and `isOrtho` identifies the type:

```@repl tensors_ortho
using TensND, Tensors
ℬ = CanonicalBasis{3,Float64}() ;
t = TensOrtho(10., 8., 9., 3., 2., 4., 2.5, 3., 1.5, ℬ) ;
show(stdout, t)
KM_material(t)
maximum(abs.(getarray(t) ⊡ getarray(inv(t)) - getarray(tensId4(Val(3), Val(Float64)))))
isOrtho(t)
isTI(t)
isISO(t)
```

### TensTI (2nd-order transversely isotropic)

A 2nd-order transversely isotropic tensor with symmetry axis ``\mathbf{n}`` is decomposed as:

```math
\mathbf{A} = a\,\mathbf{n}_T + b\,\mathbf{n}_n
```

where ``\mathbf{n}_n = \mathbf{n}\otimes\mathbf{n}`` and ``\mathbf{n}_T = \mathbf{1} - \mathbf{n}_n``.
The scalar ``a`` is the transverse coefficient (in the plane ``\perp\mathbf{n}``) and ``b`` is the
axial coefficient (along ``\mathbf{n}``).

The type `TensTI{order,T,N}` mirrors the parametric design of `TensISO{order,dim,T,N}`.
For order 2, `N=2` and `data = (a, b)`.  When `a = b`, the tensor is isotropic.

```@repl tensors_ti
using TensND, LinearAlgebra
n = [0., 0., 1.] ;
A = TensTI{2}(5.0, 8.0, n)
getarray(A)
tr(A)
isISO(A)
isTI(A)
inv(A).data
```

`TensTI{2}` supports all standard operations: addition, subtraction, scalar multiplication,
inversion.  Two TI tensors with the **same axis** can be combined:

```@repl tensors_ti
B = TensTI{2}(3.0, 2.0, n) ;
(A + B).data
(A - B).data
(2.0 * A).data
```

When the axis is not aligned with ``\mathbf{e}_3``, the full 3×3 matrix correctly reflects the
off-diagonal structure:

```@repl tensors_ti
n45 = [1/√2, 0., 1/√2] ;
C = TensTI{2}(3.0, 7.0, n45) ;
getarray(C)
```

### Symmetry class predicates

The three predicates `isISO`, `isTI`, `isOrtho` form a consistent hierarchy across all specialized
tensor types.  Any value that is not a recognized tensor type returns `false` for all three:

```@repl tensors_preds
using TensND, Tensors
𝕀, 𝕁, 𝕂 = ISO(Val(3), Val(Float64)) ;
n = 𝐞(3) ;
L = TensWalpole(2., 1., 0.5, 3., 4., n) ;
A2 = TensTI{2}(5.0, 8.0, n) ;
ℬ = CanonicalBasis{3,Float64}() ;
t = TensOrtho(10., 8., 9., 3., 2., 4., 2.5, 3., 1.5, ℬ) ;
(isISO(𝕀), isTI(𝕀),  isOrtho(𝕀))
(isISO(L),  isTI(L),  isOrtho(L))
(isISO(A2), isTI(A2), isOrtho(A2))
(isISO(t),  isTI(t),  isOrtho(t))
```

## Projection onto symmetry subspaces

Given an arbitrary tensor (2nd or 4th order), one can project it onto the closest tensor
with a prescribed symmetry class (ISO, TI, ORTHO).  The projection minimises the Frobenius
distance ``\lVert B - A \rVert`` over all tensors ``B`` of the target symmetry.

The function `proj_tens` provides this projection.  It returns a 3-tuple `(B, d, drel)`:

- `B`: the projected tensor
- `d`: absolute Frobenius distance ``\lVert B - A \rVert``
- `drel`: relative distance ``d / \lVert A \rVert``

### Fixed-axis TI projection (order 4)

Project a 4th-order tensor onto the TI subspace with a given axis ``\mathbf{n}``.
The result is a `TensWalpole{T,5}`:

```@repl tensors_proj
using TensND, LinearAlgebra
n = [0., 0., 1.] ;
C = tensTI(10., 3., 2.5, 12., 2., n) ;
B, d, drel = proj_tens(:TI, getarray(C), n) ;
drel < 1e-12
B isa TensWalpole
```

### Fixed-axis TI projection (order 2)

For a 2nd-order tensor, the projection returns a `TensTI{2}`:

```@repl tensors_proj
A = Float64[5 1 0; 1 5 0; 0 0 8] ;
B2, d2, drel2 = proj_tens(:TI, A, n) ;
B2
B2.data
```

### Fixed-frame ORTHO projection (order 4)

Project onto the orthotropic subspace with a given material frame:

```@repl tensors_proj
frame = CanonicalBasis{3,Float64}() ;
t = TensOrtho(10., 8., 12., 3., 2.5, 1.5, 2., 3., 3.5, frame) ;
Bo, do_, drelo = proj_tens(:ORTHO, getarray(t), frame) ;
drelo < 1e-12
Bo isa TensOrtho
```

### Fixed-frame ORTHO projection (order 2)

```@repl tensors_proj
M = Float64[5 1 2; 1 8 3; 2 3 12] ;
Bm, dm, drelm = proj_tens(:ORTHO, M, frame) ;
Bm
dm > 0
```

### Best symmetry detection

The function `best_sym_tens` tries symmetries from the most restrictive to the least
(ISO → TI → ORTHO) and returns the first whose relative projection error falls below
a threshold ``\varepsilon`` (default `1e-6`).  Pass an axis or frame for fixed-basis detection:

```@repl tensors_proj
C_ti = tensTI(10., 3., 2.5, 12., 2., n) ;
_, _, _, sym = best_sym_tens(C_ti, n; proj = (:TI, :ORTHO)) ;
sym
```

### Rotation-optimized projection (requires NLopt)

When no axis or frame is provided, `proj_tens` optimises the rotation angles to find the
best approximation.  This requires the `NLopt` package (declared as a weak dependency):

```julia
using NLopt   # triggers the TensNDNLoptExt extension

# Optimise the TI axis automatically
B_opt, d_opt, drel_opt = proj_tens(:TI, getarray(C))

# Optimise the ORTHO frame automatically
B_ort, d_ort, drel_ort = proj_tens(:ORTHO, getarray(C))
```

The optimiser uses a two-pass strategy (global + local refinement) inspired by the
ECHOES library, with gradients computed via ForwardDiff.
