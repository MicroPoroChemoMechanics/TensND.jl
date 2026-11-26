# [Tensors of the plane, of space and of any dimension](@id th-dimensions)

!!! tip "See also"
    **Manual**: [Two dimensions, and any dimension](@ref man-dimensions)

    **Tutorials**: [The symmetry classes of a plane section](@ref tut-plane-classes) · [Elastic Green's functions](@ref tut-green-function)

Most of the theory pages are written in dimension 3. This page says what
depends on the dimension ``d``: the Kelvin–Mandel notation, the isotropic
class and its moduli, then, for ``d = 2``, the symmetry classes of the plane,
how a tensor of space is restricted to one of its planes, and how a tensor of
the plane is projected onto its classes. The calls are in the Manual,
[Two dimensions, and any dimension](@ref man-dimensions).

## What the dimension changes

A tensor of order ``p`` in dimension ``d`` has ``d^p`` components. A symmetric
tensor of order 2 has ``n = d(d+1)/2`` independent ones, and a tensor of order
4 with the minor symmetries is a linear map on these ``n`` components:

| | ``d = 2`` | ``d = 3`` | ``d`` |
| :--- | :---: | :---: | :---: |
| components of a symmetric tensor of order 2 | 3 | 6 | ``d(d+1)/2`` |
| constants of a stiffness (major and minor symmetries) | 6 | 21 | ``n(n+1)/2`` |
| isotropic stiffnesses | 2 | 2 | 2 |
| parameters of a rotation | 1 | 3 | ``d(d-1)/2`` |

The algebra of the tensors (products, contractions, changes of basis) and the
isotropic class are written for any ``d``. The classes with a frame, the
harmonics and the curvilinear systems are those of the plane and of space.

## Kelvin–Mandel notation in dimension ``d``

The components of a symmetric tensor ``\boldsymbol{\varepsilon}`` are listed
with the diagonal first, then the off-diagonal pairs from the last column back,
each off-diagonal one weighted by ``\sqrt2``:

```math
\underline{\varepsilon} = \bigl(\varepsilon_{11}, \dots, \varepsilon_{dd},\;
\sqrt2\,\varepsilon_{d-1,d}, \dots, \sqrt2\,\varepsilon_{12}\bigr) ,
\qquad
\underline{\varepsilon}\cdot\underline{\varepsilon}' = \boldsymbol{\varepsilon}:\boldsymbol{\varepsilon}' .
```

The weights make the list an orthonormal system of components: a double
contraction is a scalar product, and a tensor ``\mathbb{C}`` of order 4 with the
minor symmetries is the ``n\times n`` matrix ``C_{IJ} = w_I w_J\,C_{ijkl}``,
with ``I = (i,j)``, ``J = (k,l)`` and ``w = 1`` or ``\sqrt2``. Products and
inverses of tensors are those of their matrices
([Kelvin–Mandel notation](@ref th-kelvin-mandel)).

```@example dims
using TensND, LinearAlgebra
mandel_pairs(Val(2)), mandel_pairs(Val(3))
```

```@example dims
mandel_length(4), mandel_pairs(Val(4))
```

A rotation ``\boldsymbol{R}`` of the basis, whose columns are the new vectors,
acts on the Kelvin–Mandel components by an orthogonal matrix
``\boldsymbol{Q}`` of size ``n``,

```math
Q_{IJ} = \frac{w_I w_J}{2}\,\bigl(R_{ik}R_{jl} + R_{il}R_{jk}\bigr) ,
\qquad \mathbb{C}' = \boldsymbol{Q}^{\mathsf T}\,\mathbb{C}\,\boldsymbol{Q} ,
```

which is `mandel_rotation`. In the plane it depends on one angle:

```@example dims
using SymPyPythonCall
θ = symbols("θ", real = true)
mandel_rotation(θ)
```

## Isotropic tensors and their moduli in dimension ``d``

The identity ``\mathbb{I}`` of the symmetric tensors splits into the two
projectors

```math
\mathbb{J} = \frac{1}{d}\,\boldsymbol{1}\otimes\boldsymbol{1} ,
\qquad \mathbb{K} = \mathbb{I} - \mathbb{J} ,
```

on the spherical tensors and on the deviators, of traces ``1`` and ``n - 1``,
and an isotropic stiffness is

```math
\mathbb{C} = d\,k\,\mathbb{J} + 2\mu\,\mathbb{K}
= \lambda\,\boldsymbol{1}\otimes\boldsymbol{1} + 2\mu\,\mathbb{I} ,
\qquad \lambda = k - \frac{2\mu}{d} .
```

``k`` relates the traces, ``\operatorname{tr}\boldsymbol{\sigma} = d\,k\,
\operatorname{tr}\boldsymbol{\varepsilon}``, and ``\mu`` the deviators. Young's
modulus and Poisson's ratio of the medium of dimension ``d`` are read on a
uniaxial stress ``\sigma_{11}``: with the compliance
``\mathbb{S} = \mathbb{J}/(d\,k) + \mathbb{K}/(2\mu)``,

```math
\varepsilon_{11} = \Bigl(\frac{1}{d^2 k} + \frac{d-1}{2\mu\,d}\Bigr)\sigma_{11} ,
\qquad
\varepsilon_{22} = \Bigl(\frac{1}{d^2 k} - \frac{1}{2\mu\,d}\Bigr)\sigma_{11} ,
```

so that, with ``E_d = \sigma_{11}/\varepsilon_{11}`` and
``\nu_d = -\varepsilon_{22}/\varepsilon_{11}``,

```math
\nu_d = \frac{d\,k - 2\mu}{d(d-1)\,k + 2\mu} = \frac{\lambda}{(d-1)\lambda + 2\mu} ,
\qquad
E_d = 2\mu\,(1 + \nu_d) = \frac{2\mu\,d^2 k}{d(d-1)\,k + 2\mu} .
```

For ``d = 3`` these are the usual ones. The builders and the readers take the
dimension as `Val(d)`:

```@example dims
k, μ = symbols("k μ", positive = true)
ℂ₄ = iso_stiffness(Val(4), k, μ)
simplify.(young_poisson(ℂ₄))
```

### The two pairs ``(E, \nu)`` of a tensor of the plane

For ``d = 2`` the formulas above give the moduli of the **plane medium**,

```math
E_2 = \frac{4k\mu}{k + \mu} , \qquad \nu_2 = \frac{k - \mu}{k + \mu} .
```

A tensor of the plane is also the stiffness **in plane strain** of a material of
space, of moduli ``(E, \nu)``: the strain ``\varepsilon_{33}`` is blocked, the
in-plane ``\lambda`` and ``\mu`` are those of the material, and
``k = \lambda + \mu``. The same tensor then has the two pairs, related by

```math
E_2 = \frac{E}{1 - \nu^2} , \qquad \nu_2 = \frac{\nu}{1 - \nu} ,
```

which are the moduli a plane-stress analysis would attribute to the plane-strain
problem. The keyword `convention` chooses the pair, `:plane_strain` by default
in the plane; ``k``, ``\lambda``, ``\mu`` and ``M = \lambda + 2\mu`` do not
depend on it.

```@example dims
E, ν = symbols("E ν", positive = true)
ℂ₂ = iso_stiffness_young_poisson(Val(2), E, ν)                 # plane strain of the material (E, ν)
simplify.(young_poisson(ℂ₂; convention = :intrinsic))          # (E₂, ν₂) of the plane medium
```

## The symmetry classes of the plane

Let the frame turn by an angle ``\theta``. The six constants of a stiffness of
the plane split into parts that turn at different rates: two invariants, and
two pairs that turn by ``2\theta`` and by ``4\theta``,

```math
\begin{aligned}
d_{12} &= C_{1112} + C_{2212} , & \delta &= \tfrac12\,(C_{1111} - C_{2222}) ,\\
a &= C_{1112} - C_{2212} , & b &= \tfrac14\,(C_{1111} + C_{2222} - 2C_{1122} - 4C_{1212}) ,
\end{aligned}
```

```math
\begin{aligned}
d_{12}' &= d_{12}\cos2\theta - \delta\sin2\theta , & \delta' &= \delta\cos2\theta + d_{12}\sin2\theta ,\\
a' &= a\cos4\theta - b\sin4\theta , & b' &= b\cos4\theta + a\sin4\theta .
\end{aligned}
```

They are checked on symbols, the components in the turned frame being those of
``\boldsymbol{Q}^{\mathsf T}\mathbb{C}\,\boldsymbol{Q}``:

```@example dims
C₁₁₁₁, C₂₂₂₂, C₁₁₂₂, C₁₂₁₂, C₁₁₁₂, C₂₂₁₂ = symbols("C1111 C2222 C1122 C1212 C1112 C2212", real = true)
r = sqrt(Sym(2))
M = [C₁₁₁₁ C₁₁₂₂ r*C₁₁₁₂; C₁₁₂₂ C₂₂₂₂ r*C₂₂₁₂; r*C₁₁₁₂ r*C₂₂₁₂ 2C₁₂₁₂]
Q = mandel_rotation(θ)
M′ = transpose(Q) * M * Q
a, b = C₁₁₁₂ - C₂₂₁₂, (C₁₁₁₁ + C₂₂₂₂ - 2C₁₁₂₂ - 4C₁₂₁₂) / 4
simplify((M′[1, 3] - M′[2, 3]) / r - (a * cos(4θ) - b * sin(4θ)))
```

A symmetry class is the set of tensors left unchanged by a group of rotations
and mirror lines; in the plane there are four, according to which of the two
pairs vanish:

| class | group | condition | constants | type |
| :--- | :--- | :--- | :---: | :--- |
| isotropic | all rotations | ``d_{12} = \delta = a = b = 0`` | 2 | `TensISO{4,2}` |
| square | of the square, 4 mirror lines | ``d_{12} = \delta = 0`` | 3 | `TensCubic{2}` |
| orthotropic | of the rectangle, 2 mirror lines | ``d_{12}' = a' = 0`` in a frame | 4 | `TensOrtho{4,2}` |
| general | none | — | 6 | `Tens` |

The orthotropic class has the Kelvin–Mandel matrix, in its frame,

```math
\begin{pmatrix} C_{11} & C_{12} & 0\\ C_{12} & C_{22} & 0\\ 0 & 0 & 2C_{66}\end{pmatrix} ,
\qquad C_{66} = C_{1212} ,
```

and the square class is the case ``C_{11} = C_{22}``. As the cubic class of
space, of which it is the restriction to a face of the cube
([Cubic symmetry](@ref th-cubic)), it is spanned by three orthogonal
projectors, here of trace 1 each,

```math
\mathbb{C} = \alpha\,\mathbb{J} + \beta\,\mathbb{E} + \gamma\,\mathbb{T} ,
\qquad \alpha = C_{11} + C_{12} ,\quad \beta = C_{11} - C_{12} ,\quad \gamma = 2C_{44} ,
```

so that sums, products and inverses act on ``(\alpha, \beta, \gamma)`` term by
term; the isotropic class is ``\beta = \gamma``.

```@example dims
c₁₁, c₁₂, c₄₄ = symbols("c11 c12 c44", positive = true)
ℚ = cubic_from_constants(c₁₁, c₁₂, c₄₄, Basis(θ))
simplify.(mandel_material(ℚ)), get_data(inv(ℚ))
```

!!! note "The plane has no transversely isotropic class"
    Transverse isotropy is invariance under the rotations about an axis. An
    axis normal to the plane gives all the rotations of the plane, the isotropic
    class; an axis in the plane gives the mirror line it carries and the one
    normal to it, the orthotropic class. `TensTI` is three-dimensional, and a
    transversely isotropic material seen in a plane that contains its axis is
    an orthotropic tensor of the plane.

## A tensor of space in one of its planes

Let ``\boldsymbol{e}_k`` be normal to a plane. In **plane strain** the strains
that carry the index ``k`` vanish, and the stiffness of the plane is the block
of the stiffness on the in-plane components, ``C_{ab}``. In **plane stress**
the stresses that carry ``k`` vanish; with the Kelvin–Mandel matrix split into
the in-plane components ``p`` and the others ``o``,

```math
\underline{\sigma}_o = 0
\;\Longrightarrow\;
\underline{\varepsilon}_o = -\,\mathbb{C}_{oo}^{-1}\,\mathbb{C}_{op}\,\underline{\varepsilon}_p ,
\qquad
\mathbb{C}' = \mathbb{C}_{pp} - \mathbb{C}_{po}\,\mathbb{C}_{oo}^{-1}\,\mathbb{C}_{op} .
```

When the plane is a plane of symmetry the out-of-plane shears are uncoupled
and only the normal component is eliminated,
``C'_{ab} = C_{ab} - C_{ak}C_{kb}/C_{kk}``. An isotropic material gives an
isotropic tensor of the plane in both cases:

```@example dims
ℂ₃ = iso_stiffness_young_poisson(E, ν)
simplify.(mandel(plane_strain_stiffness(ℂ₃)))
```

```@example dims
simplify.(mandel(plane_stress_stiffness(ℂ₃)))
```

An orthotropic material gives an orthotropic tensor of the plane in each of its
three planes of symmetry, a cubic one a square tensor in a face of its cube, and
a transversely isotropic one an isotropic tensor in the plane normal to its
axis and an orthotropic one in a plane that contains it.

## [Projections in the plane](@id th-projection-plane)

The closest tensor of a class to a given ``\mathbb{C}`` is its orthogonal
projection onto the class in a frame, at the angle ``\theta`` that makes the
distance least ([Projection onto a symmetry class](@ref th-projection)). In the
frame at ``\theta`` the projection drops the parts the class does not have, so
that the squared distances are, from the table above,

```math
d_\text{square}^2(\theta) = 2\,(d_{12}^2 + \delta^2) + 2\,a'(\theta)^2 ,
\qquad
d_\text{ortho}^2(\theta) = 2\,\bigl(d_{12}'(\theta)^2 + a'(\theta)^2\bigr) ,
```

``d_{12}^2 + \delta^2`` being invariant. Both are compared below with the
distance of the projection in a frame at an arbitrary angle:

```@example dims
using Random
𝔸 = rand(Xoshiro(1), Tens{4, 2}); m = mandel(𝔸); φ = 0.3
d₁₂ⁿ, δⁿ = (m[1, 3] + m[2, 3]) / √2, (m[1, 1] - m[2, 2]) / 2
aⁿ, bⁿ = (m[1, 3] - m[2, 3]) / √2, (m[1, 1] + m[2, 2] - 2m[1, 2] - 2m[3, 3]) / 4
a′, d₁₂′ = aⁿ * cos(4φ) - bⁿ * sin(4φ), d₁₂ⁿ * cos(2φ) - δⁿ * sin(2φ)
(project(𝔸, :CUBIC, Basis(φ))[2]^2, 2 * (d₁₂ⁿ^2 + δⁿ^2) + 2 * a′^2)
```

```@example dims
(project(𝔸, :ORTHO, Basis(φ))[2]^2, 2 * (d₁₂′^2 + a′^2))
```

In space the frame is found by a search over three angles; in the plane it is
known in closed form.

**Square class.** The distance is least where ``a' = 0``:

```math
\theta = \tfrac14\,\operatorname{atan2}(a, b) ,
```

up to an eighth of a turn, the diagonals of a square being mirror lines too.

**Orthotropic class.** With ``z = e^{4i\theta}``, ``c_1 = (d_{12} + i\delta)^2``
and ``c_2 = (a + ib)^2``,

```math
d_\text{ortho}^2 = |c_1| + |c_2| + \operatorname{Re}\bigl(c_1 z + c_2 z^2\bigr) ,
```

whose stationary points on the unit circle are roots of

```math
2c_2\,z^4 + c_1\,z^3 - \bar c_1\,z - 2\bar c_2 = 0 .
```

The least of them is taken: an exact global minimum, where a search would need
a starting point.

```@example dims
𝕆 = TensOrtho(10.0, 6.0, 3.0, 2.0, Basis(0.7))          # orthotropic, in a frame at 0.7 rad
t = Tens(get_array(𝕆))                                  # the same tensor, without its class
𝔹, dist, = project(t, :ORTHO)
dist, get_data(𝔹)
```

## Average over the rotations of the plane

The average of a tensor over all the rotations is what remains of it for a
uniform distribution of orientations. In space it is isotropic, even for a
tensor without the major symmetry. In the plane the rotations alone leave one
more tensor invariant: the deviators form a plane, which a rotation of the
frame turns by ``2\theta``, and the quarter turn ``\mathbb{W}`` of that plane
commutes with it. Hence

```math
\langle\mathbb{C}\rangle = \alpha\,\mathbb{J} + \beta\,\mathbb{K} + \omega\,\mathbb{W} ,
\qquad
\omega = \tfrac12\,\bigl(\underline e_1^{\mathsf T}\mathbb{C}\,\underline e_2 - \underline e_2^{\mathsf T}\mathbb{C}\,\underline e_1\bigr) ,
```

with ``\underline e_1 = (1,-1,0)/\sqrt2`` and ``\underline e_2 = (0,0,1)`` the
two deviators. ``\mathbb{W}`` has no major symmetry, so that ``\omega = 0`` for a
stiffness; it does not vanish for a concentration tensor. A mirror line
changes its sign: the average over the rotations *and* the reflections is the
isotropic one. `transverse_isotropify` is the first, `isotropify` the second;
at order 2 the same holds with the quarter turn of the plane itself.

## What is three-dimensional by nature

The vector product, the three Euler angles of a frame, the transversely
isotropic class and its Walpole basis, the cubic group of the cube, the
spherical, spheroidal, ellipsoidal and toroidal coordinates and their harmonics
exist in dimension 3 only. In the plane their counterparts are the single angle
of a frame, the square class and the polar coordinates; beyond dimension 3 a
frame is given by the matrix of its vectors, and the Cartesian system is the
coordinate system.
