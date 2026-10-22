# [Symbolic helpers and array algebra](@id api-symbolic)

Element-wise symbolic operations — no-ops on numeric types — and the raw-array
products the whole library is built on. Usage:
[Symbolic and numeric](@ref man-symbolic-numeric); theory:
[Tensor algebra](@ref th-tensor-algebra).

```@docs
tsimplify
tfactor
tsubs
tdiff
ttrigsimp
texpand_trig
tlimit
TensND.ApproxType
TensND.is_hard_numeric
carries_partials
primal
primal_type
dual_depth
dual_channels
float_like
isparallel
contract
qcontract
dotdot
otimesu
otimesul
sotimes
sboxtimes
isidentity
```
