# Matrix Functions API

```@meta
CurrentModule = MatrixFuns
```

## Matrix-variable functions `f(A::AbstractMatrix)::AbstractMatrix`
```@docs
MatrixFuns.mat_fun
```

## Divided Differences of `f(x::Real)::Number`
```@docs
MatrixFuns.div_diff
```

## Fréchet derivatives of `f(A::AbstractMatrix)::AbstractMatrix`
```@docs
MatrixFuns.mat_fun_frechet
```

## Implementation details
```@docs
MatrixFuns.get_splittings
MatrixFuns.split_cluster
MatrixFuns.split_by_sep
MatrixFuns.get_spread
MatrixFuns.checkspread
MatrixFuns.reorder_schur
MatrixFuns.get_swappings
MatrixFuns.atomic_block_fun!
MatrixFuns.atomic_block_fun
MatrixFuns.taylor_coeffs
MatrixFuns.parlett_recurrence
MatrixFuns.block_parlett_recurrence
MatrixFuns.DD_tensor
```