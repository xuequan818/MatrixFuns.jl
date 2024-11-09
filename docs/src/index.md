# MatrixFuns

A Julia package for computing scalar functions of matrix variables and their Fréchet derivatives. The matrix functions computation (for arbitrary square matrices) is based on the [Schur-Parlett algorithm]( https://doi.org/10.1137/S0895479802410815) (with improvements). The higher order Fréchet derivatives (for Hermitian matrices) are formulated similarly to the [Daleskii-Krein theorem](https://www.ams.org/books/trans2/047/), where the [divided differences](https://en.wikipedia.org/wiki/Divided_differences) are calculated accurately by the [Opitz' formula](https://www.emis.de/journals/SAT/papers/2/). In particular, `MatrixFuns` supports the computation of discontinuous functions. 

---
## Installation
`MatrixFuns` is currently an unregistered package and therefore needs to be downloaded or cloned to the user's local computer first, and then installed by running

```julia
julia> cd("your-local-path/MatrixFuns.jl")

julia> using Pkg

julia> Pkg.activate(".")

julia> Pkg.instantiate()
```