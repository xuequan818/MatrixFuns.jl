[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://xuequan818.github.io/MatrixFuns.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://xuequan818.github.io/MatrixFuns.jl/dev/)
[![Build Status](https://github.com/xuequan818/MatrixFuns.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/xuequan818/MatrixFuns.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/xuequan818/MatrixFuns.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/xuequan818/MatrixFuns.jl)

# MatrixFuns.jl
A Julia package for computing scalar functions of matrix variables and their Fréchet derivatives. The matrix functions computation (for arbitrary square matrices) is based on the [Schur-Parlett algorithm](https://doi.org/10.1137/S0895479802410815) (with improvements)a. The higher order Fréchet derivatives (for Hermitian matrices) are formulated similarly to the [Daleckii-Krein theorem](https://www.ams.org/books/trans2/047/), where the [divided differences](https://en.wikipedia.org/wiki/Divided_differences) are calculated accurately by the [Opitz' formula](https://www.emis.de/journals/SAT/papers/2/). In particular, `MatrixFuns` supports the computation of discontinuous functions. 

## Examples
For smooth functions, such as `exp`:
```julia
julia> using MatrixFuns

julia> using LinearAlgebra

julia> A = [1 1 0; 0 1.1 1; 0 0 1.2]
3×3 Matrix{Float64}:
 1.0  1.0  0.0
 0.0  1.1  1.0
 0.0  0.0  1.2

julia> mat_fun(exp, A) # computed by Schur-Parlett
3×3 Matrix{Float64}:
 2.71828  2.85884  1.50334
 0.0      3.00417  3.15951
 0.0      0.0      3.32012

julia> X = 0.05 * Symmetric(A) # generate a symmetric matrix
3×3 Symmetric{Float64, Matrix{Float64}}:
 0.05  0.05   0.0
 0.05  0.055  0.05
 0.0   0.05   0.06

julia> a = eigvals(X)
3-element Vector{Float64}:
 -0.01588723439378913
  0.055000000000000014
  0.12588723439378913

julia> div_diff(exp, a) # returns the 2nd order divided difference
0.5284915575854794

julia> hs = map(x->x*X, [1, 2])
2-element Vector{Symmetric{Float64, Matrix{Float64}}}:
 [0.05 0.05 0.0; 0.05 0.05500000000000001 0.05; 0.0 0.05 0.06]
 [0.1 0.1 0.0; 0.1 0.11000000000000001 0.1; 0.0 0.1 0.12]

julia> mat_fun_frechet(exp, X, hs) # returns the 2nd order Fréchet derivative d^2exp(X)hs_1hs_2
3×3 Matrix{Float64}:
 0.0110862   0.0119138  0.00588556
 0.0119138   0.0181632  0.0130909
 0.00588556  0.0130909  0.0135867
```
For special functions, such as the error function `erf`: 
```julia
julia> using SpecialFunctions

julia> mat_fun(erf, A)
3×3 Matrix{Float64}:
 0.842701  0.375043  -0.369768
 0.0       0.880205   0.301089
 0.0       0.0        0.910314
```
For singular functions, such as [Fermi-Dirac](https://en.wikipedia.org/wiki/Fermi%E2%80%93Dirac_statistics) functions with temperatures close to 0, user can set a smaller `scale` to reduce the spread of each block to avoid large Taylor expansion errors near the singularities, and can also customize `color` function to avoid the block's spectral range including singularities.
```julia
julia> μ = 1.15 

julia> f(x) = 1/(1+exp(1000*(x-μ))); # Fermi-Dirac function with temperature equal to 1e-3

julia> color(x) = x < μ ? 1 : 2; # use `color` to avoid singularities.

julia> mat_fun(f, A; scale=1/1000, color)
3×3 Matrix{Float64}:
 1.0  0.0  -50.0
 0.0  1.0  -10.0
 0.0  0.0    1.92875e-22
```

For piecewise functions consisting of continuous intervals, user can customize `color` function to first split the eigenvalues into different continuous intervals.
```julia
julia> x0 = 0.051 # discontinuous point

julia> g(x) = x < x0 ? sin(x) : (x+1)^2;

julia> color(x) = x < x0 ? 1 : 2; # specify a different color for each continuous interval

julia> mat_fun_frechet(g, X, hs; color)
3×3 Matrix{Float64}:
 0.019713    0.0213782  0.00975085
 0.0213782   0.0316017  0.0233283
 0.00975085  0.0233283  0.0241837
```

For more details, please see the [documentation](https://xuequan818.github.io/MatrixFuns.jl/dev/).

## Installation

```julia
julia> using Pkg

julia> Pkg.add("MatrixFuns")
```