[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://xuequan818.github.io/MatrixFuns.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://xuequan818.github.io/MatrixFuns.jl/dev/)
[![Build Status](https://github.com/xuequan818/MatrixFuns.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/xuequan818/MatrixFuns.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/xuequan818/MatrixFuns.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/xuequan818/MatrixFuns.jl)

# MatrixFuns.jl
A Julia package for computing scalar functions of matrix variables and their Fréchet derivatives. The matrix functions computation (for arbitrary square matrices) is based on the [Schur-Parlett algorithm]( https://doi.org/10.1137/S0895479802410815) (with improvements). The higher order Fréchet derivatives (for Hermitian matrices) are formulated similarly to the [Daleskii-Krein theorem](https://www.ams.org/books/trans2/047/), where the [divided differences](https://en.wikipedia.org/wiki/Divided_differences) are calculated accurately by the [Opitz' formula](https://www.emis.de/journals/SAT/papers/2/). In particular, `MatrixFuns` supports the computation of discontinuous functions. 

## Examples
For smooth functions like `exp`:
```julia
julia> using LinearAlgebra

julia> using MatrixFuns

julia> A = [1 1 0; 0 1.1 1; 0 0 1.2]
3×3 Matrix{Float64}:
 1.0  1.0  0.0
 0.0  1.1  1.0
 0.0  0.0  1.2

julia> F0 = exp(A) # computed by Julia exponential matrix function
3×3 Matrix{Float64}:
 2.71828  2.85884  1.50334
 0.0      3.00417  3.15951
 0.0      0.0      3.32012

julia> F1 = mat_fun(exp, A) # computed by Schur-Parlett
3×3 Matrix{Float64}:
 2.71828  2.85884  1.50334
 0.0      3.00417  3.15951
 0.0      0.0      3.32012

julia> norm(F0-F1) # error of Schur-Parlett
2.920878410678069e-14

julia> X = 0.05 * Symmetric(A)# generate a symmetric matrix
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

For singular functions, such as Fermi-Dirac functions with temperatures close to 0, users can set a smaller `scale` to reduce the spread of each block to avoid large Taylor expansion errors near the singularities.
```julia
julia> f(x;μ) = 1/(1+exp(1000*(x-μ))); # Fermi-Dirac function with temperature equal 1e-3.

julia> f1(x) = f(x;μ=1.15);

julia> f2(x) = f(x;μ=0.05);

julia> ## returns the wrong result for default `scale`

julia> mat_fun(f1, A) 
3×3 Matrix{Float64}:
 1.0  6.48449e6  -9.72675e7
 0.0  6.4845e5   -1.2969e7
 0.0  0.0        -6.48449e5

julia> div_diff(f2, a)
-5.0306612072297075e147

julia> mat_fun_frechet(f2, X, hs)
3×3 Matrix{Float64}:
 1.84628e115  8.22808e114  -2.08877e115
 2.55222e115  9.24584e114  -3.37442e115
 1.80227e115  4.61108e114  -2.77086e115

julia> ## returns the correct result for smaller `scale`

julia> scale = 0.01;

julia> mat_fun(f1, A; scale)
3×3 Matrix{Float64}:
 1.0  0.0  -50.0
 0.0  1.0  -10.0
 0.0  0.0    1.92875e-22

julia> div_diff(f2, a; scale)
98.17057693049055

julia> mat_fun_frechet(f2, X, hs; scale)
3×3 Matrix{Float64}:
  19.7425    1.97425   -19.7425
   1.97425   0.197425   -1.97425
 -19.7425   -1.97425    19.7425
```

For piecewise functions consisting of continuous intervals, users can customize `color` function to first split the eigenvalues into different continuous intervals.
```julia
julia> x0 = 0.051 # discontinuous point

julia> g(x) = x < x0 ? sin(x) : (x+1)^2;

julia> color(x) = x < x0 ? 1 : 2; # specify a different color for each continuous interval

julia> B = triu(X)
3×3 Matrix{Float64}:
 0.05  0.05   0.0
 0.0   0.055  0.05
 0.0   0.0    0.06

julia> mat_fun(g, B; sep=0.001) # reference
3×3 Matrix{Float64}:
 0.0499792  10.6305   -52.6235
 0.0         1.11302    0.10575
 0.0         0.0        1.1236

julia> mat_fun(g, B) # returns the wrong result
3×3 Matrix{Float64}:
 1.1025  0.10525  0.0025
 0.0     1.11302  0.10575
 0.0     0.0      1.1236

julia> mat_fun(g, B; color) # returns the correct result
3×3 Matrix{Float64}:
 0.0499792  10.6305   -52.6235
 0.0         1.11303    0.10575
 0.0         0.0        1.1236
```

Specially, for constant functions, or piecewise functions consisting of constant intervals, users can further set `sep=Inf` to only split the eigenvalues by `color`, which will be more efficient.
```julia
julia> heaviside(x) = x < x0 ? 1 : 0 # Heaviside step function

julia> ## set `checknative=false` to avoid the cost of checking function nativity

julia> @time mat_fun(heaviside, B; color, checknative=false) 
  0.000099 seconds (100 allocations: 6.422 KiB)
3×3 Matrix{Float64}:
 1.0  -10.0  50.0
 0.0    0.0   0.0
 0.0    0.0   0.0

julia> @time mat_fun(heaviside, B; color, sep=Inf, checknative=false)
  0.000087 seconds (87 allocations: 5.641 KiB)
3×3 Matrix{Float64}:
 1.0  -10.0  50.0
 0.0    0.0   0.0
 0.0    0.0   0.0
```
For more details, please see the [documentation]().

## Installation
`MatrixFuns` is currently an unregistered package and therefore needs to be downloaded or cloned to the user's local computer first, and then installed by running

```julia
julia> cd("your-local-path/MatrixFuns.jl")

julia> using Pkg

julia> Pkg.activate(".")

julia> Pkg.instantiate()
```