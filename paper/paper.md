---
title: 'MatrixFuns.jl: Matrix functions in Julia'
tags:
  - Matrix functions
  - Fréchet derivatives
  - Divided differences
  - Julia
authors:
 - name: Xue Quan
   orcid: 0000-0003-3349-9517
   affiliation: "1, 2"
 - name: Atonie Levitt
   orcid: 0000-0002-3999-0289
   affiliation: 1
affiliations:
 - name: Laboratoire de Mathématiques D'Orsay, Université Paris-Saclay
   index: 1
 - name: School of Mathematical Sciences, Beijing Normal University
   index: 2
date: 18 November 2024
bibliography: paper.bib
---

# Summary

The computation of matrix-variable functions (i.e., $f:\mathbb{C}^{n\times n}\to\mathbb{C}^{n\times n}$, $A\mapsto f(A)$) and their Fréchet derivatives plays a crucial role in electronic structure calculations, especially within density functional perturbation theory and response theory. 
Although it is trivial to compute $f(A)$ for Hermitian $A$, the Fréchete derivatives requires an accurate evaluation of the higher-order divided differences, which by Opitz's formula[@deBoor2005] is equivalent to the exact computation of $f(A)$ for non-normal $A$. 

In this work, we develop [MatrixFuns.jl](https://github.com/xuequan818/MatrixFuns.jl) a Julia package[@julia] to provide the robust computation of matrix functions for arbitrary square matrices and higher-order Fréchet derivatives for Hermitian matrices. 
The computation of matrix functions is based on the Schur-Parlett algorithm[@DaviesHigham03,HighamMohy10], but with some improvements so that it can also support discontinuous functions.  

# Statement of need
MatrixFuns.jl aims to provide convincing computations for general matrix functions and arbitrary-order Fréchet derivatives (including divided differences) in Julia.
Although Julia provides some native matrix functions, the choice is limited. 
Moreover, while some Julia packages offer tools for calculating Fréchet derivatives and divided differences, these are typically restricted to first-order computations. 
In MATLAB, the [funm](https://www.mathworks.com/help/symbolic/sym.funm.html) function supports the evaluation of general matrix functions but does not handle discontinuous functions. Additionally, MATLAB also lacks functions to compute higher-order Fréchet derivatives.


# Example
We first show how to use MatrixFuns.jl to compute the matrix functions, divided differences, and Fréchet derivatives for smooth functions such as `exp`.
```julia
using MatrixFuns

A = [-0.1 1.0 0.0; 0.0 -0.05 1.0; 0.0 0.0 0.01];

mat_fun(exp, A) # returns the matrix function exp(A)
3×3 Matrix{Float64}:
 0.904837  0.92784   0.477323
 0.0       0.951229  0.980346
 0.0       0.0       1.01005

div_diff(exp, -0.1, -0.05, 0.01) # returns the second-order divided difference exp[-0.1,-0.05,0.01]
0.47732345844677654

H = 0.5 * (A + A'); # generates a Hermitian matrix

hs = map(i -> i * H, [1, 2]);

mat_fun_frechet(exp, H, hs) # returns the second-order Fréchet derivative d^2exp(H)hs[1]hs[2]
3×3 Matrix{Float64}:
 0.519468  0.347941  0.55445
 0.347941  1.10871   0.46992
 0.55445   0.46992   0.610653
 ```

 
In addition to the usual smooth functions, MatrixFuns.jl can also support special functions and discontinuous functions. Here, we use the error function `erf` and the sign function `sign` to show how it can be usded to handle functions with different smoothness.
```julia
using MatrixFuns, SpecialFunctions

A = [-0.1 1.0 0.0; 0.0 -0.05 1.0; 0.0 0.0 0.01];

mat_fun(erf, A) # smooth function
3×3 Matrix{Float64}:
 -0.112463   1.12182   0.0524648
  0.0       -0.056372  1.12759
  0.0        0.0       0.0112834

mat_fun(x -> erf(500x), A; scale=0.1) # singular function
3×3 Matrix{Float64}:
 -1.0   0.0  303.03
  0.0  -1.0   33.3333
  0.0   0.0    1.0

mat_fun(sign, A; sep=Inf, color=x->Int(sign(x))) # discontinuous function
3×3 Matrix{Float64}:
 -1.0   0.0  303.03
  0.0  -1.0   33.3333
  0.0   0.0    1.0
```

# Reference
