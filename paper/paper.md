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
 - name: Antoine Levitt
   orcid: 0000-0002-3999-0289
   affiliation: 1
affiliations:
 - name: Laboratoire de Mathématiques d'Orsay, Université Paris-Saclay
   index: 1
 - name: School of Mathematical Sciences, Beijing Normal University
   index: 2
date: 5 June 2025
bibliography: paper.bib
---

# Summary

The computation of matrix functions (i.e., $f(A)$ for $A$ a $n \times n$ matrix and $f : \mathbb{C} \to \mathbb{C}$) and their Fréchet derivatives plays a crucial role in many fields of science [@Higham2008], and in particular in electronic structure calculations within density functional theory and response calculations. For Hermitian $A$, computing $f(A)$ can be done efficiently and stably by diagonalization. In the non-normal case, however, diagonalization is unstable and alternative schemes have to be used. Even in the Hermitian case, the evaluation of Fréchet derivatives requires (high-order) divided differences, which by Opitz's formula [@deBoor2005] is equivalent to the exact computation of $f(A)$ for non-normal $A$.

In this work, we develop `MatrixFuns.jl` a Julia package [@julia] to provide the robust computation of matrix functions for arbitrary square matrices and higher-order Fréchet derivatives for Hermitian matrices. This package is tailored towards high accuracy with relatively small matrices and relatively complicated functions $f$. Our work is based on the Schur-Parlett algorithm [@DaviesHigham03;@HighamMohy10], with the following modifications:
* It supports functions that are discontinuous, or have sharp variations.
* It does not require the computation of arbitrary-order derivatives of $f$.
* It exploits existing special-purpose methods for computing matrix functions (e.g. for functions involving exponentials or logarithms) when they exist.


# Statement of need
`MatrixFuns.jl` aims to provide high-accuracy computations for general matrix functions and arbitrary-order Fréchet derivatives (including divided differences) in Julia. Julia provides some native matrix functions, but the choice is limited to a few functions for which special-purpose algorithms exist (e.g., exponentials, logarithms, matrix powers...). There are no dedicated functions in Julia for computing Fréchet derivatives and divided differences; some Julia packages offer tools for their computation (e.g., `ChainRules.jl` [@ChainRules.jl], `DFTK.jl` [@DFTKjcon]...), but are typically limited to first order.

# Methods
## Matrix functions
The basic principle of the Schur-Parlett algorithm is as follows. First, one performs a Schur decomposition to reduce to the case of an upper triangular matrix. Then, one uses the Parlett recursion, which for a block matrix $A = \begin{pmatrix}A_{11}&A_{12}\\0&A_{22}\end{pmatrix}$ expresses $B = f(A)$ as $B_{11}=f(A_{11})$, $B_{22}=f(A_{22})$ and $B_{12}$ given by the solution of the Sylvester equation $A_{11} B_{12} - B_{12} A_{22} = B_{11}A_{12}-A_{12}B_{22}$. In principle, this can be used to compute $f(A)$ by a recursion, but the Sylvester equation becomes ill-conditioned when $A_{11}$ and $A_{22}$ do not have well-separated eigenvalues. In this case, one can use Taylor series, as proposed in [@DaviesHigham03;@HighamMohy10], but this has the disadvantage of requiring arbitrarily many derivatives of $f$, which might be impractical in some applications (e.g. when the function is not analytic, or has sharp variations).

Our algorithm attempts to find a partition of the eigenvalues of $A$ (computed using a Schur decomposition) into blocks that are well-separated. The diagonal blocks are then computed using Taylor series, and the Parlett recursion is used to fill out the off-diagonal blocks. The partitioning aims to find small blocks (so that low-order Taylor series can be used) that are well-separated (so that the Parlett recursion is well-conditioned)

To find the partition, we start by partitioning the set of eigenvalues $\Lambda$ into disjoint clusters $\Lambda_i$ such that the distance between two such clusters is at least ${\rm sep}$, where ${\rm sep}$ is a user-definable parameter. We then check if the partition is acceptable by estimating the error in all the clusters; if the estimated error is acceptable, we accept the partition; if not, we split the unacceptable clusters further by applying the partitioning algorithm recursively to each unacceptable $\Lambda_i$. We estimate the error in a cluster $\Lambda_i$ of diameter $d_i$ as ${\textrm{err}}_i=(\frac{d_i}{\textrm{scale}})^{\textrm{max}\_\textrm{deg}+1}$. We accept a cluster if $\textrm{err}_i < \varepsilon/\textrm{sep}$. This choice is made to balance the error originating from the Taylor expansion within a cluster $\textrm{err}_i$ with the error incurred by the use of the Parlett recursion $\varepsilon/\textrm{sep}$.

Therefore, our algorithm has the following parameters: 
* ${\rm scale}$, the characteristic scale of variations of $f$, set to $1$ by default.
* ${\rm max}\_{\rm deg}$, the order of the Taylor series used, which should be set by the user according to the regularity of the function under consideration and the feasibility of computing high-order derivatives (computed automatically using `TaylorSeries.jl` [@TaylorSeries.jl] and `Arblib.jl` [@Arblib.jl], where the latter is faster in calculating much larger orders and supports some special functions from `SpecialFunctions.jl` [@SpecialFunctions.jl]). By default, set to a large value.
* ${\rm sep}$, the initial separation distance, set to ${\rm 0.1*scale}$ by default following [@DaviesHigham03;@HighamMohy10].
* $\varepsilon$, the target accuracy, set to machine accuracy by default.

In the case where Julia natively supports the computation of $f(A)$ (as determined by trying to compute `f(ones(1,1))` and catching any resulting error), we use them instead of Taylor series to compute diagonal blocks. In the error estimate, we consider ${\rm max\_deg}=\infty$, and therefore use a partition with maximal diameter ${\rm scale}$. We partition the eigenvalues rather than simply call the native $f(A)$, because $f$ can still have sharp variations, which would cause inaccuracies in $f(A)$. For example:
```julia
f(x) = I/(I+exp(50*x));

A = [-0.1 10.0 0.0; 0.0 1 5.0; 0.0 0.0 -0.11];

f(A) # native call
3×3 Matrix{Float64}:
 0.993307  -9.03006      -3.33299e7
 0.0        1.92875e-22  -4.48617
 0.0        0.0           0.99593

mat_fun(f, A; scale=1/50) # Schur-Parlett
3×3 Matrix{Float64}:
 0.993307  -9.03006      -28.8619
 0.0        1.92875e-22   -4.48617
 0.0        0.0            0.99593
```

For discontinuous functions, or functions with sharp variations, our algorithm takes as input a color mapping ${\textrm{color}}:\mathbb{C}\to\mathbb{Z}, \lambda\mapsto a$, and makes sure that all the eigenvalues inside a cluster have the same color. This ensures that Talor expansions are not used across the discontinuity boundaries.

## Fréchet derivatives
For a Hermitian $A\in\mathbb{C}^{n\times n}$, denote the eigenpairs by $\{(\lambda_i,v_i)\}$. The $N$-th order Fréchet derivative expresses the variation of $f(A)$ with respect to a set of variations $H_1, \dots, H_N$, and is given by (see the documentation of `MatrixFuns.jl` for details)
$$
{{\rm d}}^{N}f(A)H_1\cdots H_N =\sum_{i_0,\cdots,i_{N}=1}^nv_{i_0}\Bigg(\sum_{p\in\mathcal{P}_N}(H_{p(1)})_{i_0,i_1}\cdots (H_{p(N)})_{i_{N-1},i_{N}}\Bigg){f[\lambda_{i_0},\cdots,\lambda_{i_{N}}]}v_{i_{N}}^*,
$$
where $(H_{p(k)})_{i,j}=v_i^*H_{p(k)}v_j$ and $p\in\mathcal{P}_N$ is an arbitrary permutation of $\{1,\cdots,N\}$. The higher-order divided differences $f[x_0, \dots, x_N]$ defined recursively by
$$
f[x_0, \dots, x_N]
    	= \begin{cases} 
    		(f[x_0,\dots,x_{N-1}]-f[x_1,\dots,x_{N}])/(x_0-x_N), &{\rm if}\,\,x_0\neq x_N,\\ 
    		\frac{\partial}{\partial z}f[z,x_1, \dots,x_{N-1}]{\big|}_{z=x_0}, & {\rm if}\,\,x_0= x_N.
    	\end{cases}
$$
The naive evaluation of this recurrence formula is prone to numerical stabilities. Instead, we compute the divided differences using Opitz's formula
$$
f\left(\begin{bmatrix}
    		x_0&1&&\\
    		&x_1&\ddots&\\
    		& &\ddots&1\\
    		&&&x_N
    	\end{bmatrix}\right) = 
    	\begin{bmatrix}
    		f[x_0] & f[x_0,x_1]& \cdots& f[x_0,\dots,x_N]\\
    		&f[x_1]&\ddots&\vdots\\
    		&&\ddots &f[x_{N-1},x_N]\\
    		&&&f[x_N]
    	\end{bmatrix}.
$$
Therefore, the key point in evaluating the Fréchet derivative reduces to computing matrix functions for upper triangular matrices.

# Examples
We first show how to use MatrixFuns.jl to compute the matrix functions, divided differences, and Fréchet derivatives for smooth functions such as `exp`.
```julia
using MatrixFuns

A = [-0.1 1.0 0.0; 0.0 -0.05 1.0; 0.0 0.0 0.01];

mat_fun(exp, A) # returns exp(A)
3×3 Matrix{Float64}:
 0.904837  0.92784   0.477323
 0.0       0.951229  0.980346
 0.0       0.0       1.01005

div_diff(exp, -0.1, -0.05, 0.01) # returns exp[-0.1,-0.05,0.01]
0.47732345844677654

H = 0.5 * (A + A'); # generates a Hermitian matrix

hs = map(i -> i * H, [1, 2]);

mat_fun_frechet(exp, H, hs) # returns d^2exp(H)hs[1]hs[2]
3×3 Matrix{Float64}:
 0.519468  0.347941  0.55445
 0.347941  1.10871   0.46992
 0.55445   0.46992   0.610653
 ```


In addition to the usual smooth functions, MatrixFuns.jl can also support special functions and discontinuous functions. Here, we use the error function `erf` and the sign function `sign` to show how it can be used to handle functions with different smoothness.
```julia
using MatrixFuns, SpecialFunctions

A = [-0.1 1.0 0.0; 0.0 -0.05 1.0; 0.0 0.0 0.01];

mat_fun(erf, A) # smooth function
3×3 Matrix{Float64}:
 -0.112463   1.12182   0.0524648
  0.0       -0.056372  1.12759
  0.0        0.0       0.0112834

mat_fun(x -> erf(500x), A; scale=1/500, color=x->x<0 ? 1 : 2) # singular function
3×3 Matrix{Float64}:
 -1.0   0.0  303.03
  0.0  -1.0   33.3333
  0.0   0.0    1.0

mat_fun(sign, A; color=x->Int(sign(x))) # discontinuous function with smooth branches
3×3 Matrix{Float64}:
 -1.0   0.0  303.03
  0.0  -1.0   33.3333
  0.0   0.0    1.0
```

# Reference
