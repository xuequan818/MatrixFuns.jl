# Divided difference

Given a function $f:\Omega\subset\mathbb{R}\to\mathbb{C}$ and $N+1$ data points $x_0, \dots, x_N\in\Omega$, the higher-order divided differences of $f$ are defined as 
```math
f[x_0, \dots, x_N]
= \begin{cases} 
(f[x_0,\dots,x_{N-1}]-f[x_1,\dots,x_N])/(x_0-x_N), &{\rm if}\,\,x_0\neq x_N, \\ 
\frac{\partial}{\partial z}f[z,x_1, \dots,x_{N-1}]{\big|}_{z=x_0}, & {\rm if}\,\,x_0= x_N,
\end{cases}
```
where $f[x_i] := f(x_i)$.

The divided-difference table of $f$ is an upper triangular matrix
```math
T_f(x_0,\dots,x_N) = 
\begin{bmatrix}
f[x_0] & f[x_0,x_1]& \cdots& f[x_0,\dots,x_N]\\
&f[x_1]&\ddots&\vdots\\
&&\ddots &f[x_{N-1},x_N]\\
&&&f[x_N]
\end{bmatrix}.
```
By the [Opitz' formula](https://www.emis.de/journals/SAT/papers/2/), $T_f$ can be created by 
```math
T_f(x_0,\dots,x_N) = f(J),
```
where 
```math
J = 
\begin{bmatrix}
x_0&1&&\\
&x_1&\ddots&\\
& &\ddots&1\\
&&&x_N
\end{bmatrix}.
```
Therefore, the accurate computation of the higher-order divided differences $f[x_0,\dots,x_N]$ is equivalent to the accurate computation of the matrix function $f(J)$.