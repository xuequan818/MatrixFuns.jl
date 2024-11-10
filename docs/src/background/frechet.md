# Fréchet derivative 

## Contour integral form
Let $\mathcal{H}:=\mathbb{C}^{N\times N}_{\rm herm}$ be the vector space of $N\times N$ Hermitian matrices. For $H\in\mathcal{H}$ and $h_1,...h_n\in\mathcal{H}$, and $f$ is an $n$ times continuously differentiable function on a subset of $\mathbb{C}$ containing the spectrum of $H+t_1h_1+\cdots + t_nh_n$, the $n$-th order Fréchet derivative of $f(H)$ is
```math
\begin{align*}
    &{{\rm d}}^{n}f(H)h_1\cdots h_n =\\ &\frac{1}{2\pi i}\oint_\mathcal{C} f(z) \sum_{p\in\mathcal{P}_n}(z-H)^{-1}h_{p(1)}(z-H)^{-1}\cdots(z-H)^{-1}h_{p(n)}(z-H)^{-1}\; dz,
\end{align*}
```
where $\mathcal{C}$ is a contour in the complex plane enclosing all the eigenvalues of $H$, and $p\in\mathcal{P}_n$ is an arbitrary permutation of $\{1,\cdots,n\}$. This can be proved by induction.

For $n = 1$, we have 
```math
\begin{align*}
    {\rm d} f(H)h&=\lim_{t\to 0} \frac{f(H+th)-f(H)}{t}\\&=\frac{1}{2\pi i}\oint_\mathcal{C} f(z)\lim_{t\to 0}\frac{(z-H-th)^{-1}-(z-H)^{-1}}{t}\;dz.
\end{align*}
```
Note that
```math
\begin{align*}
    (z-H-th)^{-1}  &= (I-t(z-H)^{-1}h)^{-1}(z-H)^{-1} \\&= (I+t(z-H)^{-1}h+O(t^2))(z-H)^{-1},
\end{align*}
```
we have 
```math
\begin{align*}
    {\rm d} f(H)h=\frac{1}{2\pi i}\oint_\mathcal{C} f(z)(z-H)^{-1}h(z-H)^{-1}\;dz,
\end{align*}
```
which satisfies the formula.

Assume the $n-1$-th order derivative satisfies the formula. Then we have the $n$-th order derivative
```math
\begin{align*}
    {{\rm d}}^{n}f(H)h_1\cdots h_n &=\lim_{t\to 0} \frac{{{\rm d}}^{n-1}f(H+th_n)h_1\cdots h_{n-1}-{{\rm d}}^{n-1}f(H)h_1\cdots h_{n-1}}{t}\\
    &=\frac{1}{2\pi i}\oint_\mathcal{C} f(z)\lim_{t\to 0}\sum_{p\in\mathcal{P}_{n-1}}\frac{1}{t}\Big((z-H-th_n)^{-1}h_{p(1)}\cdots h_{p(n-1)}(z-H-th_n)^{-1}\\
    &\qquad -(z-H)^{-1}h_{p(1)}\cdots h_{p(n-1)}(z-H)^{-1}\Big)\;dz.
\end{align*}
```
Similarly, we have
```math 
\begin{align*}
   (z&-H-th_n)^{-1}h_{p(1)}\cdots h_{p(n-1)}(z-H-th_n)^{-1}\\
    &=(I+t(z-H)^{-1}h_n)(z-H)^{-1}h_{p(1)}\cdots h_{p(n-1)}(z-H)^{-1}(I+th_n(z-H)^{-1}) + O(t^2)\\
    &=(z-H)^{-1}h_{p(1)}\cdots h_{p(n-1)}(z-H)^{-1} \\
    &\quad+ t\Big((z-H)^{-1}h_n(z-H)^{-1}h_{p(1)}\cdots h_{p(n-1)}(z-H)^{-1}\\
    &\qquad\quad+(z-H)^{-1}h_{p(1)}(z-H)^{-1}h_n\cdots h_{p(n-1)}(z-H)^{-1}\\
    &\qquad\quad+(z-H)^{-1}h_{p(1)}(z-H)^{-1}h_{p(2)}\cdots h_n(z-H)^{-1}\Big) + O(t^2).
\end{align*}
```
Therefore, we can obtain the formula.

## Divided difference form
Let $H=\Phi \Lambda\Phi^{-1}=\sum_i^N\lambda_i\phi_i\phi_i^{-1}$, where $\phi_i$ is the $i$-th column of $\Phi$ and $\phi_i^{-1}$ is the $i$-th row of $\Phi^{-1}$, and assume there is no degeneration. Then for $p\in\mathcal{P}_n$ we have
```math 
\begin{align*}
    (z-&H)^{-1}h_{p(1)}(z-H)^{-1}\cdots(z-H)^{-1}h_{p(n)}(z-H)^{-1}\\
    &=\sum_{i_0,\cdots,i_{n}=1}^N\phi_{i_0}(h_{p(1)})_{i_0,i_1}\cdots (h_{p(n)})_{i_{n-1},i_n}\phi_{i_{n}}^{-1}(z-\lambda_{i_0})^{-1}\cdots (z-\lambda_{i_{n}})^{-1},
\end{align*}
```
where $(h_{p(k)})_{i,j}=\phi_i^{-1}h_{p(k)}\phi_j$. We let
```math 
    (z-\lambda_{i_0})^{-1}\cdots (z-\lambda_{i_{n}})^{-1} = \sum_{k=0}^{n}C_k (z-\lambda_{i_k})^{-1},
```
then
```math 
\sum_{k=0}^{n}C_k \prod_{\ell\neq k}(z-\lambda_{i_\ell})=1.
```
Let $z=\lambda_{i_k}$, we can obtain
```math 
C_k = \frac{1}{\prod_{\ell\neq k}(\lambda_{i_k}-\lambda_{i_\ell})}.
```
Therefore, we have
```math 
\begin{align*}
    \frac{1}{2\pi i}\oint_\mathcal{C} f(z)(z-\lambda_{i_0})^{-1}\cdots (z-\lambda_{i_{n}})^{-1}\;dz = \sum_{k=0}^{n}\frac{f(\lambda_{i_k})}{\prod_{\ell\neq k}(\lambda_{i_k}-\lambda_{i_\ell})}=f[\lambda_{i_0},\cdots,\lambda_{i_{n}}].
\end{align*}
```
Finally, we obtain
```math 
\begin{equation}
    {{\rm d}}^{n}f(H)h_1\cdots h_n =\sum_{i_0,\cdots,i_{n}=1}^N\phi_{i_0}\Bigg(\sum_{p\in\mathcal{P}_n}(h_{p(1)})_{i_0,i_1}\cdots (h_{p(n)})_{i_{n-1},i_{n}}\Bigg)f[\lambda_{i_0},\cdots,\lambda_{i_{n}}]\phi_{i_{n}}^{-1},
\end{equation}
```
which can be also written as 
```math 
\begin{equation}
(\Phi^{-1}[{{\rm d}}^{n}f(H)h_1\cdots h_n]\Phi)_{k\ell}=\sum_{i_1,\cdots,i_{n-1}=1}^N\Bigg(\sum_{p\in\mathcal{P}_n}(h_{p(1)})_{k,i_1}\cdots (h_{p(n)})_{i_{n-1},\ell}\Bigg)f[\lambda_k,\lambda_{i_1},\cdots,\lambda_{i_{n-1}},\lambda_\ell].
\end{equation}
```

## Array operations
Use array operations to efficiently compute the Fréchet derivative. For simplicity, just consider the no permutation case and define
```math
(F_n)_{kℓ}:=∑_{i_1,⋯,i_{n-1}=1}^N(h_1)_{k,i_1}⋯ (h_n)_{i_{n-1},ℓ}Λ^{0,1,…,n-1,n}_{k,i_1,…,i_{n-1},ℓ},
```
where $Λ^{0,…,n}_{i_0,…,i_n} := f[λ_{i_0},⋯,λ_{i_n}]$. It is immediately to obtain that 
```math
F_1 =  h_1 ∘ Λ^{0,1},
```
and 
```math
F_2 = ∑_{i=1}^N (\mathfrak{h}^{1,2} ∘ Λ^{0,1,2})_{:,i,:},
```
 with $\mathfrak{h}^{1,2}_{:,i,:} := (h_1)_{:,i}(h_2)_{i,:}$. For $n ≥ 3$, first compute 
```math
\mathfrak{F}^{0,2,…,n}_{:,:,j_3,…,j_n} := ∑_{i=1}^N (\mathfrak{h}^{1,2} ∘ Λ^{0,1,…,n}_{:,:,:,j_3,…,j_n})_{:,i,:}
```
and permute such that the $0$-dimension is at the end $\mathfrak{F}^{2,…,n,0}$. Then for $m ≥ 3$, there is the recursion
```math
\mathfrak{F}^{m,…,n,0}_{:,j_m,…,j_n} = ∑_{i=1}^N (h_m ∘ \mathfrak{F}^{m-1,…,n,0}_{:,:,j_m,…,j_n})_{i,:}
```
and $F_n = (\mathfrak{F}^{n,0})^T$.
