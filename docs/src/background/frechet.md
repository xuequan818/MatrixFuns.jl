# Fréchet derivative 

## Contour integral form
Let $A\in \mathbb{C}^{n\times n}$ be a Hermitian matrix, $H_1, \dots, H_N$ be a set of Hermitian variations, and $f$ be an $N$ times continuously differentiable function on a subset of $\mathbb{C}$ containing the spectrum of $A+t_1H_1+\cdots + t_NH_N$, the $N$-th order Fréchet derivative of $f(A)$ is
```math
\begin{align*}
    &{{\rm d}}^{N}f(A)H_1\cdots H_N =\\ &\frac{1}{2\pi i}\oint_\mathcal{C} f(z) \sum_{p\in\mathcal{P}_N}(z-A)^{-1}H_{p(1)}(z-A)^{-1}\cdots(z-A)^{-1}H_{p(N)}(z-A)^{-1}\; dz,
\end{align*}
```
where $\mathcal{C}$ is a contour in the complex plane enclosing all the eigenvalues of $A$, and $p\in\mathcal{P}_N$ is an arbitrary permutation of $\{1,\cdots,N\}$. This can be proved by induction.

For $N = 1$, we have 
```math
\begin{align*}
    {\rm d} f(A)H&=\lim_{t\to 0} \frac{f(A+tH)-f(A)}{t}\\&=\frac{1}{2\pi i}\oint_\mathcal{C} f(z)\lim_{t\to 0}\frac{(z-A-tH)^{-1}-(z-A)^{-1}}{t}\;dz.
\end{align*}
```
Note that
```math
\begin{align*}
    (z-A-tH)^{-1}  &= (I-t(z-A)^{-1}H)^{-1}(z-A)^{-1} \\&= (I+t(z-A)^{-1}H+O(t^2))(z-A)^{-1},
\end{align*}
```
we have 
```math
\begin{align*}
    {\rm d} f(A)H=\frac{1}{2\pi i}\oint_\mathcal{C} f(z)(z-A)^{-1}H(z-A)^{-1}\;dz,
\end{align*}
```
which satisfies the formula.

Assume the $(N-1)$-th order derivative satisfies the formula. Then we have the $N$-th order derivative
```math
\begin{align*}
    {{\rm d}}^{N}f(A)H_1\cdots H_N &=\lim_{t\to 0} \frac{{{\rm d}}^{N-1}f(A+tH_N)H_1\cdots H_{N-1}-{{\rm d}}^{N-1}f(A)H_1\cdots H_{N-1}}{t}\\
    &=\frac{1}{2\pi i}\oint_\mathcal{C} f(z)\lim_{t\to 0}\sum_{p\in\mathcal{P}_{N-1}}\frac{1}{t}\Big((z-A-tH_N)^{-1}H_{p(1)}\cdots H_{p(N-1)}(z-A-tH_N)^{-1}\\
    &\qquad -(z-A)^{-1}H_{p(1)}\cdots H_{p(N-1)}(z-A)^{-1}\Big)\;dz.
\end{align*}
```
Similarly, we have
```math 
\begin{align*}
   (z&-H-tH_N)^{-1}H_{p(1)}\cdots H_{p(N-1)}(z-A-tH_N)^{-1}\\
    &=(I+t(z-A)^{-1}H_N)(z-A)^{-1}H_{p(1)}\cdots H_{p(N-1)}(z-A)^{-1}(I+tH_N(z-A)^{-1}) + O(t^2)\\
    &=(z-A)^{-1}H_{p(1)}\cdots H_{p(N-1)}(z-A)^{-1} \\
    &\quad+ t\Big((z-A)^{-1}H_N(z-A)^{-1}H_{p(1)}\cdots H_{p(N-1)}(z-A)^{-1}\\
    &\qquad\quad+(z-A)^{-1}H_{p(1)}(z-A)^{-1}H_N\cdots H_{p(N-1)}(z-A)^{-1}\\
    &\qquad\quad+\cdots+(z-A)^{-1}H_{p(1)}(z-A)^{-1}H_{p(2)}\cdots H_N(z-A)^{-1}\Big) + O(t^2).
\end{align*}
```
This completes the proof.

## Divided difference form
Let $\{(\lambda_i,v_i)\}$ be the eigenpairs of $A$, and assume there is no degeneration. Then for $p\in\mathcal{P}_N$, we have
```math 
\begin{align*}
    (z-&A)^{-1}H_{p(1)}(z-A)^{-1}\cdots(z-A)^{-1}H_{p(N)}(z-A)^{-1}\\
    &=\sum_{i_0,\cdots,i_N=1}^nv_{i_0}(H_{p(1)})_{i_0,i_1}\cdots (H_{p(N)})_{i_{N-1},i_N}v_{i_N}^*(z-\lambda_{i_0})^{-1}\cdots (z-\lambda_{i_N})^{-1},
\end{align*}
```
where $(H_{p(k)})_{i,j}=v_i^*H_{p(k)}v_j$. Let
```math 
    (z-\lambda_{i_0})^{-1}\cdots (z-\lambda_{i_N})^{-1} = \sum_{k=0}^{N}C_k (z-\lambda_{i_k})^{-1},
```
then we have
```math 
\sum_{k=0}^{N}C_k \prod_{\ell\neq k}(z-\lambda_{i_\ell})=1.
```
Let $z=\lambda_{i_k}$, we can obtain
```math 
C_k = \frac{1}{\prod_{\ell\neq k}(\lambda_{i_k}-\lambda_{i_\ell})}.
```
Therefore, we have
```math 
\begin{align*}
    \frac{1}{2\pi i}\oint_\mathcal{C} f(z)(z-\lambda_{i_0})^{-1}\cdots (z-\lambda_{i_N})^{-1}\;dz = \sum_{k=0}^{N}\frac{f(\lambda_{i_k})}{\prod_{\ell\neq k}(\lambda_{i_k}-\lambda_{i_\ell})}=f[\lambda_{i_0},\cdots,\lambda_{i_N}].
\end{align*}
```
Finally, we obtain
```math 
\begin{equation}
    {{\rm d}}^{N}f(A)H_1\cdots H_N =\sum_{i_0,\cdots,i_N=1}^nv_{i_0}\Bigg(\sum_{p\in\mathcal{P}_N}(H_{p(1)})_{i_0,i_1}\cdots (H_{p(N)})_{i_{N-1},i_N}\Bigg)f[\lambda_{i_0},\cdots,\lambda_{i_N}]v_{i_N}^*,
\end{equation}
```
which can be also written as 
```math 
\begin{equation}
(U^*[{{\rm d}}^{N}f(A)H_1\cdots H_N]U)_{k\ell}=\sum_{i_1,\cdots,i_{N-1}=1}^n\Bigg(\sum_{p\in\mathcal{P}_N}(H_{p(1)})_{k,i_1}\cdots (H_{p(N)})_{i_{N-1},\ell}\Bigg)f[\lambda_k,\lambda_{i_1},\cdots,\lambda_{i_{N-1}},\lambda_\ell],
\end{equation}
```
where $U=(v_1,\cdots,v_n)$.

## Array operations
Use array operations to efficiently compute the Fréchet derivative. For simplicity, just consider the no permutation case and define
```math
(F_N)_{kℓ}:=∑_{i_1,⋯,i_{N-1}=1}^n(H_1)_{k,i_1}⋯ (H_N)_{i_{N-1},ℓ}Λ^{0,1,…,N-1,N}_{k,i_1,…,i_{N-1},ℓ},
```
where $Λ^{0,…,N}_{i_0,…,i_N} := f[λ_{i_0},⋯,λ_{i_N}]$. It is immediately to obtain that 
```math
F_1 =  H_1 ∘ Λ^{0,1},
```
and 
```math
F_2 = ∑_{i=1}^n (\mathfrak{H}^{1,2} ∘ Λ^{0,1,2})_{:,i,:},
```
 with $\mathfrak{H}^{1,2}_{:,i,:} := (H_1)_{:,i}(H_2)_{i,:}$. For $N ≥ 3$, first compute 
```math
\mathfrak{F}^{0,2,…,N}_{:,:,j_3,…,j_N} := ∑_{i=1}^n (\mathfrak{H}^{1,2} ∘ Λ^{0,1,…,N}_{:,:,:,j_3,…,j_N})_{:,i,:}
```
and permute such that the $0$-dimension is at the end $\mathfrak{F}^{2,…,N,0}$. Then for $M ≥ 3$, there is the recursion
```math
\mathfrak{F}^{M,…,N,0}_{:,j_M,…,j_N} = ∑_{i=1}^N (H_M ∘ \mathfrak{F}^{M-1,…,N,0}_{:,:,j_M,…,j_N})_{i,:}
```
and $F_N = (\mathfrak{F}^{N,0})^T$.

## Adjoint
Denote the first order Fréchet derivative by $\dot{F}= U(\Lambda^{0,1} ∘ U^* \dot{A} U) U^*$, by
```math 
\begin{align*}
  {\rm Re}\langle\overline{F},\dot{F}\rangle&={\rm Re}\langle\overline{F}, U(\Lambda^{0,1} ∘ U^* \dot{A} U) U^*\rangle\\
    &={\rm Re}\langle U^*\overline{F} U,\Lambda^{0,1} ∘ U^* \dot{A} U\rangle\\
    &={\rm Re}\langle U^*\overline{F} U ∘(\Lambda^{0,1})^*,  U^* \dot{A} U\rangle\\
    &={\rm Re}\langle U( U^*\overline{F} U ∘(\Lambda^{0,1})^*) U^*, \dot{A}\rangle\\
    &={\rm Re}\langle\overline{A},\dot{A}\rangle,
\end{align*}
```
we have the adjoint 
```math
\overline{A} =  U( U^*\overline{F} U ∘(\Lambda^{0,1})^*) U^*.
```



