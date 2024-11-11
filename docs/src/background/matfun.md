# Computing Matrix Functions

[Matrix functions](https://www.statlect.com/matrix-algebra/matrix-function) are scalar functions that map $\mathbb{C}^{n\times n}$ to $\mathbb{C}^{n\times n}$. A general approach to compute $f(A)$ for $A\in\mathbb{C}^{n\times n}$ is to use the similarity transformation $A=ZBZ^{-1}$ and then $f(A)=Zf(B)Z^{-1}$, where $f(B)$ is easy to compute. Specially, if $A$ is diagonalizable, i.e., $B={\rm diag}(b_1,\dots,b_n)$, then $f(B)={\rm diag}\big(f(b_1),\dots,f(b_n)\big)$ is trivially computed. However, this approach is numerically unstable when there are errors in evaluating $f(B)$ and the condition number $\kappa(Z)=\|Z\|\|Z^{-1}\|$ is large. 

For general $A$, a robust approach is to employ the Schur decomposition $A=QTQ^*$, where $Q$ is unitary and $T$ is upper triangular, and then calculate $f(T)$ by Parlett recurrence.

## Parlett recurrence
Since $T$ is upper triangular, $F:=f(T)$ is also upper triangular with $f_{ii}=f(t_{ii})$, and $F$ commutes with $T$. The standard [Parlett recurrence](https://www.sciencedirect.com/science/article/pii/0024379576900185) comes from $FT=TF$:
```math
f_{ij} = t_{ij}\frac{f_{ii}-f_{jj}}{t_{ii}-t_{jj}} + \sum_{k=i+1}^{j-1}\frac{f_{ik}t_{kj}-t_{ik}f_{kj}}{t_{ii}-t_{jj}},\qquad i<j.
```
Unfortunately, this recurrence breaks down when $t_{ii}=t_{jj}$ for some $i\neq j$. The [block Parlett recurrence](https://www2.eecs.berkeley.edu/Pubs/TechRpts/1974/28788.html) is further proposed that regards $T$ as a block upper triangular $T=(T_{ij})$, then $F=(F_{ij})$ has the same block structure,
```math
T_{ii}F_{ij}-F_{ij}T_{jj} = F_{ii}T_{ij}-T_{ij}F_{jj} + \Bigg(\sum_{k=i+1}^{j-1}F_{ik}T_{kj}-T_{ik}F_{kj}\Bigg),\qquad i<j.
```
This Sylvester equation is nonsingular when $T_{ii}$ and $T_{jj}$ have no common eigenvalues.

## Schur-Parlett with improvements
The [Schur-Parlett](https://doi.org/10.1137/S0895479802410815) algorithm is inspired by the block Parlett recurrence, which has two key parts: reordering and blocking of the Schur factor $T$, and computation of the atomic block $f(T_{jj})$. Here we will focus on the first part, for more details on the atomic block computation, please see Section 2 of [![DOI](https://img.shields.io/badge/DOI-10.1137/S0895479802410815-blue)](https://doi.org/10.1137/S0895479802410815).


Let $\widetilde{T}=U^*TU=(\widetilde{T}_{ij})$ be the reordered upper triangular matrix, where $U$ is unitary. The splitting strategy requires that the spectra of the diagonal blocks satisfy: 
* $\min\big\{|\lambda -\mu| : \lambda\in \Lambda(\widetilde{T}_{ii}),\, \mu\in \Lambda(\widetilde{T}_{jj}),\, i\neq j\big\}>\delta$; 
* for $\widetilde{T}_{ii}\in\mathbb{C}^{m\times m}$ with $m>1$, $\forall λ \in \widetilde{T}_{ii}$, $\exists μ ∈ \widetilde{T}_{ii}$ and $μ ≠ λ$, s.t. |$λ - μ| ≤ \delta$.

Here, $\delta>0$ is a tolerance. The second condition can easily lead to large blocks, which destabilizes the atomic block computation based on Taylor expansion. Let $\Delta:=\max\big\{|\lambda -\mu| : \lambda,\,\mu\in \Lambda(\widetilde{T}_{ii})\big\}$ be the spread of the block $\widetilde{T}_{ii}$, $N$ be the maximum Taylor series order, and $\alpha$ be the scaling of the Talyor series error. We split the large block with smaller $\delta$ until 
```math
\bigg(\frac{\Delta}{\alpha}\bigg)^{N+1} \leq \frac{\varepsilon}{\delta},
```
where the left side is the Taylor expansion error and the right side is the splitting error. This condition reduces to $\Delta < \alpha$ when $N=\infty$. Note that the scaling $\alpha$ depends on the smoothness of $f$ in the convex sets containing $ \Lambda(\widetilde{T}_{ii})$. 

Additionally, in order to deal with discontinuous functions, we also use a color mapping $\mathfrak{c}: \mathbb{{C}} \to \mathbb{Z}$ so that eigenvalues from different continuous intervals are not split together. For example, consider the Heaviside step function $H(x)=\pmb{1}_{x\geq 0}$, the color mapping can be defined as 
```math
	\mathfrak{c}(x) = \begin{cases} 
a, & x \geq 0, \\ 
b, & x < 0, 
\end{cases}
```
where $a,b\in\mathbb{Z}$ and $a\neq b$.


The splitting strategy maps each eigenvalue $\lambda_i$ of $T$ to an integer $q_i$, $1\leq q_i \leq n$. The remainning problem is to find a series of swaps to convert $q=(q_1,\dots,q_n)$ to a confluent permutation, i.e., any repeated $q_i$ are next to each other. Instead of ordering $q$ in ascending average index, here we only sort $q_i$ that are not confluent in descending order, e.g., 
```math
q=(1,4,2,1,2,3,3,2) \to (1,1,4,2,2,3,3,2) \to(2,2,2,1,1,4,3,3).
```
The reordering operation is implemented by `ordschur!` in `LinearAlgebra`.