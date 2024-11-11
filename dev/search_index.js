var documenterSearchIndex = {"docs":
[{"location":"api/#Matrix-Functions-API","page":"API","title":"Matrix Functions API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"CurrentModule = MatrixFuns","category":"page"},{"location":"api/#Matrix-variable-functions-f(A::AbstractMatrix)::AbstractMatrix","page":"API","title":"Matrix-variable functions f(A::AbstractMatrix)::AbstractMatrix","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MatrixFuns.mat_fun","category":"page"},{"location":"api/#MatrixFuns.mat_fun","page":"API","title":"MatrixFuns.mat_fun","text":"mat_fun(f, A; sep, max_deg, scale, color, ε, checknative)\n\nCompute the matrix-variable function f(A). Return a matrix.\n\nf a general scalar function, and called as f(x).\n\nA an arbitrary square matrix.\n\nsep the initial separation distance, to split the eigenvalues of A into clusters that satisfy (i) min{|λ - μ|: λ ∈ clp, μ ∈ clq, p ≠ q} > sep.\t\t (ii) for clp with |clp| > 1, ∀ λ ∈ clp, ∃ μ ∈ clp and μ ≠ λ, s.t. |λ - μ| ≤ sep.  Here cl_p is the cluster with index p. When sep=Inf, the eigenvalues are only split by color.\n\nmax_deg the maximum Taylor series order for the diagonal blocks computation.\n\ntol_taylor the termination tolerance for evaluating the Taylor series of diagonal blocks.\n\nscale the scaling of the Talyor series error, is used to control the spread of each splitting cluster. When scale=Inf, split the eigenvalues only once with sep.\n\ncolor(x::Number)::Integer an index mapping function.  This is to ensure all eigenvalues within a cluster have the same color index.  By default, the color function maps x to 1.  For discontinuous functions, users can assign a distinct color to each continuous interval.  E.g., for the sign function, users can customize the color mapping using  color = x -> Int(real(sign(x)))\n\nε the default error.\n\nchecknative check f is a Julia native function.\n\n\n\n\n\n","category":"function"},{"location":"api/#Divided-Differences-of-f(x::Real)::Number","page":"API","title":"Divided Differences of f(x::Real)::Number","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MatrixFuns.div_diff","category":"page"},{"location":"api/#MatrixFuns.div_diff","page":"API","title":"MatrixFuns.div_diff","text":"div_diff(f, x::AbstractVector; kwargs...)\ndiv_diff(f, x::Tuple; kwargs...)\ndiv_diff(f, x...; kwargs...)\n\nReturn the divided difference f[x_0,x_1,...,x_n], assuming f is called as f(x).  See mat_fun for a description of possible keyword arguments.\n\n\n\n\n\n","category":"function"},{"location":"api/#Fréchet-derivatives-of-f(A::AbstractMatrix)::AbstractMatrix","page":"API","title":"Fréchet derivatives of f(A::AbstractMatrix)::AbstractMatrix","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MatrixFuns.mat_fun_frechet","category":"page"},{"location":"api/#MatrixFuns.mat_fun_frechet","page":"API","title":"MatrixFuns.mat_fun_frechet","text":"mat_fun_frechet(f, eigs, Ψ::AbstractMatrix, h::Vector{AbstractMatrix}; kwargs...)\nmat_fun_frechet(f, H::AbstractMatrix, h::Vector{AbstractMatrix}; kwargs...)\n\nReturn the n-th order Fréchet derivative d^nf(H)h[1]…h[n], assuming f is called as f(x). See mat_fun for a description of possible keyword arguments.\n\n\n\n\n\n","category":"function"},{"location":"api/#Implementation-details","page":"API","title":"Implementation details","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MatrixFuns.get_splittings\nMatrixFuns.split_cluster\nMatrixFuns.split_by_sep\nMatrixFuns.get_spread\nMatrixFuns.checkspread\nMatrixFuns.reorder_schur\nMatrixFuns.get_swappings\nMatrixFuns.atomic_block_fun!\nMatrixFuns.atomic_block_fun\nMatrixFuns.taylor_coeffs\nMatrixFuns.parlett_recurrence\nMatrixFuns.block_parlett_recurrence\nMatrixFuns.DD_tensor","category":"page"},{"location":"api/#MatrixFuns.get_splittings","page":"API","title":"MatrixFuns.get_splittings","text":"Construct the splitting cluster assignments vector z  (z_i is the index of the splitting cluster for the i-th data points (eigenvalues)). \n\n\n\n\n\n","category":"function"},{"location":"api/#MatrixFuns.split_cluster","page":"API","title":"MatrixFuns.split_cluster","text":"Split the same color mapping index into one cluster\n\n\n\n\n\n","category":"function"},{"location":"api/#MatrixFuns.split_by_sep","page":"API","title":"MatrixFuns.split_by_sep","text":"split_by_sep(pts::AbstractVector{<:Real}, δ::Real)\nsplit_by_sep(pts::AbstractVector{<:Complex}, δ::Real)\n\nSplit the data points (eigenvalues) by separation parameter δ.\n\n\n\n\n\n","category":"function"},{"location":"api/#MatrixFuns.get_spread","page":"API","title":"MatrixFuns.get_spread","text":"Find the spread of the data points (eigenvalues).\n\n\n\n\n\n","category":"function"},{"location":"api/#MatrixFuns.checkspread","page":"API","title":"MatrixFuns.checkspread","text":"Check that the spread error is smaller than the splitting error.\n\n\n\n\n\n","category":"function"},{"location":"api/#MatrixFuns.reorder_schur","page":"API","title":"MatrixFuns.reorder_schur","text":"Reorder the Schur decomposition using ordschur!.\n\n\n\n\n\n","category":"function"},{"location":"api/#MatrixFuns.get_swappings","page":"API","title":"MatrixFuns.get_swappings","text":"Find the swap strategy that converts an integer sequence              into a confluent sequence that the repeated indexes are next to each other. E.g., (1,4,2,1,2,3,3,2) -> (1,1,4,2,2,3,3,2) -> (2,2,2,1,1,4,3,3)\n\n\n\n\n\n","category":"function"},{"location":"api/#MatrixFuns.atomic_block_fun!","page":"API","title":"MatrixFuns.atomic_block_fun!","text":"atomic_block_fun!(f, F, A; max_deg, tol_taylor, checknative)\n\nImplement the atomic block computation f(A), and overwriting the values of F.\n\n\n\n\n\n","category":"function"},{"location":"api/#MatrixFuns.atomic_block_fun","page":"API","title":"MatrixFuns.atomic_block_fun","text":"atomic_block_fun(f, A; kwargs...)\n\nImplement the atomic block computation f(A).\n\n\n\n\n\n","category":"function"},{"location":"api/#MatrixFuns.taylor_coeffs","page":"API","title":"MatrixFuns.taylor_coeffs","text":"taylor_coeffs(f, x0, ordr, ordl::Integer=0)\n\nCompute Taylor coefficients of f at x0 for orders from ordl to ordr. \n\n\n\n\n\n","category":"function"},{"location":"api/#MatrixFuns.parlett_recurrence","page":"API","title":"MatrixFuns.parlett_recurrence","text":"Implement the standard Parlett recurrence. This is for the case where all eigenvalues are separated from each other.\n\n\n\n\n\n","category":"function"},{"location":"api/#MatrixFuns.block_parlett_recurrence","page":"API","title":"MatrixFuns.block_parlett_recurrence","text":"Implement the block Parlett recurrence. This is for the case where there are eigenvalues are close.\n\n\n\n\n\n","category":"function"},{"location":"api/#MatrixFuns.DD_tensor","page":"API","title":"MatrixFuns.DD_tensor","text":"Generate the divided difference tensor  (DD_F)_{i_0,...,i_n} = f[λ_{i_0},..., λ_{i_n}]. By the permutation symmetry of the divided difference,  only need to calculate the irreducible vals.\n\n\n\n\n\n","category":"function"},{"location":"background/divdiff/#Divided-difference","page":"Divided difference","title":"Divided difference","text":"","category":"section"},{"location":"background/divdiff/","page":"Divided difference","title":"Divided difference","text":"Given a function fOmegasubsetmathbbRtomathbbC and n+1 data points x_0 dots x_ninOmega, the higher-order divided differences of f are defined as ","category":"page"},{"location":"background/divdiff/","page":"Divided difference","title":"Divided difference","text":"fx_0 dots x_n\n= begincases \n(fx_0dotsx_n-1-fx_1dotsx_n)(x_0-x_n) rm ifx_0neq x_n  \nfracpartialpartial zfzx_1 dotsx_n-1big_z=x_0  rm ifx_0= x_n\nendcases","category":"page"},{"location":"background/divdiff/","page":"Divided difference","title":"Divided difference","text":"where fx_i = f(x_i).","category":"page"},{"location":"background/divdiff/","page":"Divided difference","title":"Divided difference","text":"The divided-difference table of f is an upper triangular matrix","category":"page"},{"location":"background/divdiff/","page":"Divided difference","title":"Divided difference","text":"T_f(x_0dotsx_n) = \nbeginbmatrix\nfx_0  fx_0x_1 cdots fx_0dotsx_n\nfx_1ddotsvdots\nddots fx_n-1x_n\nfx_n\nendbmatrix","category":"page"},{"location":"background/divdiff/","page":"Divided difference","title":"Divided difference","text":"By the Opitz' formula, T_f can be created by ","category":"page"},{"location":"background/divdiff/","page":"Divided difference","title":"Divided difference","text":"T_f(x_0dotsx_n) = f(J)","category":"page"},{"location":"background/divdiff/","page":"Divided difference","title":"Divided difference","text":"where ","category":"page"},{"location":"background/divdiff/","page":"Divided difference","title":"Divided difference","text":"J = \nbeginbmatrix\nx_01\nx_1ddots\n ddots1\nx_n\nendbmatrix","category":"page"},{"location":"background/divdiff/","page":"Divided difference","title":"Divided difference","text":"Therefore, the accurate computation of the higher-order divided differences fx_0dotsx_n is equivalent to the accurate computation of the matrix function f(J).","category":"page"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"For smooth functions like exp:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> using LinearAlgebra\n\njulia> using MatrixFuns\n\njulia> A = [1 1 0; 0 1.1 1; 0 0 1.2]\n3×3 Matrix{Float64}:\n 1.0  1.0  0.0\n 0.0  1.1  1.0\n 0.0  0.0  1.2\n\njulia> F0 = exp(A) # computed by Julia exponential matrix function\n3×3 Matrix{Float64}:\n 2.71828  2.85884  1.50334\n 0.0      3.00417  3.15951\n 0.0      0.0      3.32012\n\njulia> F1 = mat_fun(exp, A) # computed by Schur-Parlett\n3×3 Matrix{Float64}:\n 2.71828  2.85884  1.50334\n 0.0      3.00417  3.15951\n 0.0      0.0      3.32012\n\njulia> norm(F0-F1) # error of Schur-Parlett\n2.920878410678069e-14\n\njulia> X = 0.05 * Symmetric(A) # generate a symmetric matrix\n3×3 Symmetric{Float64, Matrix{Float64}}:\n 0.05  0.05   0.0\n 0.05  0.055  0.05\n 0.0   0.05   0.06\n\njulia> a = eigvals(X)\n3-element Vector{Float64}:\n -0.01588723439378913\n  0.055000000000000014\n  0.12588723439378913\n\njulia> div_diff(exp, a) # returns the 2nd order divided difference\n0.5284915575854794\n\njulia> hs = map(x->x*X, [1, 2])\n2-element Vector{Symmetric{Float64, Matrix{Float64}}}:\n [0.05 0.05 0.0; 0.05 0.05500000000000001 0.05; 0.0 0.05 0.06]\n [0.1 0.1 0.0; 0.1 0.11000000000000001 0.1; 0.0 0.1 0.12]\n\njulia> mat_fun_frechet(exp, X, hs) # returns the 2nd order Fréchet derivative d^2exp(X)hs_1hs_2\n3×3 Matrix{Float64}:\n 0.0110862   0.0119138  0.00588556\n 0.0119138   0.0181632  0.0130909\n 0.00588556  0.0130909  0.0135867","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"For special functions like error function erf: ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> using SpecialFunctions\n\njulia> mat_fun(erf, A)\n3×3 Matrix{Float64}:\n 0.842701  0.375043  -0.369768\n 0.0       0.880205   0.301089\n 0.0       0.0        0.910314","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"For singular functions, such as Fermi-Dirac functions with temperatures close to 0, users can set a smaller scale to reduce the spread of each block to avoid large Taylor expansion errors near the singularities.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> f(x;μ) = 1/(1+exp(1000*(x-μ))); # Fermi-Dirac function with temperature equal 1e-3.\n\njulia> f1(x) = f(x;μ=1.15);\n\njulia> f2(x) = f(x;μ=0.05);\n\njulia> ## returns the wrong result for default `scale`\n\njulia> mat_fun(f1, A) \n3×3 Matrix{Float64}:\n 1.0  6.48449e6  -9.72675e7\n 0.0  6.4845e5   -1.2969e7\n 0.0  0.0        -6.48449e5\n\njulia> div_diff(f2, a)\n-5.0306612072297075e147\n\njulia> mat_fun_frechet(f2, X, hs)\n3×3 Matrix{Float64}:\n 1.84628e115  8.22808e114  -2.08877e115\n 2.55222e115  9.24584e114  -3.37442e115\n 1.80227e115  4.61108e114  -2.77086e115\n\njulia> ## returns the correct result for smaller `scale`\n\njulia> scale = 0.01;\n\njulia> mat_fun(f1, A; scale)\n3×3 Matrix{Float64}:\n 1.0  0.0  -50.0\n 0.0  1.0  -10.0\n 0.0  0.0    1.92875e-22\n\njulia> div_diff(f2, a; scale)\n98.17057693049055\n\njulia> mat_fun_frechet(f2, X, hs; scale)\n3×3 Matrix{Float64}:\n  19.7425    1.97425   -19.7425\n   1.97425   0.197425   -1.97425\n -19.7425   -1.97425    19.7425","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"For piecewise functions consisting of continuous intervals, users can customize color function to first split the eigenvalues into different continuous intervals.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> x0 = 0.051 # discontinuous point\n\njulia> g(x) = x < x0 ? sin(x) : (x+1)^2;\n\njulia> color(x) = x < x0 ? 1 : 2; # specify a different color for each continuous interval\n\njulia> B = triu(X)\n3×3 Matrix{Float64}:\n 0.05  0.05   0.0\n 0.0   0.055  0.05\n 0.0   0.0    0.06\n\njulia> mat_fun(g, B; sep=0.001) # reference computed by standard Parlett recurrence\n3×3 Matrix{Float64}:\n 0.0499792  10.6305   -52.6235\n 0.0         1.11302    0.10575\n 0.0         0.0        1.1236\n\njulia> mat_fun(g, B) # returns the wrong result\n3×3 Matrix{Float64}:\n 1.1025  0.10525  0.0025\n 0.0     1.11302  0.10575\n 0.0     0.0      1.1236\n\njulia> mat_fun(g, B; color) # returns the correct result\n3×3 Matrix{Float64}:\n 0.0499792  10.6305   -52.6235\n 0.0         1.11303    0.10575\n 0.0         0.0        1.1236","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Specially, for constant functions, or piecewise functions consisting of constant intervals, users can further set sep=Inf to only split the eigenvalues by color, which will be more efficient.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> heaviside(x) = x < x0 ? 1 : 0 # Heaviside step function\n\njulia> ## set `checknative=false` to avoid the cost of checking function nativity\n\njulia> @time mat_fun(heaviside, B; color, checknative=false) \n  0.000099 seconds (100 allocations: 6.422 KiB)\n3×3 Matrix{Float64}:\n 1.0  -10.0  50.0\n 0.0    0.0   0.0\n 0.0    0.0   0.0\n\njulia> @time mat_fun(heaviside, B; color, sep=Inf, checknative=false)\n  0.000087 seconds (87 allocations: 5.641 KiB)\n3×3 Matrix{Float64}:\n 1.0  -10.0  50.0\n 0.0    0.0   0.0\n 0.0    0.0   0.0","category":"page"},{"location":"background/frechet/#Fréchet-derivative","page":"Fréchet derivative","title":"Fréchet derivative","text":"","category":"section"},{"location":"background/frechet/#Contour-integral-form","page":"Fréchet derivative","title":"Contour integral form","text":"","category":"section"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"Let mathcalH=mathbbC^Ntimes N_rm herm be the vector space of Ntimes N Hermitian matrices. For HinmathcalH and h_1h_ninmathcalH, and f is an n times continuously differentiable function on a subset of mathbbC containing the spectrum of H+t_1h_1+cdots + t_nh_n, the n-th order Fréchet derivative of f(H) is","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"beginalign*\n    rm d^nf(H)h_1cdots h_n = frac12pi ioint_mathcalC f(z) sum_pinmathcalP_n(z-H)^-1h_p(1)(z-H)^-1cdots(z-H)^-1h_p(n)(z-H)^-1 dz\nendalign*","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"where mathcalC is a contour in the complex plane enclosing all the eigenvalues of H, and pinmathcalP_n is an arbitrary permutation of 1cdotsn. This can be proved by induction.","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"For n = 1, we have ","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"beginalign*\n    rm d f(H)h=lim_tto 0 fracf(H+th)-f(H)t=frac12pi ioint_mathcalC f(z)lim_tto 0frac(z-H-th)^-1-(z-H)^-1tdz\nendalign*","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"Note that","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"beginalign*\n    (z-H-th)^-1  = (I-t(z-H)^-1h)^-1(z-H)^-1 = (I+t(z-H)^-1h+O(t^2))(z-H)^-1\nendalign*","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"we have ","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"beginalign*\n    rm d f(H)h=frac12pi ioint_mathcalC f(z)(z-H)^-1h(z-H)^-1dz\nendalign*","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"which satisfies the formula.","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"Assume the n-1-th order derivative satisfies the formula. Then we have the n-th order derivative","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"beginalign*\n    rm d^nf(H)h_1cdots h_n =lim_tto 0 fracrm d^n-1f(H+th_n)h_1cdots h_n-1-rm d^n-1f(H)h_1cdots h_n-1t\n    =frac12pi ioint_mathcalC f(z)lim_tto 0sum_pinmathcalP_n-1frac1tBig((z-H-th_n)^-1h_p(1)cdots h_p(n-1)(z-H-th_n)^-1\n    qquad -(z-H)^-1h_p(1)cdots h_p(n-1)(z-H)^-1Big)dz\nendalign*","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"Similarly, we have","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"beginalign*\n   (z-H-th_n)^-1h_p(1)cdots h_p(n-1)(z-H-th_n)^-1\n    =(I+t(z-H)^-1h_n)(z-H)^-1h_p(1)cdots h_p(n-1)(z-H)^-1(I+th_n(z-H)^-1) + O(t^2)\n    =(z-H)^-1h_p(1)cdots h_p(n-1)(z-H)^-1 \n    quad+ tBig((z-H)^-1h_n(z-H)^-1h_p(1)cdots h_p(n-1)(z-H)^-1\n    qquadquad+(z-H)^-1h_p(1)(z-H)^-1h_ncdots h_p(n-1)(z-H)^-1\n    qquadquad+(z-H)^-1h_p(1)(z-H)^-1h_p(2)cdots h_n(z-H)^-1Big) + O(t^2)\nendalign*","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"Therefore, we can obtain the formula.","category":"page"},{"location":"background/frechet/#Divided-difference-form","page":"Fréchet derivative","title":"Divided difference form","text":"","category":"section"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"Let H=Phi LambdaPhi^-1=sum_i^Nlambda_iphi_iphi_i^-1, where phi_i is the i-th column of Phi and phi_i^-1 is the i-th row of Phi^-1, and assume there is no degeneration. Then for pinmathcalP_n we have","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"beginalign*\n    (z-H)^-1h_p(1)(z-H)^-1cdots(z-H)^-1h_p(n)(z-H)^-1\n    =sum_i_0cdotsi_n=1^Nphi_i_0(h_p(1))_i_0i_1cdots (h_p(n))_i_n-1i_nphi_i_n^-1(z-lambda_i_0)^-1cdots (z-lambda_i_n)^-1\nendalign*","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"where (h_p(k))_ij=phi_i^-1h_p(k)phi_j. We let","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"    (z-lambda_i_0)^-1cdots (z-lambda_i_n)^-1 = sum_k=0^nC_k (z-lambda_i_k)^-1","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"then","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"sum_k=0^nC_k prod_ellneq k(z-lambda_i_ell)=1","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"Let z=lambda_i_k, we can obtain","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"C_k = frac1prod_ellneq k(lambda_i_k-lambda_i_ell)","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"Therefore, we have","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"beginalign*\n    frac12pi ioint_mathcalC f(z)(z-lambda_i_0)^-1cdots (z-lambda_i_n)^-1dz = sum_k=0^nfracf(lambda_i_k)prod_ellneq k(lambda_i_k-lambda_i_ell)=flambda_i_0cdotslambda_i_n\nendalign*","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"Finally, we obtain","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"beginequation\n    rm d^nf(H)h_1cdots h_n =sum_i_0cdotsi_n=1^Nphi_i_0Bigg(sum_pinmathcalP_n(h_p(1))_i_0i_1cdots (h_p(n))_i_n-1i_nBigg)flambda_i_0cdotslambda_i_nphi_i_n^-1\nendequation","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"which can be also written as ","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"beginequation\n(Phi^-1rm d^nf(H)h_1cdots h_nPhi)_kell=sum_i_1cdotsi_n-1=1^NBigg(sum_pinmathcalP_n(h_p(1))_ki_1cdots (h_p(n))_i_n-1ellBigg)flambda_klambda_i_1cdotslambda_i_n-1lambda_ell\nendequation","category":"page"},{"location":"background/frechet/#Array-operations","page":"Fréchet derivative","title":"Array operations","text":"","category":"section"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"Use array operations to efficiently compute the Fréchet derivative. For simplicity, just consider the no permutation case and define","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"(F_n)_kℓ=_i_1i_n-1=1^N(h_1)_ki_1 (h_n)_i_n-1ℓΛ^01n-1n_ki_1i_n-1ℓ","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"where Λ^0n_i_0i_n = fλ_i_0λ_i_n. It is immediately to obtain that ","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"F_1 =  h_1  Λ^01","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"and ","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"F_2 = _i=1^N (mathfrakh^12  Λ^012)_i","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"with mathfrakh^12_i = (h_1)_i(h_2)_i. For n  3, first compute ","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"mathfrakF^02n_j_3j_n = _i=1^N (mathfrakh^12  Λ^01n_j_3j_n)_i","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"and permute such that the 0-dimension is at the end mathfrakF^2n0. Then for m  3, there is the recursion","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"mathfrakF^mn0_j_mj_n = _i=1^N (h_m  mathfrakF^m-1n0_j_mj_n)_i","category":"page"},{"location":"background/frechet/","page":"Fréchet derivative","title":"Fréchet derivative","text":"and F_n = (mathfrakF^n0)^T.","category":"page"},{"location":"#MatrixFuns","page":"Home","title":"MatrixFuns","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A Julia package for computing scalar functions of matrix variables and their Fréchet derivatives. The matrix functions computation (for arbitrary square matrices) is based on the Schur-Parlett algorithm (with improvements). The higher order Fréchet derivatives (for Hermitian matrices) are formulated similarly to the Daleskii-Krein theorem, where the divided differences are calculated accurately by the Opitz' formula. In particular, MatrixFuns supports the computation of discontinuous functions. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"MatrixFuns is currently an unregistered package and therefore needs to be downloaded or cloned to the user's local computer first, and then installed by running","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> cd(\"your-local-path/MatrixFuns.jl\")\n\njulia> using Pkg\n\njulia> Pkg.activate(\".\")\n\njulia> Pkg.instantiate()","category":"page"},{"location":"background/matfun/#Computing-Matrix-Functions","page":"Computing Matrix Functions","title":"Computing Matrix Functions","text":"","category":"section"},{"location":"background/matfun/","page":"Computing Matrix Functions","title":"Computing Matrix Functions","text":"Matrix functions are scalar functions that map mathbbC^ntimes n to mathbbC^ntimes n. A general approach to compute f(A) for AinmathbbC^ntimes n is to use the similarity transformation A=ZBZ^-1 and then f(A)=Zf(B)Z^-1, where f(B) is easy to compute. Specially, if A is diagonalizable, i.e., B=rm diag(b_1dotsb_n), then f(B)=rm diagbig(f(b_1)dotsf(b_n)big) is trivially computed. However, this approach is numerically unstable when there are errors in evaluating f(B) and the condition number kappa(Z)=ZZ^-1 is large. ","category":"page"},{"location":"background/matfun/","page":"Computing Matrix Functions","title":"Computing Matrix Functions","text":"For general A, a robust approach is to employ the Schur decomposition A=QTQ^*, where Q is unitary and T is upper triangular, and then calculate f(T) by Parlett recurrence.","category":"page"},{"location":"background/matfun/#Parlett-recurrence","page":"Computing Matrix Functions","title":"Parlett recurrence","text":"","category":"section"},{"location":"background/matfun/","page":"Computing Matrix Functions","title":"Computing Matrix Functions","text":"Since T is upper triangular, F=f(T) is also upper triangular with f_ii=f(t_ii), and F commutes with T. The standard Parlett recurrence comes from FT=TF:","category":"page"},{"location":"background/matfun/","page":"Computing Matrix Functions","title":"Computing Matrix Functions","text":"f_ij = t_ijfracf_ii-f_jjt_ii-t_jj + sum_k=i+1^j-1fracf_ikt_kj-t_ikf_kjt_ii-t_jjqquad ij","category":"page"},{"location":"background/matfun/","page":"Computing Matrix Functions","title":"Computing Matrix Functions","text":"Unfortunately, this recurrence breaks down when t_ii=t_jj for some ineq j. The block Parlett recurrence is further proposed that regards T as a block upper triangular T=(T_ij), then F=(F_ij) has the same block structure,","category":"page"},{"location":"background/matfun/","page":"Computing Matrix Functions","title":"Computing Matrix Functions","text":"T_iiF_ij-F_ijT_jj = F_iiT_ij-T_ijF_jj + Bigg(sum_k=i+1^j-1F_ikT_kj-T_ikF_kjBigg)qquad ij","category":"page"},{"location":"background/matfun/","page":"Computing Matrix Functions","title":"Computing Matrix Functions","text":"This Sylvester equation is nonsingular when T_ii and T_jj have no common eigenvalues.","category":"page"},{"location":"background/matfun/#Schur-Parlett-with-improvements","page":"Computing Matrix Functions","title":"Schur-Parlett with improvements","text":"","category":"section"},{"location":"background/matfun/","page":"Computing Matrix Functions","title":"Computing Matrix Functions","text":"The Schur-Parlett algorithm is inspired by the block Parlett recurrence, which has two key parts: reordering and blocking of the Schur factor T, and computation of the atomic block f(T_jj). Here we will focus on the first part, for more details on the atomic block computation, please see Section 2 of (Image: DOI).","category":"page"},{"location":"background/matfun/","page":"Computing Matrix Functions","title":"Computing Matrix Functions","text":"Let widetildeT=U^*TU=(widetildeT_ij) be the reordered upper triangular matrix, where U is unitary. The splitting strategy requires that the spectra of the diagonal blocks satisfy: ","category":"page"},{"location":"background/matfun/","page":"Computing Matrix Functions","title":"Computing Matrix Functions","text":"separation between blocks: minbiglambda -mu  lambdain Lambda(widetildeT_ii) muin Lambda(widetildeT_jj) ineq jbigdelta,\nseparation within blocks: for widetildeT_iiinmathbbC^mtimes m with m1, forall λ in widetildeT_ii, exists μ  widetildeT_ii and μ  λ, s.t. |λ - μ  delta.","category":"page"},{"location":"background/matfun/","page":"Computing Matrix Functions","title":"Computing Matrix Functions","text":"Here, delta0 is a splitting tolerance. The second condition can easily lead to large blocks, which destabilizes the atomic block computation based on Taylor expansion. Let Delta=maxbiglambda -mu  lambdamuin Lambda(widetildeT_ii)big be the spread of the block widetildeT_ii, N be the maximum Taylor series order, and alpha be the scaling of the Talyor series error. We split the large block with smaller delta until ","category":"page"},{"location":"background/matfun/","page":"Computing Matrix Functions","title":"Computing Matrix Functions","text":"bigg(fracDeltaalphabigg)^N+1 leq fracvarepsilondelta","category":"page"},{"location":"background/matfun/","page":"Computing Matrix Functions","title":"Computing Matrix Functions","text":"where the left side is the Taylor expansion error and the right side is the splitting error. This condition reduces to Delta  alpha when N=infty. Note that the scaling alpha depends on the smoothness of f in the convex sets containing $ \\Lambda(\\widetilde{T}_{ii})$. ","category":"page"},{"location":"background/matfun/","page":"Computing Matrix Functions","title":"Computing Matrix Functions","text":"Additionally, in order to deal with discontinuous functions, we also use a color mapping mathfrakc mathbbC to mathbbZ so that eigenvalues from different continuous intervals are not split together. For example, consider the Heaviside step function H(x)=pmb1_xgeq 0, the color mapping can be defined as ","category":"page"},{"location":"background/matfun/","page":"Computing Matrix Functions","title":"Computing Matrix Functions","text":"\tmathfrakc(x) = begincases \na  x geq 0  \nb  x  0 \nendcases","category":"page"},{"location":"background/matfun/","page":"Computing Matrix Functions","title":"Computing Matrix Functions","text":"where abinmathbbZ and aneq b.","category":"page"},{"location":"background/matfun/","page":"Computing Matrix Functions","title":"Computing Matrix Functions","text":"The splitting strategy maps each eigenvalue lambda_i of T to an integer q_i, 1leq q_i leq n. The remainning problem is to find a series of swaps to convert q=(q_1dotsq_n) to a confluent permutation, i.e., any repeated q_i are next to each other. Instead of ordering q in ascending average index, here we only sort q_i that are not confluent in descending order, e.g., ","category":"page"},{"location":"background/matfun/","page":"Computing Matrix Functions","title":"Computing Matrix Functions","text":"q=(14212332) to (11422332) to(22211433)","category":"page"},{"location":"background/matfun/","page":"Computing Matrix Functions","title":"Computing Matrix Functions","text":"The reordering operation is implemented by ordschur! in LinearAlgebra.","category":"page"}]
}
