using Test
using MatrixFunDiff
using Combinatorics
using LinearAlgebra

function frechet_naive(f::Function, eigs::Vector{T}, Ψ::Matrix{T}, E::Vector{Matrix{T}}) where {T}
	N = length(eigs)
    order = length(E)
    DF = MatrixFunDiff.DD_tensor(f, eigs, order + 1)
    E = map(x -> Ψ' * x * Ψ, E)

    pert = collect(permutations(1:order))
	val = zeros(T, N, N)
	for i = 1:N, j = 1:N
        ktr = ones(Int, order - 1)
		for k = 1:N^(order-1)
			kind = vcat(i,ktr,j)
			Eval = 0.
			for p in pert
				Eval += prod(l->E[p[l]][kind[l], kind[l+1]], 1:order)
			end
			val[i,j] += Eval * DF[kind...]

            ktr[1] += 1
            for ℓ = 1:order-2
                if ktr[ℓ] == N+1
                    ktr[ℓ] = 1
                    ktr[ℓ+1] += 1
                end
            end
		end
	end

	return Ψ * val * Ψ'
end

f(x) = exp(0.1/(1+x^2))
N = 10
A = rand(N, N)
A = (A + A') / 2
eigs, Ψ = eigen(A)
for order in (2,3,4)
	E = [rand(N, N) for i = 1:order]
	E = [(x + x') / 2 for x in E]

	fd_test = frechet_naive(f, eigs, Ψ, E);
	fd_val = frechet_matrix_fun(f, eigs, Ψ, E);
	@test norm(fd_test-fd_val) < 1e-10
end