module FrechetDerivative

using Test
import MatrixFuns: mat_fun_frechet, DD_tensor
using LinearAlgebra
using Combinatorics

# compute the frechet derivative by looping elements
function frechet_loop(DD_F::AbstractArray, Ψ::AbstractMatrix,
                      h::Vector{TT}) where {TT<:AbstractMatrix}
    N = size(Ψ, 1)
    order = length(h)
    hs = map(x -> inv(Ψ) * x * Ψ, h)

    pert = collect(permutations(1:order))
    T = promote_type(eltype(Ψ), eltype(TT))
    val = zeros(T, N, N)
    kind = ones(Int, order+1)
    ktr = ones(Int, order-1)
    for i = 1:N, j = 1:N
        fill!(ktr, 1)
        kind[1] = i
        kind[end] = j
        @inbounds for k = 1:N^(order-1)
            copyto!(kind, 2, ktr, 1, order-1)
            hval = zero(T)
            for p in pert
                hval += prod(l -> hs[p[l]][kind[l], kind[l+1]], 1:order)
            end
            val[i, j] += hval * DD_F[kind...]

            if length(ktr) > 0
                ktr[1] += 1
                for ℓ = 1:order-2
                    if ktr[ℓ] == N + 1
                        ktr[ℓ] = 1
                        ktr[ℓ+1] += 1
                    end
                end
            end
        end
    end

    return Ψ * val * inv(Ψ)
end

const N = 10
const f(x) = 1/(1+exp(x))

@testset "Real symmetric matrix" begin
    X = rand(N, N);
    H = 0.5 * (X + X');
	eigs, Ψ = eigen(H);
	for order in 1:4
		h = [rand(N, N) for i = 1:order]
		hs = [0.5 * (x + x') for x in h]
        DD_F = DD_tensor(f, eigs, order)
        fd_test = frechet_loop(DD_F, Ψ, hs)
        fd = mat_fun_frechet(DD_F, Ψ, hs)
        @test isapprox(fd_test, fd)
	end
end

@testset "Complex hermitian matrix" begin
    T = ComplexF64
    X = rand(T, N, N)
    H = 0.5 * (X + X')
    eigs, Ψ = eigen(H)
    for order in 1:4
        h = [rand(T, N, N) for i = 1:order]
        hs = [0.5 * (x + x') for x in h]
        DD_F = DD_tensor(f, eigs, order)
        fd_test = frechet_loop(DD_F, Ψ, hs)
        fd = mat_fun_frechet(DD_F, Ψ, hs)
        @test isapprox(fd_test, fd)
    end
end

@testset "Diagonalizable matrix" begin
    D = rand(N) .- 0.5
    P = rand(N,N)
    H = P * Diagonal(D) * inv(P)
    eigs, Ψ = eigen(H)
    for order in 1:4
        h = [rand(N, N) for i = 1:order]
        hs = [0.5 * (x + x') for x in h]
        DD_F = DD_tensor(f, eigs, order)
        fd_test = frechet_loop(DD_F, Ψ, hs)
        fd = mat_fun_frechet(DD_F, Ψ, hs)
        @test isapprox(fd_test, fd)
    end
end

end