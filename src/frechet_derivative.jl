"""
    mat_fun_frechet(f, eigs, Ψ::AbstractMatrix, h::Vector{AbstractMatrix})
    mat_fun_frechet(f, H::AbstractMatrix, h::Vector{AbstractMatrix})

Return the n-th order Fréchet derivative `d^nf(H)h[1]…h[n]`, assuming `f` is called as `f(x)`.
"""
@inline function mat_fun_frechet(f::Function, eigs::Vector{Float64},
                                 Ψ::AbstractMatrix, h::Vector{V}; 
                                 kwargs...) where {V<:AbstractMatrix}
    N = length(eigs)
    order = length(h)
    DD_F = DD_tensor(f, eigs, order; kwargs...)
    map!(x -> inv(Ψ) * x * Ψ, h, h)

    # F_1 =  h_1 ∘ Λ^{0,1}
    if order == 1
        return Ψ * (h[1] .* DD_F) * inv(Ψ)
    end

    N0 = N^(order - 2)
    DD_F = reshape(DD_F, N, N, N, N0)

    T = promote_type(eltype(Ψ), eltype(V))
    val = zeros(T, N, N)
    h12 = zeros(T, N, N, N)

    # loop for the permutations
    pert = collect(permutations(1:order))
    @views for p in pert 
        hp = h[p]
        # compute {h}^{1,2}_{:,i,:} := (h_1)_{:,i}(h_2)_{i,:}
        @. h12 = zero(T)
        for i = 1:N
            h12[:, i, :] = hp[1][:, i] * transpose(hp[2][i, :])
        end

        # compute {F}^{0,2,…,n}_{:,:,j_3,…,j_n} := ∑_{i=1}^N ({h}^{1,2} ∘ Λ^{0,1,…,n}_{:,:,:,j_3,…,j_n})_{:,i,:}
        hF = zeros(T, N, N, N0)
        for i = 1:N0
            hF[:, :, i] = dropdims(sum(h12 .* DD_F[:, :, :, i], dims=2); dims=2)
        end
        hF = reshape(hF, ntuple(x -> N, order))
        hF = permutedims(hF, tuple(2:order..., 1))

        # compute the recursion 
        # {F}^{m,…,n,0}_{:,j_m,…,j_n} = ∑_{i=1}^N (h_m ∘ {F}^{m-1,…,n,0}_{:,:,j_m,…,j_n})_{i,:}
        for k = 3:order
            Nk = N^(order + 1 - k)
            DD_Fk = reshape(hF, N, N, Nk)
            hF = zeros(T, N, Nk)
            for i = 1:Nk
                hF[:, i] = dropdims(sum(hp[k] .* DD_Fk[:, :, i], dims=1), dims=1)
            end
        end

        val += hF 
    end

    Ψ * transpose(val) * inv(Ψ)
end

@inline function mat_fun_frechet(f::Function, H::AbstractMatrix, h::Vector{V};                     
                                 kwargs...) where {V<:AbstractMatrix}
    # compute the full eigen decomposition
    # H = Ψ * diagm(eigs) * inv(Ψ)
    eigs, Ψ = eigen(H)

    # compute the frechet derivative by eigen pairs
    mat_fun_frechet(f, eigs, Ψ, h; kwargs...)
end

# Generate the divided difference tensor 
# DD_F = f[λ_i0, λ_i2, ..., λ_in]
# By the permutation symmetry of the divided difference, 
# only need to calculate the irreducible vals
function DD_tensor(f::Function, eigs::Vector{T}, 
                   order::Integer; kwargs...) where {T}
    N = length(eigs)
    dim = order + 1
    eigs_take(ind) = map(x -> eigs[x], ind)

    DD_F_sym_index = find_full_indices(N, dim)
    DD_F_sym_val = @. div_diff(f, eigs_take(DD_F_sym_index); kwargs...)
    DD_F = SymmetricTensor(DD_F_sym_val, Val(N), Val(dim))

    return Array(DD_F)
end