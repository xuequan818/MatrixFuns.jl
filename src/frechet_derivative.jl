"""
    mat_fun_frechet(f, eigs, Ψ::AbstractMatrix, h::Vector{AbstractMatrix}; kwargs...)
    mat_fun_frechet(f, H::AbstractMatrix, h::Vector{AbstractMatrix}; kwargs...)

Return the n-th order Fréchet derivative `d^nf(H)h[1]…h[n]`, assuming `f` is called as `f(x)`.
See `mat_fun` for a description of possible keyword arguments.
"""
@inline function mat_fun_frechet(f::Function, eigs::AbstractVector{<:Real},
                                 Ψ::AbstractMatrix, h::Vector{TT};
                                 kwargs...) where {TT<:AbstractMatrix}
    N = length(eigs)
    order = length(h)
    DD_F = DD_tensor(f, eigs, order; kwargs...)
    h = map(x -> inv(Ψ) * x * Ψ, h)

    # F_1 =  h_1 ∘ Λ^{0,1}
    if order == 1
        TD, Th = eltype(DD_F), eltype(TT)
        V = promote_type(TD, Th)
        if TD != V
            DD_F = V.(DD_F)
        end
        @. DD_F *= h[1]
        return Ψ * DD_F * inv(Ψ)
    end

    N0 = N^(order - 2)
    DD_F = reshape(DD_F, N, N, N, N0)

    T = promote_type(eltype(Ψ), eltype(TT))
    val = zeros(T, N, N)
    h12 = zeros(T, N, N, N)

    # loop for the permutations
    pert = collect(permutations(1:order))
    hFinit = zeros(T, N, N, N0)
    @views for p in pert
        hp = h[p]
        # compute {h}^{1,2}_{i,k,j} := (h_1)_{i,k}*(h_2)_{k,j}
        fill!(h12, 0)
        for l = 1:N
            for i = 1:N, j = 1:N
                @inbounds h12[i, l, j] = hp[1][i, l] * hp[2][l, j]
            end
        end

        # compute {F}^{0,2,…,n}_{:,:,j_3,…,j_n} := ∑_{i=1}^N ({h}^{1,2} ∘ Λ^{0,1,…,n}_{:,:,:,j_3,…,j_n})_{:,i,:}
        for i = 1:N0
            tensormultsumdims2!(hFinit[:, :, i], h12, DD_F[:, :, :, i])
        end
        hF = reshape(hFinit, ntuple(x -> N, order))
        hF = permutedims(hF, tuple(2:order..., 1))

        # compute the recursion 
        # {F}^{m,…,n,0}_{:,j_m,…,j_n} = ∑_{i=1}^N (h_m ∘ {F}^{m-1,…,n,0}_{:,:,j_m,…,j_n})_{i,:}
        for k = 3:order
            Nk = N^(order + 1 - k)
            DD_Fk = reshape(hF, N, N, Nk)
            hF = zeros(T, N, Nk)
            for i = 1:Nk
                matmultsumdims1!(hF[:, i], hp[k], DD_Fk[:, :, i])
            end
        end

        mul!(val, true, hF, true, true) # val += hF
    end

    Ψ * transpose(val) * inv(Ψ)
end

@inline function mat_fun_frechet(f::Function, H::AbstractMatrix, 
                                 h::Vector{<:AbstractMatrix}; kwargs...)          
    # compute the full eigen decomposition
    eigs, Ψ = eigen(H)
    # compute the frechet derivative by eigen pairs
    mat_fun_frechet(f, eigs, Ψ, h; kwargs...)
end

"""
Generate the divided difference tensor 
`(DD_F)_{i_0,...,i_n} = f[λ_{i_0},..., λ_{i_n}]`.
By the permutation symmetry of the divided difference, 
only need to calculate the irreducible vals.
"""
function DD_tensor(f::Function, eigs::AbstractVector{<:Real}, 
                   order::Integer; kwargs...)
    N = length(eigs)
    dim = order + 1
    eigs_take(ind) = map(x -> eigs[x], ind)

    DD_F_sym_index = find_full_indices(N, dim)
    DD_F_sym_val = @. div_diff(f, eigs_take(DD_F_sym_index); kwargs...)
    DD_F = SymmetricTensor(DD_F_sym_val, Val(N), Val(dim))

    return Array(DD_F)
end

# Calculate the tensor-tensor(3d) elementwise product and 
# do the in place sum along the second dimension
function tensormultsumdims2!(result, A, B)
    @assert ndims(A) == 3
    is, js, ks = axes(A)
    if (is != axes(result, 1)) || (ks != axes(result, 2)) || ((is, js, ks) != axes(B))
        error("mismatched array sizes")
    end
    for k in ks
        for j in js
            for i in is
                @inbounds c, r = A[i, j, k], B[i, j, k]
                @inbounds result[i, k] = (j == first(js)) ? c * r : muladd(c, r, result[i, k])
            end
        end
    end
end

# Calculate the matrix-matrix elementwise product and 
# do the in place sum along the first dimension
function matmultsumdims1!(result, A, B)
    @assert ndims(A) == 2
    is, js = axes(A)
    if (js != eachindex(result)) || ((is, js) != axes(B))
        error("mismatched array sizes")
    end
    for j in js
        for i in is
            @inbounds c, r = A[i, j], B[i, j]
            @inbounds result[j] = (i == first(is)) ? c * r : muladd(c, r, result[j])
        end
    end
end
