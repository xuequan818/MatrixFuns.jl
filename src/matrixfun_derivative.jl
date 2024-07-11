# Generate the divided difference tensor 
# DF = f[λ_i0, λ_i1, ..., λ_im]
# By the symmetry of the divided difference, 
# we just calculate the irreducible vals
function DD_tensor(f::Function, eigs::Vector{T}, dim::Int64) where {T}
    N = length(eigs)
    DF_sym_index = find_full_indices(N, dim)
    eigs_take(ind) = map(x -> eigs[x], ind)
    DF_sym_val = @. divided_difference(f, eigs_take(DF_sym_index))
    DF = SymmetricTensor(DF_sym_val, Val(N), Val(dim))
end

function frechet_matrix_fun(f::Function, eigs::Vector{T}, Ψ::Matrix{T}, 
                            E::Vector{Matrix{T}}) where {T}
    N = length(eigs)
    order = length(E)
    N0 = N^(order - 2)
    DF = DD_tensor(f, eigs, order + 1)
    DF = reshape(DF, N, N, N, N0)
    E = map(x -> Ψ' * x * Ψ, E)

    if order == 1
        @. DF = DF .* E[1]
        return Ψ * DF * Ψ'
    end

    pert = collect(permutations(1:order))
    val = zeros(T, N, N)
    Einit = zeros(T, N, N, N)
    for p in pert
        @views Ep = E[p]

        @. Einit = zero(T)
        @views for i = 1:N
            Einit[:, i, :] = Ep[1][:, i] * Ep[2][i, :]'
        end

        EF = zeros(T, N, N, N0)
        @views for i = 1:N0
            EF[:, :, i] = dropdims(sum(Einit .* DF[:, :, :, i], dims=2); dims=2)
        end
        EF = reshape(EF, ntuple(x -> N, order))
        EF = permutedims(EF, tuple(2:order..., 1))

        for k = 3:order
            Nk = N^(order + 1 - k)
            DFk = reshape(EF, N, N, Nk)
            EF = zeros(T, N, Nk)
            @views for i = 1:Nk
                EF[:, i] = sum(Ep[k] .* DFk[:, :, i], dims=1)'
            end
        end
        val += EF
    end

    Ψ * val * Ψ'
end