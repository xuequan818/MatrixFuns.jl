# Implemention of the parlett recurrence in      
# DOI: https://doi.org/10.1017/S0962492910000036 

# Implement the standard Parlett recurrence (Algorithm 4.2).
# This is for the case where all eigenvalues are separated from each other.
function parlett_recurrence(f::Function, T::AbstractMatrix)
    @assert istriu(T)
    n = size(T, 1)

    Λ = diag(T)
    F = diagm(elem_fun!(f, Λ, Λ))
    @views for j = 2:n, i = j-1:-1:1
        # TODO: increase numerical stability
        k = i+1:j-1
        Y = T[i, j] * (F[i, i] - F[j, j]) + 
            (dot_noconj(T[k, j], F[i, k]) - dot_noconj(T[i, k], F[k, j]))
        F[i, j] = Y / (T[i, i] - T[j, j])
    end
 
    return F
end

dot_noconj(x::AbstractVector{T}, y::AbstractVector) where {T<:Real} = dot(x,y)
dot_noconj(x::AbstractVector{T}, y::AbstractVector{T}) where {T<:Complex} = BLAS.dotu(x, y)
dot_noconj(x::Number, y::Number) = x*y

# Implement the block Parlett recurrence (Algorithm 4.3).
# This is for the case where there are eigenvalues are close.
function block_parlett_recurrence(f::Function, T::AbstractMatrix,
    block::Vector{UnitRange{Int}};
    kwargs...)
    @assert istriu(T)

    F = fill!(similar(T, typeof(f(T[1]))), 0)
    @views for (j, jb) in enumerate(block)
        Tjj = T[jb, jb]
        Fjj = F[jb, jb]
        atomic_block_fun!(f, Fjj, Tjj; kwargs...)
        for i = j-1:-1:1
            ib = block[i]
            Tii = T[ib, ib]
            Fii = F[ib, ib]
            atomic_block_fun!(f, Fii, Tii; kwargs...)
            # Fij = Fii*Tij - Tij*Fjj
            Tij = T[ib, jb]
            Fij = F[ib, jb]
            mul!(Fij, Fii, Tij)
            mul!(Fij, Tij, Fjj, -1, true)
            # Fij = Fij + (Fik*Tkj - Tik*Fkj)
            for k = i+1:j-1
                kb = block[k]
                mul!(Fij, F[ib, kb], T[kb, jb], 1, true)
                mul!(Fij, T[ib, kb], F[kb, jb], -1, true)
            end

            # solve Tii*Fij - Fij*Tjj = Y
            if length(ib) == 1 && length(jb) == 1
                rdiv!(Fij, Tii[1] - Tjj[1])
            else
                _, scale = LAPACK.trsyl!('N', 'N', Tii, Tjj, Fij, -1)
                rdiv!(Fij, scale)
            end
        end
    end

    return F
end
