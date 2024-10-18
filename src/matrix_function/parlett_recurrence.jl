# Implemention of the parlett recurrence in      
# DOI: https://doi.org/10.1017/S0962492910000036 

# Implement the standard Parlett recurrence (Algorithm 4.2).
# This is for the case where all eigenvalues are separated from each other.
function parlett_recurrence(f::Function, T::AbstractMatrix, Λ)
    @assert istriu(T)
    n = size(T, 1)

    F = diagm(dot_fun(f, Λ))
    for j = 2:n, i = j-1:-1:1
        Y = T[i, j] * (F[i, i] - F[j, j]) 
        for k = i+1:j-1
            Y += (F[i, k] * T[k, j] - T[i, k] * F[k, j])
        end
        F[i, j] = Y / (T[i, i] - T[j, j])
    end
 
    return F
end

# Implement the block Parlett recurrence (Algorithm 4.3).
# This is for the case where there are eigenvalues are close.
function block_parlett_recurrence(f::Function, T::AbstractMatrix,
   	 							  block::Vector{UnitRange{Int}}; 
								  kwargs...)
    @assert istriu(T)
	
    F = fill!(similar(T, typeof(f(T[1]))), 0)
    @views for (j, jb) in enumerate(block)
        Tjj = T[jb, jb]
        F[jb, jb] = atomic_block_fun(f, Tjj; kwargs...)
        for i = j-1:-1:1
            ib = block[i]
            Tii = T[ib, ib]
            F[ib, ib] = atomic_block_fun(f, Tii; kwargs...)
            Y = F[ib, ib] * T[ib, jb] - T[ib, jb] * F[jb, jb]
            for k = i+1:j-1
                kb = block[k]
                Y += (F[ib, kb] * T[kb, jb] - T[ib, kb] * F[kb, jb])
            end

            # solve Tii*Fij + Fij*(-Tjj) + (-Y) = 0
            if length(ib) == 1 && length(jb) == 1
                F[ib, jb] = Y ./ (Tii - Tjj)
            else
                F[ib, jb] = sylvester(Tii, -Tjj, -Y)
            end
        end
    end

    return F
end
