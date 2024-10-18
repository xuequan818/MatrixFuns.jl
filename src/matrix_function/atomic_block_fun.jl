
function atomic_block_fun(f::Function, A::AbstractMatrix, ::Val{N},
                          max_deg::Int) where {N}
    @assert istriu(A)

    n = size(A, 1)
    tol = cbrt(eps())^2
    Λ = diag(A)
    σ = sum(Λ) / n
    M = copy(A - σ * I)
    
    triu!(A, 1)
    @. A = abs(A)
    μ = norm((I - A) \ ones(n), Inf)

    ft = taylor_coeffs(f, σ, max_deg)
    F = ft[1] * Matrix{eltype(A)}(I, n, n)
    P0 = copy(M)
    P1 = similar(M)
    for s = 1:max_deg
        mul!(F, ft[s+1], P0, true, true) # F += ft[s] * P0
        normF = norm(F, Inf)
        mul!(P1, P0, M)
        if abs(ft[s+1]) * norm(P0, Inf) ≤ tol * normF # check convergence
            # Δ = max_{λ} max_{0≤r≤n-1} |f^{(s+r)} (λ)| / r!
            #   = max_{λ} max_{0≤r≤n-1} (|f^{(s+r)} (λ)| / (s+r)!) * ((s+r)!/r!)
            Δ = 0.0
            for λ in Λ
                dfsr = abs.(taylor_coeffs(f, λ, s + n - 1)[s+1:end])
                for (r, vr) in enumerate(dfsr)
                    dfsr[r] = facm2n(vr, s + r - 1, r - 1)
                end
                Δλ = maximum(dfsr)
                (Δλ > Δ) && (Δ = Δλ)
            end
            μ * Δ * norm(P1, Inf) ≤ tol * normF && break
        end
        copy!(P0, P1)
    end

    return F
end
atomic_block_fun(f::Function, A::AbstractMatrix, ::Val{1}, max_deg) = dot_fun(f, A)

function atomic_block_fun(f::Function, A::AbstractMatrix; 
                          max_deg=100, isnative = native(f))
    isnative ? f(A) : atomic_block_fun(f, A, Val(size(A,1)), max_deg)
end

@inline function taylor_coeffs(f::Function, x0::Number, order::Int)
	try
        order ≤ 500 ? coef_by_ts(f, x0, order) : coef_by_arb(f, x0, order)
	catch e
		if isa(e, MethodError)
			coef_by_arb(f, x0, order)
		else
			throw(e)
		end
	end
end

@inline function coef_by_ts(f::Function, x0::Number, order::Int)
    f(x0 + Taylor1(typeof(x0), order)).coeffs
end

@inline function coef_by_ts(f::typeof(erf), x0::Number, order::Int)
    coef_by_ts_special_fun(f, x -> 2 / sqrt(pi) * exp(-x^2), x0, order)
end

@inline function coef_by_ts(f::typeof(erfc), x0::Number, order::Int)
    coef_by_ts_special_fun(f, x -> -2 / sqrt(pi) * exp(-x^2), x0, order)
end

@inline function coef_by_ts_special_fun(f::Function, df::Function, x0::Number, order::Int)
    c0 = f(x0)
    iszero(order) && return [c0]

    cn = df(x0 + Taylor1(typeof(x0), order - 1)).coeffs
    for (i, ci) in enumerate(cn)
        cn[i] = ci / i
    end
    pushfirst!(cn, c0)
end

@inline function coef_by_arb(f::Function, x0::Number, order::Int)
    farb = f(ArbSeries((x0, 1), degree=order))
    T = typeof(x0)
    c = map(i->convert(T, farb[i]), 0:order)
end

# Accurately compute `x * n! / m!` for large n
@inline facm2n(x, n, m) = exp(log(x) - logfactorial(m) + logfactorial(n))

# Compute f.(A)
function dot_fun!(f::Function, A::AbstractArray)
    try
        @. A = f(A)
    catch e
        if isa(e, DomainError)
            @. A = f(complex(A))
        else
            throw(e)
        end
    end
end
dot_fun(f, A) = dot_fun!(f, copy(A))
