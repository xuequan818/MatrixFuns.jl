## Implement the atomic block evaluation(Algorithm 4.4).
function atomic_block_fun!(f::Function, F::AbstractMatrix,
                           A::AbstractMatrix, ::Val{N},
                           max_deg::Integer, tol) where {N}
    @assert istriu(A)

    n = size(A, 1)
    Λ = diag(A)
    σ = sum(Λ) / n

    # If f(x) = c, f(A) = c*I
    checkconstf, fσ = constant(f, σ)
    if checkconstf
        copyto!(F, fσ*I)
        return F
    end

    ft = taylor_coeffs(f, σ, max_deg)
    if !isfinite(ft[end])
        max_iter = findfirst(!isfinite, ft) - 2
    else
        max_iter = findlast(!iszero, ft) - 1
    end

    copyto!(F, ft[1]*I)
    M = copy(A - σ * I)
    P0 = copy(M)
    P1 = similar(M)
    sflag = 1
    sext = max(1, ceil(Int,log(cbrt(eps(typeof(tol))),tol))) * n
    Δs = nothing
    μ = compute_μ(A; n)
    for s = 1:max_iter
        mul!(F, ft[s+1], P0, true, true) # F += ft[s] * P0
        normF = norm(F, Inf)
        mul!(P1, P0, M)

        # check convergence
        if abs(ft[s+1]) * norm(P0, Inf) ≤ tol * normF
            if isone(sflag)
                rem_iter = max_iter - s
                (rem_iter < sext) && (sext = rem_iter)
                Δs = compute_Δ(f, Λ, s, n, sext)
            end
            μ * Δs[sflag] * norm(P1, Inf) ≤ tol * normF && break
            sflag += 1
        end
        (sflag == sext+2) && (sflag = 1)
        copy!(P0, P1)
    end

    return F
end

function atomic_block_fun!(f::Function, F::AbstractMatrix,
                           A::AbstractMatrix, ::Val{1}, max_deg::Integer, tol)
    elem_fun!(f, F, A)
end

function atomic_block_fun!(f::Function, F::AbstractMatrix, 
                           A::AbstractMatrix; max_deg=250,
                           tol_tay=cbrt(eps())^2,
                           checknative=native(f))
    if checknative
        fA = f(A)
        typeof(fA) <: Number ? copyto!(F, fA*I) : copy!(F, fA)
    else
        atomic_block_fun!(f, F, A, Val(size(A,1)), max_deg, tol_tay)               
    end
    F
end
function atomic_block_fun(f, A; kwargs...) 
    T = promote_type(eltype(A), typeof(f(A[1])))
    F = fill!(similar(A,T), 0)
    atomic_block_fun!(f, F, A; kwargs...)
end

# Compute Taylor coefficients of `f` at `x0` for orders from `ordl` to `ordr`. 
@inline function taylor_coeffs(f::Function, x0::Number, ordr::Integer, ordl::Integer=0)
    @assert 0 ≤ ordl ≤ ordr
    try
        ordr ≤ 500 ? coef_by_ts(f, x0, ordr, ordl) : coef_by_arb(f, x0, ordr, ordl)
    catch e
        if isa(e, MethodError)
            coef_by_arb(f, x0, ordr, ordl)
        else
            throw(e)
        end
    end
end

# ordl = 0
@inline function coef_by_ts(f::Function, x0::Number, 
                            ordr::Integer, ordl::Val{0})
    f(x0 + Taylor1(typeof(x0), ordr)).coeffs
end
# ordl ≥ 1
@inline function coef_by_ts(f::Function, x0::Number, 
                            ordr::Integer, ordl::Val{N}) where {N}
    f(x0 + Taylor1(typeof(x0), ordr)).coeffs[N+1:end]
end
@inline function coef_by_ts(f::Function, x0::Number, ordr::Integer, ordl::Integer)
    coef_by_ts(f, x0, ordr, Val(ordl))
end

# Special functions
@inline function coef_by_ts(f::typeof(erf), x0::Number, ordr::Integer, ordl::Integer)
    coef_by_ts_special_fun(f, x -> 2/sqrt(pi)*exp(-x^2), x0, ordr, Val(ordl))
end
@inline function coef_by_ts(f::typeof(erfc), x0::Number, ordr::Integer, ordl::Integer)
    coef_by_ts_special_fun(f, x -> -2/sqrt(pi)*exp(-x^2), x0, ordr, Val(ordl))
end

# ordl = 0
@inline function coef_by_ts_special_fun(f::Function, df::Function, x0::Number, 
                                        ordr::Integer, ordl::Val{0})
    c0 = f(x0)
    iszero(ordr) && return [c0]
    isone(ordr)  && return [c0, df(x0)]

    cn = df(x0 + Taylor1(typeof(x0), ordr-1)).coeffs
    @. cn /= 1:ordr
    pushfirst!(cn, c0)
    return cn
end
# ordl ≥ 1
@inline function coef_by_ts_special_fun(f::Function, df::Function, x0::Number,
                                        ordr::Integer, ordl::Val{N}) where {N}
    isone(ordr) && return [df(x0)]
    cn = df(x0 + Taylor1(typeof(x0), ordr-1)).coeffs[N:end]
    @. cn /= N:ordr
end

@inline function coef_by_arb(f::Function, x0::Number, ordr::Integer, ordl::Integer)
    farb = f(ArbSeries((x0, 1), degree=ordr))
    T = float(typeof(x0))
    c = map(i->convert(T, farb[i]), ordl:ordr)
end

# Check `f` is a constant function
const TAYLOR0 = Taylor1(0)

function constant(f, x)
    try 
        fx = f(x + TAYLOR0)
        typeof(fx) <: Taylor1 ? (false, nothing) : (true, fx)
    catch
        (false, nothing)
    end
end

# compute μ = ‖y‖_∞, where y solves (I - |N|)y = e 
# and N is the strictly upper triangular of A
function compute_μ(A::AbstractMatrix{TT}; n=size(A, 1)) where {TT}
    T = real(TT)
    As = triu(A, 1)
    @. As = -abs(As)
    As = convert(Matrix{T}, As)
    mul!(As, true, I, true, true)
    y = LAPACK.trtrs!('U', 'N', 'N', As, ones(n))
    μ = norm(y, Inf)
end

# Δ = max_{λ,0≤r≤n-1}|f^{(s+r+1)}(λ)|/r!
#   = max_{λ,0≤r≤n-1}(|f^{(s+r+1)}(λ)|/(s+r+1)!)*((s+r+1)!/(r+1)!)
function compute_Δ(f::Function, Λ::AbstractVector, s, n, m=10)
    T = eltype(real(f(Λ[1])))
    dfmat = zeros(T, n + m, length(Λ))
    # compute th taylor coefficients
    @views for (i, λi) in enumerate(Λ)
        dfmat[:, i] = abs.(taylor_coeffs(exp, λi, s + m + n, s + 1))
    end
    # normalization
    @views for r = 0:n+m-1
        @. dfmat[r+1, :] = x_factorial(dfmat[r+1, :], s + r + 1, r + 1)
    end
    # find the maximum values
    @views Δ1 = maximum(dfmat[1:n, :])
    @views Δs = maximum!(ones(m), dfmat[n+1:end, :])
    for (j, Δj) in enumerate(Δs)
        Δj < Δ1 && (Δs[j] = Δ1)
    end

    pushfirst!(Δs, Δ1)
end

# Accurately compute `x * n! / m!` for large n
@inline x_factorial(x, n, m) = exp(log(x) - logfactorial(m) + logfactorial(n))

## Compute f.(A)
function elem_fun!(f::Function, F::AbstractArray, A::AbstractArray)
    @assert size(F) == size(A)
    try
        @. F = f(A)
    catch e
        if isa(e, DomainError)
            (eltype(F) <: Complex) || (F = complex.(F))
            @. F = f(complex(A))
        else
            throw(e)
        end
    end
end

function elem_fun(f, A)
    T = promote_type(eltype(A), typeof(f(A[1])))
    F = fill!(similar(A,T), 0)
    elem_fun!(f, F, A)
end
