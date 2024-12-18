using ChainRulesCore

function ChainRulesCore.frule((_, Δf, ΔA), ::typeof(mat_fun), f::Function, A::Matrix{<:RealOrComplex}; kwargs...)
    @assert ishermitian(A)

    Ω = mat_fun(f, A; kwargs...)
    ∂Ω = mat_fun_frechet(f, A, [ΔA]; kwargs...)

    return (Ω, ∂Ω)
end

function ChainRulesCore.rrule(::typeof(mat_fun), f::Function, A::Matrix{<:RealOrComplex}; kwargs...)
	@assert ishermitian(A)

    Ω = mat_fun(f, A; kwargs...)
    function mat_fun_pullback(ΔΩ)
        ∂f = NoTangent()
        ∂A = mat_fun_frechet_adjoint(f, A, ΔΩ; kwargs...)
        return (NoTangent(), ∂f, ∂A)
	end
	
    return Ω, mat_fun_pullback
end

function mat_fun_frechet_adjoint(f::Function, A::Matrix{<:RealOrComplex}, h; kwargs...)
    eigs, Ψ = eigen(A)
    ad_DF = DD_tensor(f, eigs, 1; kwargs...)'
    invΨ = inv(Ψ)
    h = invΨ * h * Ψ

    TD, Th = eltype(ad_DF), eltype(h)
    V = promote_type(TD, Th)
    if TD != V
        ad_DF = V.(ad_DF)
    end
    @. ad_DF *= h

    return Ψ * ad_DF * invΨ
end
