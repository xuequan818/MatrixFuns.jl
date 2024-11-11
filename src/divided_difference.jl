"""
    div_diff(f, x::AbstractVector; kwargs...)
    div_diff(f, x::Tuple; kwargs...)
    div_diff(f, x...; kwargs...)

Return the divided difference `f[x_0,x_1,...,x_n]`, assuming `f` is called as `f(x)`. 
See `mat_fun` for a description of possible keyword arguments.
"""
@inline function div_diff(f::Function, x::AbstractVector{T}; 
						  kwargs...) where {T<:Real}
	N = length(x)
    if isone(N)
        return f(x[1])
    end

	if !issorted(x; rev=true)
        x = sort(x; rev=true)
    end

    if !(typeof(x) <: Vector)
        x = Array(x)
    end

    # Compute the divided difference by Opitz' formula 
    DD_table = Bidiagonal{T}(x, ones(T, N-1), :U)
    fd = mat_fun(f, DD_table; kwargs...)
    @assert istriu(fd) || norm(tril(fd, -1)) < sqrt(eps(T))

    return extract_dd(fd[1,N])
end
@inline div_diff(f, x::Tuple; kwargs...) = div_diff(f, collect(x); kwargs...)
@inline div_diff(f, x...; kwargs...) = div_diff(f, tuple(x...); kwargs...)

function extract_dd(dd::T) where {T<:Complex}
    abs(imag(dd)) < eps(float(real(T))) ? real(dd) : dd
end
extract_dd(dd::Real) = dd
