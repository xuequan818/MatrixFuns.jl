"""
    div_diff(f, x::Vector; kwargs...)
    div_diff(f, x::Tuple; kwargs...)
    div_diff(f, x...; kwargs...)

Return the divided difference `f[x_0,x_1,...,x_n]`, assuming `f` is called as `f(x)`. 
"""
@inline function div_diff(f::Function, x::Vector{T}; 
						  kwargs...) where {T<:Real}
	N = length(x)
	@assert N > 1

    sort!(x)
    DD_table = Bidiagonal{T}(x, ones(T, N-1), :U)
    fd = mat_fun(f, DD_table; kwargs...)
    @assert istriu(fd) || norm(tril(fd, -1)) < sqrt(eps(T))

	return fd[1, N]
end
@inline div_diff(f::F, x::Tuple; kwargs...) where {F} = div_diff(f, collect(x); kwargs...)
@inline div_diff(f::F, x...; kwargs...) where {F} = div_diff(f, tuple(x...); kwargs...)
