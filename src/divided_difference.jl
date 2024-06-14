@inline function DD_by_finite_dual(f::Function, xl, xr)
    FD = FiniteDual(xl, xr, one(xl))
    return extract_DD(f(FD))
end

@inline extract_DD(x::FiniteDual) = x.l, x.r, x.Δ
@inline extract_DD(x::Number) = x, x, zero(x)

@inline function divided_difference(f::Function, x, y)
	fx1, fy1, Δ1 = DD_by_finite_dual(f, x, y) # branches of x, can trust fx1 but not fy1
    fy2, fx2, Δ2 = DD_by_finite_dual(f, y, x) # branches of y, can trust fy2 but not fx2
	
    # are x and y in the same branch?
    if fx1 == fx2 && fy1 == fy2
        # same branch
        return Δ1
    else
        return (fy2 - fx1) / (y - x)
    end
end

@inline function divided_difference(f::Function, x::NTuple{N,V}) where {N,V<:Real}
    @assert N > 1 "DD_by_finite_dual only supports more than two x-points"

    df(xl::NTuple{1,V}, xr) = divided_difference(f, xl[1], xr)
    _dfuncs = Function[df]
    for i = 3:N
        iN = i - 1
        dfi(xl::NTuple{iN,V}, xr) where {iN} = divided_difference(x -> _dfuncs[i-2](xl[1:end-1], x), xl[end], xr)
        push!(_dfuncs, dfi)
    end

    return _dfuncs[end](x[1:end-1], x[end])
end
@inline function divided_difference(f::Function, x...)
    V = promote_type(typeof.(x)...)
    N = length(x)
    x = convert(NTuple{N,V}, tuple(x...))
    divided_difference(f, x)
end