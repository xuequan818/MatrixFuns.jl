@inline function divided_difference(f::Function, x::Vector)
    FD = f(FiniteDual(x))
    return FD.DD_table[1, FD.points_num]
end
@inline divided_difference(f::Function, x::Tuple) = divided_difference(f, collect(x))
@inline divided_difference(f::Function, x...) = divided_difference(f, tuple(x...))