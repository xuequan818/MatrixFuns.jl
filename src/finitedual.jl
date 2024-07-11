import Base: +, -, *, /, ^, >, <, >=, <=, ==
import Base: exp, log
import Base: sin, cos

##############
# FiniteDual #
##############
using LinearAlgebra
struct FiniteDual{T} <: Number
	points_num::Int
    DD_table::UpperTriangular{T}
end

################
# Constructors #
################
@inline FiniteDual(x::UpperTriangular) = FiniteDual(size(x,1), x)

@inline function FiniteDual(x::AbstractMatrix)
	if !istriu(x)
		error("Cannot create a divided difference table")
    end
	
	FiniteDual(UpperTriangular(x))
end

@inline FiniteDual(x::Vector{T}) where {T<:Number} = FiniteDual(Bidiagonal(x, ones((length(x) - 1)), :U))

@inline FiniteDual(x::T, points_num::Int) where {T<:Number}= FiniteDual(x * I(points_num))

#####################
# Generic Functions #
#####################

const AMBIGUOUS_TYPES = (AbstractFloat, Irrational, Integer, Rational, Real)

@inline function build_op_table(x::FiniteDual{T}, op_val::BitVector) where {T} 
    N = x.points_num
    points_val = diag(x.DD_table)

    if sum(op_val) in (0, N)
        return FiniteDual(T(op_val[1]), N) 
	else
        op_table = diagm(T.(op_val))
        for k = 1:N-1
            indk = diagind(op_table, k)
            valk = diag(op_table, k - 1)
            valkdiff = valk[2:end] - valk[1:end-1]
            ptsdiff = points_val[k+1:N] - points_val[1:N-k]
            diffBit = (valkdiff .!= 0)
            op_table[indk[diffBit]] = valkdiff[diffBit] ./ ptsdiff[diffBit]
        end
        return FiniteDual(op_table)
	end
end

for op in [:>, :<, :(==), :(>=), :(<=)]
    for R in AMBIGUOUS_TYPES
        @eval @inline $op(x::FiniteDual, y::$R) = build_op_table(x, $op.(diag(x.DD_table), y))
		@eval @inline $op(x::$R, y::FiniteDual) = build_op_table(y, $op.(x,diag(y.DD_table)))
    end
end

###################################
# General Mathematical Operations #
###################################

#################
# Special Cases #
#################

# +/- #
#-----#

for op in (:+, :-)
	@eval @inline function $op(x::FiniteDual, y::FiniteDual) 
		@assert x.points_num == y.points_num
		FiniteDual($op(x.DD_table,y.DD_table))
	end
	@eval @inline $op(x::FiniteDual, y::Number) = FiniteDual($op(x.DD_table, y * I))
	@eval @inline $op(x::Number, y::FiniteDual) = FiniteDual($op(x * I, y.DD_table))
end
@inline -(x::FiniteDual) = FiniteDual(-x.DD_table)

# * #
#---#

@inline function *(x::FiniteDual, y::FiniteDual)
    @assert x.points_num == y.points_num
    FiniteDual(x.DD_table * y.DD_table)
end
@inline *(x::FiniteDual, y::Number) = FiniteDual(y * x.DD_table)
@inline *(x::Number, y::FiniteDual) = y * x

# / #
#---#

@inline function /(x::FiniteDual, y::FiniteDual)
    @assert x.points_num == y.points_num
    FiniteDual(x.DD_table / y.DD_table)
end
@inline /(x::FiniteDual, y::Number) = FiniteDual(x.DD_table / y)
@inline /(x::Number, y::FiniteDual) = FiniteDual(x * I / y.DD_table)

@inline exp(x::FiniteDual) = FiniteDual(exp(x.DD_table))
@inline log(x::FiniteDual) = FiniteDual(Array(exp(x.DD_table)))
@inline ^(x::FiniteDual, y::FiniteDual) = exp(y * log(x))
@inline ^(x::FiniteDual, y::Int) = FiniteDual(^(x.DD_table, y))
@inline ^(x::FiniteDual, y::Number) = FiniteDual(^(x.DD_table, y))
@inline ^(x::Number, y::FiniteDual) = exp(y * log(x))

for func in (:log, :sin, :cos)
    @eval @inline $func(x::FiniteDual) = FiniteDual($func(Array(x.DD_table)))
end
