using ForwardDiff
##############
# FiniteDual #
##############

struct FiniteDual{V}
    l::V
    r::V
    Δ::V
end

################
# Constructors #
################
@inline function FiniteDual(l::V, r::V, Δ::V) where {V<:Int}
    FiniteDual(float(l), float(r), float(Δ))
end

@inline function FiniteDual(l::A, r::B, Δ::C) where {A,B,C}
    D = promote_type(A, B, C)
    return FiniteDual(convert(D, l), convert(D, r), convert(D, Δ))
end

@inline FiniteDual(l::FiniteDual, r::Number, Δ::Number) = FiniteDual(l, convert(typeof(l), r), convert(typeof(l), Δ))
@inline FiniteDual(l::Number, r::FiniteDual, Δ::Number) = FiniteDual(convert(typeof(r), l), r, convert(typeof(r), Δ))

@inline FiniteDual{V}(x::Number) where {V} = convert(FiniteDual{V}, x)
@inline Base.convert(::Type{FiniteDual{V}}, x::Number) where {V} = FiniteDual(V(x), V(x), zero(V))


####################################
# N-ary Operation Definition Tools #
####################################

macro define_binary_finitedual_op(f, xy_body, x_body, y_body)
    defs = quote
        @inline $(f)(x::$FiniteDual{Txy}, y::$FiniteDual{Txy}) where {Txy} = $xy_body
        #@inline $(f)(x::$FiniteDual{Tx}, y::$FiniteDual{Ty}) where {Tx,Ty} = Ty ≺ Tx ? $x_body : $y_body
    end
    expr = quote
        @inline $(f)(x::$FiniteDual{Tx}, y::Number) where {Tx} = $x_body
        @inline $(f)(x::Number, y::$FiniteDual{Ty}) where {Ty} = $y_body
    end
    append!(defs.args, expr.args)
    return esc(defs)
end

#####################
# Generic Functions #
#####################
@inline Base.zero(x::FiniteDual) = zero(typeof(x))
@inline Base.zero(::Type{FiniteDual{V}}) where {V} = FiniteDual(zero(V), zero(V), zero(V))

@inline Base.one(x::FiniteDual) = one(typeof(x))
@inline Base.one(::Type{FiniteDual{V}}) where {V} = one(V)

for pred in [:isless, :<, :>, :(<=), :(>=)]
    @eval begin
        @define_binary_finitedual_op(
            Base.$(pred),
            $(pred)(x.l, y.l),
            $(pred)(x.l, y),
            $(pred)(y, x.l),
        )
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

@define_binary_finitedual_op(
    Base.:+,
    FiniteDual{Txy}(x.l + y.l, x.r + y.r, x.Δ + y.Δ),
    FiniteDual{Tx}(x.l + y, x.r + y, x.Δ),
    FiniteDual{Ty}(x + y.l, x + y.r, y.Δ)
)

@define_binary_finitedual_op(
    Base.:-,
    FiniteDual{Txy}(x.l - y.l, x.r - y.r, x.Δ - y.Δ),
    FiniteDual{Tx}(x.l - y, x.r - y, x.Δ),
    FiniteDual{Ty}(x - y.l, x - y.r, -y.Δ)
)
@inline Base.:-(x::FiniteDual) = FiniteDual(-x.l, -x.r, -x.Δ)

# * #
#---#

@define_binary_finitedual_op(
    Base.:*,
    FiniteDual{Txy}(x.l * y.l, x.r * y.r, x.Δ * y.l + y.Δ * x.r),
    FiniteDual{Tx}(x.l * y, x.r * y, x.Δ * y),
    FiniteDual{Ty}(x * y.l, x * y.r, x * y.Δ)
)

# / #
#---#

@define_binary_finitedual_op(
    Base.:/,
    FiniteDual{Txy}(x.l / y.l, x.r / y.r, 0.5 * (x.Δ * (y.l + y.r) - y.Δ * (x.l + x.r)) / (y.l * y.r)),
    FiniteDual{Tx}(x.l / y, x.r / y, x.Δ / y),
    FiniteDual{Ty}(x / y.l, x / y.r, -(x * y.Δ) / (y.l * y.r))
)

# exponentiation #
#----------------#

function Base.:^(x::FiniteDual, y::T) where {T<:Number}
    if T <: Int
        vd = zero(typeof(x.l))
        for i = 1:y
            vd += x.l^(y - i) * x.r^(i - 1)
        end
        return FiniteDual(x.l^y, x.r^y, vd * x.Δ)
    else
        return exp(y * log(x))
    end
end
Base.:^(x, y::FiniteDual) = exp(y * log(x))
    
function Base.exp(x::FiniteDual)
    vl = exp(x.l)
    vr = exp(x.r)
    if vl == vr
        return FiniteDual(vl, vr, zero(vl))
    else
        return FiniteDual(vl, vr, vr * expm1(x.l - x.r) / (x.l - x.r) * x.Δ)
    end
end

Base.expm1(x::FiniteDual) = FiniteDual(expm1(x.l), expm1(x.r), (exp(x) - 1.0).Δ)

function Base.log(x::FiniteDual)
    vl = log(x.l)
    vr = log(x.r)
    if vl == vr
        return FiniteDual(vl, vr, zero(vl))
    else
        return FiniteDual(vl, vr, log1p((x.l - x.r) / x.r) / (x.l - x.r) * x.Δ)
    end
end

Base.log1p(x::FiniteDual) = FiniteDual(log1p(x.l), log1p(x.r), log(x+1.0).Δ)

# sin/cos #
#--------#

function Base.sin(x::FiniteDual)
    vl = sin(x.l)
    vr = sin(x.r)
    if vl == vr
        return FiniteDual(vl, vr, zero(vl))
    else
        return FiniteDual(vl, vr, cos((x.l + x.r) / 2) * sinc((x.l - x.r) / 2pi) * x.Δ)
    end
end

function Base.cos(x::FiniteDual)
    vl = cos(x.l)
    vr = cos(x.r)
    if vl == vr
        return FiniteDual(vl, vr, zero(vl))
    else
        return FiniteDual(vl, vr, -sin((x.l + x.r) / 2) * sinc((x.l - x.r) / 2pi) * x.Δ)
    end
end

Base.sinc(x::FiniteDual) = FiniteDual(sinc(x.l), sinc(x.r), (sin(pi*x)/(pi*x)).Δ)


###################
# Pretty Printing #
###################

function Base.show(io::IO, x::FiniteDual)
    print(io, "FiniteDual(", x.l)
    for xx in [x.r, x.Δ]
        print(io, ",", xx)
    end
    print(io, ")")
end