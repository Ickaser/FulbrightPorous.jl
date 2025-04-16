
export InterpK, get_ks
export MYEGA_k, Sigmoid_k

# ------------------------------
# Fit values directly by interpolating 
struct InterpK{itype<:DataInterpolations.AbstractInterpolation, rtype<:Number} 
    spl::itype
    scale::rtype
end
(ik::InterpK)(T) = exp(ik.spl(T))*ik.scale
function InterpK(klog::AbstractVector, T, scale)
    return InterpK(LinearInterpolation(klog, T, extrapolation=ExtrapolationType.Linear), scale)
end
function InterpK(klog::AbstractVector, T)
    return InterpK(klog, T, 1e-14u"s^-1")
end
function InterpK(klog::AbstractVector)
    T = [-45.0, -38, -32, -28]u"°C" .|> u"K"
    return InterpK(klog, T)
    # return InterpK(LinearInterpolation(klog, T_stops,extrapolation=ExtrapolationType.Linear))
end
function Base.show(io::IO, ik::InterpK)
    str = "InterpK($(ik.spl.u), $(ik.spl.t), $(ik.scale))"
    return print(io, str)
end

function get_ks(ik::InterpK)
    return ik.(ik.spl.t)
end

# ------------------------
# Fit to MYEGA model of viscosity
struct MYEGA_k{T1<:Number, T2<:Number}
    k0::T1
    K1::T2
    K2::T2
end
MYEGA_k(a) = MYEGA_k(a...)
(mk::MYEGA_k)(T) = mk.k0*exp((-mk.K1/T)*exp(mk.K2/T))

# ------------------------
# Fit to a sigmoid form
struct Sigmoid_k{T1<:Number, T2<:Number}
    k0::T1
    kf::T1
    K1::T2
    K2::T2
end
Sigmoid_k(a) = Sigmoid_k(a...)
function (pk::Sigmoid_k)(T)
    rat = pk.k0/pk.kf
    x = (T-pk.K2)/pk.K1
    # return k0/(tanh(-x)+1) # This has the right behavior...
    return pk.k0*(exp(x) + exp(-x))/((1+rat)*exp(-x)+rat*exp(x)) # but this way it can avoid numerical issues
end