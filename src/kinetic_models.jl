
export InterpK, get_ks
export MYEGA_k, Sigmoid_k

# ------------------------------
# Fit values directly by interpolating 
"""
    $(SIGNATURES)

InterpK interpolates in log space, then takes exp and applies a scale factor.
At construction, supply the a of log(k) values for (k1, k2/k1, k3/k2, k4/k3) to enforce monotonicity.
"""
struct InterpK{itype<:DataInterpolations.AbstractInterpolation, rtype<:Number} 
    spl::itype
    scale::rtype
end
(ik::InterpK)(T) = exp(ik.spl(T))*ik.scale
function InterpK(klog, T, scale)
    k_increasing = cumsum(klog)
    return InterpK(LinearInterpolation(k_increasing, T, extrapolation=ExtrapolationType.Linear), scale)
end
function InterpK(klog, T)
    return InterpK(klog, T, 1e-14u"s^-1")
end
function InterpK(klog)
    T = [-45.0, -38, -32, -28]u"°C" .|> u"K"
    return InterpK(klog, T)
    # return InterpK(LinearInterpolation(klog, T_stops,extrapolation=ExtrapolationType.Linear))
end
InterpK(;klog, T) = InterpK(klog, T)
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
MYEGA_k(;k0, K1, K2) = MYEGA_k(k0, K1, K2)
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
Sigmoid_k(;k0, kf, K1, K2) = Sigmoid_k(k0, kf, K1, K2)
Sigmoid_k(a) = Sigmoid_k(a...)
function (pk::Sigmoid_k)(T)
    e_x = exp((T-pk.K2)/pk.K1)
    return (pk.k0 + pk.kf*e_x)/(1 + e_x)
end