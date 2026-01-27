
# ---------------------------------
function sigmoid_transform_basic(k0, kf, K1, K2)
    tr = as((
        k_f = as((
            k0 = TVScale(k0) ∘ TVExp(),
            kf = TVScale(kf) ∘ TVExp(),
            K1 = TVScale(K1) ∘ TVExp(),
            K2 = TVScale(K2) ∘ TVExp(),
        )),
    ))
end
function sigmoid_transform_bounded(kf, K1; extrafac=1)
    Nmin = 1e0u"μm^-3"
    Vmin=0.5u"cm^3"
    tmin=1.0u"hr"
    k0_max = 1/Nmin/Vmin/tmin *extrafac|> u"s^-1" 
    tr = as((
        k_f = as((
            k0 = TVScale(k0_max) ∘ TVLogistic(),
            kf = TVScale(kf) ∘ TVExp(),
            K1 = TVScale(K1) ∘ TVExp(),
            K2 = TVScale(1.0u"K") ∘ TVShift(216.0) ∘ TVScale(45.0) ∘ TVLogistic() ∘TVScale(0.1),
        )),
    ))
end
function myega_transform_basic(k0, K1, K2)
    tr = as((
        k_f = as((
            k0 = TVScale(k0) ∘ TVExp(),
            K1 = TVScale(K1) ∘ TVExp(),
            K2 = TVScale(K2) ∘ TVExp(),
        )),
    ))
end
function porous_transform_bounded(Vfill)
    tr = as((
        N = ConstWrapTV() ∘ TVScale(1.0u"μm^-3") ∘ TVLogistic() ∘ TVShift(TransformVariables.logit(1e-3)),
        S = ConstWrapTV() ∘ TVScale(Vfill) ∘ TVLogistic() ∘ TVShift(TransformVariables.logit(1e-7)),
        C = TVScale(3π/sqrt(8) * 10 * sqrt(1.0u"mm")) ∘ TVLogistic() ∘ TVShift(TransformVariables.logit(1e-1)),
        Rs = TVScale(1.0u"cm^2*Torr*hr/g") ∘ TVLogistic(),
        ))
end
function porous_transform_basic(Ng, Sg, Cg, Rsg)
    tr = as((
        N = ConstWrapTV() ∘ TVScale(Ng) ∘ TVExp(),
        S = ConstWrapTV() ∘ TVScale(Sg) ∘ TVExp(),
        C = TVScale(Cg) ∘ TVExp(),
        Rs = TVScale(Rsg) ∘ TVExp(),
    ))
end

function badNS(nt) 
    if hasproperty(nt, :N) && hasproperty(nt, :S)
        return (nt.N(0) < 1/nt.S(0))
    end
    error("Unhandled nt in badNS")
end

"""
    $(SIGNATURES)

Construct a typical transform for fitting to the microcollapse model with MYEGA fitting.
"""
function mc_myega_transform_basic(Ng, Sg, Cg, Rsg, k0g, K1g, K2g)
    k = myega_transform_basic(k0g, K1g, K2g).transformations.k
    p = porous_transform_basic(Ng, Sg, Cg, Rsg).transformations
    tr = as(merge(p, (;k)))
    return tr
end
"""
    $(SIGNATURES)

Construct a typical transform for fitting to the microcollapse model with MYEGA fitting.
"""
function mc_sigmoid_transform_basic(Ng, Sg, Cg, Rsg, k0g, kfg, K1g, K2g)
    k = sigmoid_transform_basic(k0g, kfg, K1g, K2g).transforms.k
    p = porous_transform_basic(Ng, Sg, Cg, Rsg).transforms
    tr = as(merge(p, (;k)))
    return tr
end

# --------------------------------- 
# ----------------------------------------
# Fitting with microcollapse model for Rp


k_f_nt(nt::NamedTuple{(:k0, :kf, :K1, :K2)}) = Sigmoid_k(nt.k0, nt.kf, nt.K1, nt.K2)
k_f_nt(nt::NamedTuple{(:k0, :K1, :K2)}) = MYEGA_k(nt.k0, nt.K1, nt.K2)
k_f_nt(nt) = error("Mis-ordered k_f parameters. This should be handled more carefully")
