
# ---------------------------------
function sigmoid_transform_basic(k0, kf, K1, K2)
tr = as((
    k = as((
        k0 = TVScale(k0) ∘ TVExp(),
        kf = TVScale(kf) ∘ TVExp(),
        K1 = TVScale(K1) ∘ TVExp(),
        K2 = TVScale(K2) ∘ TVExp(),
        )),
    ))
end
function myega_transform_basic(k0, K1, K2)
tr = as((
    k = as((
        k0 = TVScale(k0) ∘ TVExp(),
        K1 = TVScale(K1) ∘ TVExp(),
        K2 = TVScale(K2) ∘ TVExp(),
        )),
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
