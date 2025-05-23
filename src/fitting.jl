export identify_pd_end_2der, identify_pd_end_onoff

"""
    $(SIGNATURES)
    
Identify the end of primary drying using second derivative of Pirani pressure.

# Arguments
- `data`: Table containing the Pirani pressure `pch_pir` and time `t`
- `window_width`: Width of the Savitzky-Golay filter window (default: 91)
- `tmin`: Minimum time for analysis (default: 0 hours, unitful)
- `tmax`: Maximum time for analysis (default: Inf hours, unitful)

Very simple implementation:

    pch_pir_der2 = savitzky_golay(data.pch_pir, window_width, 3, deriv=2).y
    t_end = data.t[argmax(pch_pir_der2[tmin .< data.t .< tmax])] + tmin
    return t_end
"""
function identify_pd_end_2der(data, window_width=91, tmin=0u"hr", tmax=Inf*u"hr")
    pch_pir_der2 = savitzky_golay(data.pch_pir, window_width, 3, deriv=2).y
    t_end = data.t[argmax(pch_pir_der2[tmin .< data.t .< tmax])] + tmin
    return t_end
end
"""
    $(SIGNATURES)

Identify a range of times for the end of primary drying, with the onset and offset of the
Pirani pressure curve.

Use a Savitzky-Golay filter to smooth, identify minimum of derivative, and then find linear
intersects with maximum and minimum values of pressure. Those intersects are taken as onset 
and offset.
"""
function identify_pd_end_onoff(data, window_width=91, tmin=0u"hr", tmax=Inf*u"hr")
    pch_pir_sm = savitzky_golay(data.pch_pir, window_width, 3, deriv=0).y
    dt = sum(diff(data.t)) / (length(data.t)-1)
    pch_pir_der1 = savitzky_golay(data.pch_pir, window_width, 3, deriv=1, rate=1/dt).y
    ti_inds = tmin .< data.t .< tmax
    i_mid = argmin(pch_pir_der1[ti_inds]) + searchsortedfirst(data.t, tmin)
    t_mid = data.t[i_mid]
    dp_mid = pch_pir_der1[i_mid]
    p_mid = pch_pir_sm[i_mid]
    p_max = maximum(pch_pir_sm[ti_inds])
    p_min = minimum(pch_pir_sm[ti_inds])
    t_onset = (p_max - p_mid)/dp_mid + t_mid
    t_offset = (p_min - p_mid)/dp_mid + t_mid
    return (t_onset, t_offset)
end

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

export k_f_nt

k_f_nt(nt::NamedTuple{(:k0, :kf, :K1, :K2)}) = Sigmoid_k(nt.k0, nt.kf, nt.K1, nt.K2)
k_f_nt(nt::NamedTuple{(:k0, :K1, :K2)}) = MYEGA_k(nt.k0, nt.K1, nt.K2)
k_f_nt(nt) = error("Mis-ordered k_f parameters. This should be handled more carefully")
