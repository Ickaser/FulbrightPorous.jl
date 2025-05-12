export gennsol_Rp, objn_Rp
export gen_sol_mc_arb, gennsol_mc, objn_mc

# ---------------------------------
export identify_pd_end_2der

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
gennsol_Rp(fitprm, tr, pos) = [gen_sol_pd(fitprm, tr, po) for po in pos]
gennsol_Rp(fitprm, tr, pos, fitdats) = [gen_sol_pd(fitprm, tr, po, fitdat) for (po,fitdat) in zip(pos,fitdats)]
function objn_Rp(fitprm, tpfs; tweight=1e-3)
    tr, pos, fitdats = tpfs
    sols = gennsol_Rp(fitprm, tr, pos, fitdats)
    objs = map(sols, fitdats) do sol, fitdat
        obj_expT(sol, fitdat, tweight=tweight)
    end
    return sum(objs)
end
# ----------------------------------------
# Fitting with microcollapse model for Rp

# Estimates of N and S:
# d_typical = 50u"μm"
# V_typical = 4/3*π*d_typical^3
# N_typical = 1/V_typical
# S_typical = V_typical^2*N_typical = V_typical
export k_f_nt

k_f_nt(nt::NamedTuple{(:k0, :kf, :K1, :K2)}) = Sigmoid_k(nt.k0, nt.kf, nt.K1, nt.K2)
k_f_nt(nt::NamedTuple{(:k0, :K1, :K2)}) = MYEGA_k(nt.k0, nt.K1, nt.K2)
k_f_nt(nt) = error("Mis-ordered k_f parameters. This should be handled more carefully")

function gen_sol_mc_arb(fitprmlog, tr, po; verbose=false, saveat=[])
    fitprm = transform(tr, fitprmlog)
    verbose && @info "fitprm" fitprm
    new_params = setproperties(po, fitprm)
    ϕ = 1 - new_params.csolid/new_params.ρsolution
    if new_params.S(0) < ϕ/new_params.N(0)
        @warn "S is too small: implies negative variance"
        # @reset new_params.S = ConstPhysProp(ϕ/new_params.N(0))
        return NaN
    end
    prob = ODEProblem(new_params; tspan=(0.0, 1000.0))
    sol = solve(prob, odealg_chunk4; saveat=saveat)
    return sol
end
function gen_sol_mc_arb(fitprmlog, tr, po, fitdat; verbose=false)
    return gen_sol_mc_arb(fitprmlog, tr, po; verbose=verbose, saveat=ustrip.(u"hr", fitdat.t))
end
# function gen_sol_mc_arb(fitprmlog, po, k_f_gen; verbose=false, fit_scl=fit_scl_porous, saveat=[])
#     new_params = deepcopy(po)
#     fitprm = exp.(fitprmlog) # Constrain all to be positive by taking exp
#     verbose && @info "fitprm" fitprm
#     @reset new_params.N = ConstPhysProp(fitprm[1]*fit_scl[1])
#     @reset new_params.S = ConstPhysProp(fitprm[2]*fit_scl[2])
#     @reset new_params.C  = fitprm[3]*fit_scl[3]
#     @reset new_params.Rs = fitprm[4]*fit_scl[4]
#     @reset new_params.k_f = k_f_gen(fitprmlog[5:end]) # Exp taken elsewhere
#     ϕ = 1 - new_params.csolid/new_params.ρsolution
#     if new_params.S(0) < ϕ/new_params.N(0)
#         @warn "S is too small: implies negative variance"
#         # @reset new_params.S = ConstPhysProp(ϕ/new_params.N(0))
#         return NaN
#     end
#     prob = ODEProblem(new_params; tspan=(0.0, 1000.0))
#     sol = solve(prob, Rodas3())
#     return sol
# end
# function gen_sol_mc_arb(fitprmlog, po, fitdat, k_f_gen; 
#                         verbose=false, fit_scl = fit_scl_porous)
#     # new_params = deepcopy(po)
#     # fitprm = exp.(fitprmlog) # Constrain all to be positive by taking exp
#     # verbose && @info "fitprm" fitprm
#     # @reset new_params.N = ConstPhysProp(fitprm[1]*fit_scl[1])
#     # @reset new_params.S = ConstPhysProp(fitprm[2]*fit_scl[2])
#     # @reset new_params.C  = fitprm[3]*fit_scl[3]
#     # @reset new_params.Rs = fitprm[4]*fit_scl[4]
#     # @reset new_params.k_f = k_f_gen(fitprmlog[5:end]) # Exp taken elsewhere
#     # if new_params.S(0) < ϕ/new_params.N(0)
#     #     @warn "S is too small: implies negative variance" new_params.S(0) ϕ/new_params.N(0)
#     #     new_params.S = ConstPhysProp(ϕ/new_params.N(0))
#     # end
#     # prob = ODEProblem(new_params; tspan=(0.0, 1000.0))
#     # sol = solve(prob, Rodas3(); saveat=ustrip.(u"hr", fitdat.t))
#     sol = gen_sol_mc_arb(fitprmlog, po, k_f_gen; verbose=verbose, fit_scl=fit_scl, saveat=ustrip.(u"hr", fitdat.t))
#     return sol
# end
function gennsol_mc(fitprm, tr, pos; kwargs...)
    sols = map(pos) do po
        gen_sol_mc_arb(fitprm, tr, po; kwargs...)
    end
    return sols
end
function gennsol_mc(fitprm, tr, pos, fitdats; kwargs...)
    sols = map(pos, fitdats) do po, fitdat
        gen_sol_mc_arb(fitprm, tr, po, fitdat; kwargs...)
    end
    return sols
end

function objn_mc(fitprmlog, tpf; tweight=0, verbose=false)
    rtype = eltype(fitprmlog)
    tr, pos, fitdats = tpf
    # verbose && @info "fitprm" fitprm
    if pos isa ParamObjMC && fitdats isa PrimaryDryFit # single element
        fitprm = transform(tr, fitprmlog)
        if fitprm.S(0) < (1 - fitprm.csolid/fitprm.ρsolution)/fitprm.N(0)
            @warn "S is too small: implies negative variance" fitprm.S(0) 1/fitprm.N(0)
            return rtype(NaN)
        end
        sols = gen_sol_pd(fitprmlog, tr, pos, fitdats)
        objs = [obj_expT(sols, fitdats; tweight)]
    else
        sols = gen_nsol_pd(fitprmlog, tr, pos, fitdats)
        for sol in sols
            fitprm = sol.prob.p
            if fitprm.S(0) < (1 - fitprm.csolid/fitprm.ρsolution)/fitprm.N(0)
                @warn "S is too small: implies negative variance" fitprm.S(0) 1/fitprm.N(0)
                return rtype(NaN)
            end
            if sol.retcode != ReturnCode.Terminated
                hf = length(sol)>0 ? sol[begin,end] : sol.prob.p.hf0
                @warn "unfinished drying" hf fitprmlog 
            end
        end
        objs = map(sols, fitdats) do sol, fitdat
            obj_expT(sol, fitdat; tweight)
        end
    end
    if verbose
        if any([sol.retcode != ReturnCode.Terminated for sol in sols])
            for sol in sols 
                if sol.retcode != ReturnCode.Terminated
                    @info "unfinished drying" sol[1,end] fitprmlog 
                end
            end
        end
        # @info "unfinished drying" [sol for 
    end
    verbose && @info "obj" objs fitprmlog 
    return rtype(sum(objs))
end