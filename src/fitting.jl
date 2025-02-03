export gen_sol_Rp, gennsol_Rp, obj_Rp, objn_Rp
export gen_sol_mc_arb, gennsol_mc, objn_mc


# --------------------------------- 
# Conventional fitting for Rp
const Rp_scale = (0.5u"cm^2*Torr*hr/g", 20.0u"cm*Torr*hr/g", 1.0u"1/cm")
function gen_sol_Rp(a123, po)
    new_params = @set po.Rp = RpFormFit((exp.(a123).*Rp_scale)...)
    prob = ODEProblem(new_params; tspan=(0.0, 1000.0))
    sol = solve(prob, Rodas3())
    return sol
end
function gen_sol_Rp(a123, po, fitdat)
    new_params = @set po.Rp = RpFormFit((exp.(a123).*Rp_scale)...)
    prob = ODEProblem(new_params; tspan=(0.0, 1000.0))
    sol = solve(prob, Rodas3(); saveat=ustrip.(u"hr", fitdat.t), )
    return sol
end
function obj_Rp(x, p; tweight=1e-4, verbose=false)
    sol = gen_sol_Rp(x, p[1], p[2])
    return obj_expT(sol, p[2], tweight=tweight, verbose=verbose)
end
gennsol_Rp(fitprm, pos) = [gen_sol_Rp(fitprm, po) for po in pos]
gennsol_Rp(fitprm, pos, fitdats) = [gen_sol_Rp(fitprm, po, fitdat) for (po,fitdat) in zip(pos,fitdats)]
function objn_Rp(fitprm, pfs; tweight=1e-3)
    pos, fitdats = pfs
    sols = gennsol_Rp(fitprm, pos, fitdats)
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

const fit_scl_porous = (2e6u"cm^-3", 5e-7u"cm^3", 2e-1u"m^0.5", 0.1u"cm^2*Torr*hr/g")
function gen_sol_mc_arb(fitprm, po, k_f_gen; verbose=false)
    new_params = deepcopy(po)
    # fitprm = exp.(fitprmlog) # Constrain all to be positive by taking exp
    verbose && @info "fitprm" fitprm
    @reset new_params.N = ConstPhysProp(fitprm[1]*fit_scl_porous[1])
    @reset new_params.S = ConstPhysProp(fitprm[2]*fit_scl_porous[2])
    @reset new_params.C  = fitprm[3]*fit_scl_porous[3]
    @reset new_params.Rs = fitprm[4]*fit_scl_porous[4]
    @reset new_params.k_f = k_f_gen(fitprm[5:end])
    prob = ODEProblem(new_params; tspan=(0.0, 1000.0))
    sol = solve(prob, Rodas3())
    return sol
end
function gen_sol_mc_arb(fitprm, po, fitdat, k_f_gen; verbose=false)
    new_params = deepcopy(po)
    # fitprm = exp.(fitprmlog) # Constrain all to be positive by taking exp
    verbose && @info "fitprm" fitprm
    @reset new_params.N = ConstPhysProp(fitprm[1]*fit_scl_porous[1])
    @reset new_params.S = ConstPhysProp(fitprm[2]*fit_scl_porous[2])
    @reset new_params.C  = fitprm[3]*fit_scl_porous[3]
    @reset new_params.Rs = fitprm[4]*fit_scl_porous[4]
    @reset new_params.k_f = k_f_gen(fitprm[5:end])
    prob = ODEProblem(new_params; tspan=(0.0, 1000.0))
    sol = solve(prob, Rodas3(); saveat=ustrip.(u"hr", fitdat.t))
    return sol
end
function gennsol_mc(fitprm, pos, k_f_gen)
    sols = map(pos) do po
        gen_sol_mc_arb(fitprm, po, k_f_gen)
    end
    return sols
end
function gennsol_mc(fitprm, pos, fitdats, k_f_gen)
    sols = map(pos, fitdats) do po, fitdat
        gen_sol_mc_arb(fitprm, po, fitdat, k_f_gen)
    end
    return sols
end

function objn_mc(fitprm, pfk; tweight=0, verbose=false)
    pos, fitdats, k_gen = pfk
    # verbose && @info "fitprm" fitprm
    if pos isa ParamObjMC && fitdats isa PrimaryDryFit # single element
        sols = gen_sol_mc_arb(fitprm, pos, fitdats, k_gen)
        objs = [obj_expT(sols, fitdats, tweight=tweight)]
    else
        sols = gennsol_mc(fitprm, pos, fitdats, k_gen)
        objs = map(sols, fitdats) do sol, fitdat
            obj_expT(sol, fitdat, tweight=tweight)
        end
    end
    if verbose
        if any([sol.retcode != ReturnCode.Terminated for sol in sols])
            for sol in sols 
                if sol.retcode != ReturnCode.Terminated
                    @info "unfinished drying" sol[1,end] fitprm 
                end
            end
        end
        # @info "unfinished drying" [sol for 
    end
    verbose && @info "obj" objs fitprm 
    return sum(objs)
end