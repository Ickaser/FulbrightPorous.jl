using DrWatson
@quickactivate :FulbrightPores
using ForwardDiff, SciMLSensitivity
using DataInterpolations
using UnitfulLatexify, Latexify
using PrettyTables

const NLSP = NonlinearLeastSquaresProblem
const LM = LevenbergMarquardt
const LMa = LevenbergMarquardt(autodiff=AutoFiniteDiff())
const LMf = LevenbergMarquardt(autodiff=AutoForwardDiff())
pythonplot()

tcond = (termination_condition=RelTerminationMode(),)

# ------------------
# Kv fitting to FIg 7

pKvdat = [ 50.0 4.18; 75.0 5.25;
100. 6.21; 125. 7.08; 150. 7.81;
175. 8.53; 200. 9.09; 250. 10.0; ] # WebPlotDigitized from Fig 7 of IntJPharm
pd = pKvdat[:,1]*u"mTorr"
Kvd = pKvdat[:,2]*u"J/hr/cm^2/K"
function fitKv(KKK) 
    Kc, Kp, Kd = KKK .* [u"W/m^2/K", u"W/m^2/K/Torr", u"1/Torr"]
    Kv = RpFormFit(Kc, Kp, Kd)
    err = @. Kv(pd) - Kvd
    return ustrip(sum(abs2.(err))), Kv
end
optKv = solve(NLSP((x,p,t)->fitKv(x)[1], [10.0, 1, 1]), LMf)
K_shf = fitKv(optKv.u)[2]

# --------------------------
# Cycle 1
# Labeled as case 4 in Table 8

dat1raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle1.csv"), Table)
dat1 = map(simtra_rename, dat1raw)
# dat1 = dat1[dat1.T2 .< 250u"K"]

fitdat1 = PrimaryDryFit(collect(dat1.t), (
    collect(filter(x->(!isnan(x) && x<-32.5u"°C"), dat1.T1)),
    collect(filter(!isnan, dat1.T2[dat1.t .< 60u"hr"])),
    # collect(filter(!isnan, dat1.T3)),
))

plot(fitdat1)
pch1 = RampedVariable(100u"mTorr")
@df dat1 plot(:t, :pch_cm, ylim=(0, 200))
plot!(pch1)
# There's a slight delay before shelf temperature ramps
Tsh1 = RampedVariable(uconvert.(u"K", [-40, -40, -25]u"°C"), [1u"K/minute", 0.5u"K/minute"], [23u"minute"])
@df dat1 plot(:t, :Tsh, xlim=(0, 2))
plot!(Tsh1)

# ------------------------------
# Cycle 2

dat2raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle2.csv"), Table)
dat2 = map(simtra_rename, dat2raw)
filter!(x->x.t .> 0u"hr", dat2)

begin
@df dat1 exptfplot(:t, :T1, :T2, :T3)
@df dat2 exptfplot!(:t, :T1, :T2, :T3, c=:red)
plot!(ylim=(-40, -30))
end

fitdat2 = @df dat2 PrimaryDryFit(collect(:t), (
    # collect(filter(!isnan, :T1)),
    collect(filter(!isnan, :T2[:t .< 85u"hr"])),
    collect(filter(!isnan, :T3)),))
plot(fitdat2)

pch2 = RampedVariable(100u"mTorr")
@df dat2 plot(:t, :pch_cm, ylim=(0, 200))
plot!(pch2)
Tsh2 = Tsh1 - 5u"K"
@df dat2 plot(:t, :Tsh, xlim=(0, 2))
plot!(Tsh2)

# -------------------------------  
# Cycle 3
dat3raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle3.csv"), Table)
dat3 = map(simtra_rename, dat3raw)
@df dat3 exptfplot(:t, :T1, :T2, :T3)

pch3 = RampedVariable(60u"mTorr")
@df dat3 plot(:t, :pch_cm, ylim=(0, 200))
plot!(pch3)
Tsh3 = RampedVariable(uconvert.(u"K", [-40, -25]u"°C"), 0.5u"K/minute")
@df dat3 plot(:t, :Tsh, xlim=(0, 2))
plot!(Tsh3)

# -----------------
# Cycle 5
dat5raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle5.csv"), Table)
dat5 = map(simtra_rename, dat5raw)
@df dat5 exptfplot(:t, :T1, :T2, :T3)

pch5 = RampedVariable(50u"mTorr")
@df dat5 plot(:t, :pch_cm, ylim=(0, 200))
plot!(pch5)
Tsh5 = RampedVariable(uconvert.(u"K", [-40, -40, -20]u"°C"), [1u"K/minute", 0.5u"K/minute"], [12u"minute"])
@df dat5 plot(:t, :Tsh_i, xlim=(0, 2))
plot!(Tsh5)


# -----------------
# Cycle 8, high temperature
dat8raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle8.csv"), Table)
dat8 = map(simtra_rename, dat8raw)
@df dat8 exptfplot(:t, :T1, :T2, :T3)
@df dat8 plot!(:t, :Tsh_i, xlim=(0, 20), ylim=(-40, -20))

pch8 = RampedVariable(100u"mTorr")
@df dat8 plot(:t, :pch_cm, ylim=(0, 200))
plot!(pch8)
Tsh8 = RampedVariable(uconvert.(u"K", [-40, -40, -10]u"°C"), [1u"K/minute", 0.5u"K/minute"], [4u"minute"])
@df dat8 plot(:t, :Tsh_i, xlim=(0, 2))
plot!(Tsh8)

fitdat8 = @df dat8 PrimaryDryFit(collect(:t), (
    collect(filter(!isnan, :T1)),
    collect(filter(!isnan, :T2[:t .< 17u"hr"])),
    collect(filter(!isnan, :T3[:t .< 17u"hr"])),))
plot(fitdat8)


# -------------
# Model setup

rad_i = 1.1u"cm"
rad_o = 1.2u"cm"
A_p = π*rad_i^2
A_v = π*rad_o^2
fillvol = 5u"mL"
hf0 = fillvol/A_p

Rp_g = RpFormFit(0.8u"cm^2*Torr*hr/g", 14u"cm*Torr*hr/g", 1u"1/cm")

c_solid = 0.04u"g/mL"
ρ_solution = 1u"g/mL"

pb_pikal = Vector{Any}([
    (Rp_g, hf0, c_solid, ρ_solution),
    (K_shf, A_v, A_p),
    (pch1, Tsh1),
])

pop1 = ParamObjPikal(pb_pikal)
pop2 = deepcopy(pop1)
@reset pop2.pch = pch2
@reset pop2.Tsh = Tsh2
pop8 = deepcopy(pop1)
@reset pop8.pch = pch8
@reset pop8.Tsh = Tsh8

function gen_sol_Rp(a123, po)
    Rp_un = [u"cm^2*Torr*hr/g", u"cm*Torr*hr/g", u"1/cm"]
    new_params = @set po.Rp = RpFormFit((a123.*Rp_un)...)
    prob = ODEProblem(new_params)
    sol = solve(prob, Rodas3(autodiff=false))
    return sol, new_params
end
obj_Rp = (x,p)->[obj_expT(gen_sol_Rp(x, p[1])[1], p[2])]
@time obj_Rp([0.8, 14, 1], (pop1, fitdat1))

opt_Rp1 = solve(NLSP(obj_Rp, [0.8, 14, 1], p=(pop1, fitdat1)), LMf; reltol=1e-1, show_trace=Val(true), tcond...)

sol1_Rp1, Rp1_params = gen_sol_Rp(opt_Rp1.u, pop1);
plot(u"hr", u"°C", xlim=(0, 60), ylim=(-40, -30))
modconvtplot!(sol1_Rp1)
@df dat1 exptfplot!(:t, :T1, :T2, :T3)

# -------------- Provide guess values for new parameters

# typical diameter: 50 μm
d_typical = 50u"μm"
# typical total volume: 
V_typical = 4/3*π*d_typical^3
N_typical = fillvol/V_typical

# Good guesses for Knudsen
N = RampedVariable(N_typical / hf0)
S = RampedVariable(V_typical^2*N_typical/hf0)
C = 2e-1u"m^0.5"
R0 = 0.05u"cm^2*Torr*hr/g"

# Good guesses for viscous
# N = RampedVariable(N_typical / hf0) 
# S = RampedVariable(1e0u"cm^5")
# C = 2e2u"m"
# R0 = 0.01u"cm^2*Torr*hr/g"

k0 = 1e-16u"1/s"
E = u"R"*Tsh1(0u"s")
kf = T->k0*exp(-E/(R*T))

pb_mc = (
    (Rp_g, hf0, c_solid, ρ_solution),
    (K_shf, A_v, A_p),
    (pch1, Tsh1),
    (kf, C, N, S, R0)
)
pomc1 = ParamObjMC(pb_mc)
pomc2 = deepcopy(pomc1)
@reset pomc2.pch = pch2
@reset pomc2.Tsh = Tsh2
pomc8 = deepcopy(pomc1)
@reset pomc8.pch = pch8
@reset pomc8.Tsh = Tsh8

# tprob = ODEProblem(pomc1)
# tsol = solve(tprob, Rodas3(autodiff=true))
# poreevolplot(tsol)

fit_un_arr = [u"cm^-1", u"cm^5", u"m^0.5", u"m/s", u"s^-1", u"J/mol"]
fitguess = ustrip.(fit_un_arr, [N(0), S(0), C, R0, k0, E, ])
fitguesslog = log.(fitguess)
function gen_sol_constNS_arr(fitprmlog, po; verbose=false)
    new_params = deepcopy(po)
    fitprm = exp.(fitprmlog) # Constrain all to be positive by taking exp
    verbose && @info "fitprm" fitprm
    @reset new_params.N = RampedVariable(fitprm[1]*fit_un_arr[1])
    @reset new_params.S = RampedVariable(fitprm[2]*fit_un_arr[2])
    @reset new_params.C  = fitprm[3]*fit_un_arr[3]
    @reset new_params.R0 = fitprm[4]*fit_un_arr[4]
    # Arrhenius fit
    k0 = fitprm[6]*fit_un_arr[5]
    E = fitprm[6]*fit_un_arr[6]
    @reset new_params.k_f = T->k0*exp(-E/(R*T))
    prob = ODEProblem(new_params)
    sol = solve(prob, Rodas3(autodiff=true))
    return sol, new_params
end

fit_un_arb = [u"cm^-1", u"cm^5", u"m^0.5", u"m/s", fill(u"s^-1", 4)...]
fitguess = ustrip.(fit_un_arb, [N(0), S(0), C, R0, k0, k0*1e1, k0*1e2, k0*1e4])
fitguesslog = log.(fitguess)
Tstops = [-45.0, -38, -32, -28]u"°C" .|> u"K"
function gen_sol_constNS_arb(fitprmlog, po; verbose=false)
    new_params = deepcopy(po)
    fitprm = exp.(fitprmlog) # Constrain all to be positive by taking exp
    verbose && @info "fitprm" fitprm
    @reset new_params.N = RampedVariable(fitprm[1]*fit_un_arb[1])
    @reset new_params.S = RampedVariable(fitprm[2]*fit_un_arb[2])
    @reset new_params.C  = fitprm[3]*fit_un_arb[3]
    @reset new_params.R0 = fitprm[4]*fit_un_arb[4]
    # Linear interpolant

    # kstops = fitprm[5:8]*fit_un_arb[3]
    # spl = LinearInterpolation(kstops, Tstops, extrapolate=true)
    # Lin interp in log space
    klogs = fitprmlog[5:8]
    # klogs = log.(fitprm[5:8])
    spll = LinearInterpolation(klogs, Tstops, extrapolate=true)
    spl = T->exp(spll(T))*fit_un_arb[5]
    # PCHIP interp
    # spl = PCHIPInterpolation(fitprm[5:8]*fit_un_arb[3], Tstops, extrapolate = true)
    @reset new_params.k_f = spl
    prob = ODEProblem(new_params)
    # u0 = ustrip.([u"cm", u"K", u"cm^6", NoUnits], [new_params.hf0, new_params.Tsh(0u"s"), 0.0u"cm^6", 0.0])
    # prob = ODEProblem(lyo_mc_dae_f2, u0, (0.0,200.0), new_params, callback=end_drying_callback)
    sol = solve(prob, Rodas3(autodiff=true))
    return sol, new_params
end

# tsol, tparams = gen_sol_constNS_arb(fitguesslog, pomc1); 
# modconvtplot(tsol)
# plot!(fitdat1)

obj_MC = (x,p)->[(issorted(x[5:8]) ? (obj_expT(gen_sol_constNS_arb(x, p[1])[1], p[2])) : NaN)]
@time obj_MC(fitguesslog, (pomc1, fitdat1))

@time opt_MC1 = solve(NLSP(obj_MC, fitguesslog, p=(pomc1, fitdat1)), LMf; show_trace=Val(true), maxiters=500)
obj_MC(opt_MC1.u, (pomc1, fitdat1))

exp.(opt_MC1.u) .*fit_un_arb

sol1_MC1, MC1_params = gen_sol_constNS_arb(opt_MC1.u, pomc1);


# ------ Fit to Cycle 2 and Cycle 8 data separately

opt2p = solve(NLSP(obj_Rp, [0.8, 14, 1], p=(pop2, fitdat2)), LM(); tcond..., reltol=1e-1)
opt_MC2 = solve(NLSP(obj_MC, opt_MC1.u, p=(pomc2, fitdat2)), LM(); tcond..., reltol=1e-1, maxiters=300, show_trace=Val(true), maxtime=20)

# obj8p = x->obj_expT(gen_sol_Rp(x, pop8)[1], fitdat8)
# opt8p = optimize(obj8p, opt_Rp1.minimizer, LevenbergMarquardt(), )
opt8p = solve(NLSP(obj_Rp, [0.8, 14, 1], p=(pop8,fitdat8)), LM(); tcond..., reltol=1e-1)
# Enforce monontonic k in objective
# obj8mc = x->(issorted(x[5:8]) ? obj_expT(gen_sol_constNS_arb(x, pomc8)[1], fitdat8) : NaN)
opt_MC8 = solve(NLSP(obj_MC, opt_MC2.u, p=(pomc8, fitdat8)), LMa; maxiters=1000)

sol2_Rp2, Rp2_params = gen_sol_Rp(opt2p.u, pop2);
sol8_Rp8, Rp8_params = gen_sol_Rp(opt8p.u, pop8);
sol2_MC2, MC2_params = gen_sol_constNS_arb(opt_MC2.u, pomc2);
sol8_MC8, MC8_params = gen_sol_constNS_arb(opt_MC8.u, pomc8);

# ----------- Plot separately
begin
plot(u"hr", u"°C", xlim=(0, 70), ylim=(-36, -31), legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3, marker=:square)
modconvtplot!(sol1_MC1, sol1_Rp1, labels=["MC" "Pikal"],marker=:square)
end
savefig(plotsdir("suc1_T.svg"))

begin
plot(u"hr", u"°C", xlim=(0, 125), ylim=(-36, -32), legend=:bottomright)
@df dat2 exptfplot!(:t, :T1, :T2, :T3, marker=:circle)
modconvtplot!(sol2_MC2, sol2_Rp2, labels=["MC" "Pikal"],marker=:circle)
end
savefig(plotsdir("suc2_T.svg"))

begin
blankplot_hrC(xlim=(0, 20), ylim=(-36, -28), legend=:bottomright)
@df dat8 exptfplot!(:t, :T1, :T2, :T3, marker=:circle)
modconvtplot!(sol8_MC8, sol8_Rp8, labels=["MC" "Pikal"],marker=:circle)
end
savefig(plotsdir("suc8_T.svg"))

# Cross-check fits
sol1_MC2 = gen_sol_constNS_arb(opt_MC2.minimizer, pomc1)[1]
sol2_MC1 = gen_sol_constNS_arb(opt_MC1.minimizer, pomc2)[1]

# --------------- Fit to Cycle 1 and 2 simultaneously

function gennsol_Rp(fitprm, pos...)
    rets = map(pos) do po
        gen_sol_Rp(fitprm, po)
    end
    return [ret[1] for ret in rets], [ret[2] for ret in rets]
end
function objn_Rp(fitprm, pos, fitdats)
    sols, prms = gennsol_Rp(fitprm, pos...)
    objs = map(sols, fitdats) do sol, fitdat
        obj_expT(sol, fitdat)
    end
    return sum(objs)
end
obj18p(x) = objn_Rp(x, [pop1, pop2], [fitdat1, fitdat2])
obj128p(x) = objn_Rp(x, [pop1, pop2, pop8], [fitdat1, fitdat2, fitdat8])

opt18p = optimize(obj18p, opt_Rp1.minimizer, LevenbergMarquardt(), )
opt128p = optimize(obj128p, opt_Rp1.minimizer, LevenbergMarquardt(), )

sol1_Rp18, Rp1_params18 = gen_sol_Rp(opt12p.minimizer, pop1);
sol2_Rp18, Rp2_params18 = gen_sol_Rp(opt12p.minimizer, pop2);
sol8_Rp18, Rp8_params18 = gen_sol_Rp(opt12p.minimizer, pop8);
sol1_Rp128, Rp1_params128 = gen_sol_Rp(opt128p.minimizer, pop1);
sol2_Rp128, Rp2_params128 = gen_sol_Rp(opt128p.minimizer, pop2);
sol8_Rp128, Rp8_params128 = gen_sol_Rp(opt128p.minimizer, pop8);

function gennsol_constNS(fitprm, pos...)
    rets = map(pos) do po
        gen_sol_constNS_arb(fitprm, po)
    end
    return [ret[1] for ret in rets], [ret[2] for ret in rets]
end
function objnmc(fitprm, pos, fitdats; verbose=false)
    if !issorted(fitprm[5:8])
        return NaN
    end
    sols, prms = gennsol_constNS(fitprm, pos...)
    objs = map(sols, fitdats) do sol, fitdat
        obj_expT(sol, fitdat)
    end
    return sum(objs)
end
obj12mc(x) = objnmc(x, [pomc1, pomc2], [fitdat1, fitdat2])
obj18mc(x) = objnmc(x, [pomc1, pomc8], [fitdat1, fitdat8])
obj128mc = x->objnmc(x, [pomc1, pomc2, pomc8], [fitdat1, fitdat2, fitdat8])

g12mc = (opt_MC2.minimizer .+ opt_MC1.minimizer)/2
opt12mc = optimize(obj12mc, g12mc, LevenbergMarquardt(), )
g18mc = (opt_MC8.minimizer .+ opt_MC1.minimizer)/2
opt18mc = optimize(obj18mc, g18mc, LevenbergMarquardt(), )
g128mc = (opt_MC2.minimizer .+ opt_MC1.minimizer .+ opt_MC8.minimizer)/3
opt128mc = optimize(obj128mc, g128mc, LevenbergMarquardt(), )



sol1_MC12, MC1_params12 = gen_sol_constNS_arb(opt12mc.minimizer, pomc1);
sol2_MC12, MC2_params12 = gen_sol_constNS_arb(opt12mc.minimizer, pomc2);
sol8_MC12, MC8_params12 = gen_sol_constNS_arb(opt12mc.minimizer, pomc8);

sol1_MC18, MC1_params18 = gen_sol_constNS_arb(opt18mc.minimizer, pomc1);
sol2_MC18, MC2_params18 = gen_sol_constNS_arb(opt18mc.minimizer, pomc2);
sol8_MC18, MC8_params18 = gen_sol_constNS_arb(opt18mc.minimizer, pomc8);

sol1_MC128, MC1_params128 = gen_sol_constNS_arb(opt128mc.minimizer, pomc1);
sol2_MC128, MC2_params128 = gen_sol_constNS_arb(opt128mc.minimizer, pomc2);
sol8_MC128, MC8_params128 = gen_sol_constNS_arb(opt128mc.minimizer, pomc8);

begin
pl1 = blankplot_hrC(legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3, markershape=:circle)
# plot!(fitdat1)
plot!(xlim=(0, 120), ylim=(-36, -31), )
modconvtplot!(sol1_MC12, sol1_Rp12, labels=["MC" "Pikal"], markershape=:circle)
pl2 = blankplot_hrC(legend=:bottomright)
@df dat2 exptfplot!(:t, :T1, :T2, :T3, markershape=:square)
# plot!(fitdat2)
modconvtplot!(sol2_MC12, sol2_Rp12, labels=["MC" "Pikal"], markershape=:square)
plot!(xlim=(0, 120), ylim=(-36, -31), )

plot(pl1, pl2)
end
savefig(plotsdir("suc12_T.svg"))

begin
pl1 = blankplot_hrC(legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3, markershape=:circle)
# plot!(fitdat1)
plot!(xlim=(0, 60), ylim=(-36, -31), )
modconvtplot!(sol1_MC128, sol1_Rp128, labels=["MC" "Pikal"], markershape=:circle)
pl2 = blankplot_hrC(legend=:bottomright)
@df dat2 exptfplot!(:t, :T1, :T2, :T3, markershape=:square)
# plot!(fitdat2)
modconvtplot!(sol2_MC128, sol2_Rp128, labels=["MC" "Pikal"], markershape=:square)
plot!(xlim=(0, 120), ylim=(-36, -31), )

pl8 = blankplot_hrC(legend=:bottomright)
@df dat8 exptfplot!(:t, :T1, :T2, :T3, markershape=:circle)
# plot!(fitdat8)
modconvtplot!(sol8_MC128, sol8_Rp128, labels=["MC" "Pikal"], markershape=:circle)
plot!(xlim=(0, 20), ylim=(-32, -27))
plot(pl1, pl2, pl8, layout=(3,1), size=(600, 600), legend=:outerright)
end
savefig(plotsdir("suc128_T.svg"))

klogs12 = opt12mc.minimizer[5:8]
klogs18 = opt18mc.minimizer[5:8]
kstops12 = exp.(klogs12)*fit_un_arb[5]
kstops18 = exp.(klogs18)*fit_un_arb[5]
Tcheck = range(extrema(Tstops)..., length=50)
plot(Tstops, kstops12, yscale=:log10)
plot!(Tstops, kstops18, yscale=:log10)

# Willams Landel Ferry temparature dependence:...
function k_wlf(T, p)
    k0, Tref, C1, C2 = p
    return k0*exp(C1*(T-Tref)/(T-Tref-C2))
end
wlf_un = [u"s^-1", u"K", NoUnits, u"K"]
function obj_wlf(p)
    @info "check" p
    if p[1] < 0   
        return NaN
    end
    ks = k_wlf.(Tstops, [p.*wlf_un])
    logks = log.(ustrip.(u"s^-1", ks))
    sum(abs2.(logks .- klogs12))
end
p0_wlf = [1e-16, 273-85, 9, 10]
obj_wlf(p0_wlf)
opt_wlf = optimize(obj_wlf, p0_wlf, LevenbergMarquardt(), )
plot(Tstops,  k_wlf.(Tcheck, [opt_wlf.minimizer.*wlf_un]))
plot!(Tstops,  kstops12, yscale=:log10)
plot!(Tstops,  kstops18, yscale=:log10)

# ------------------------
# MYEGA temparature dependence of viscosity
function k_myega(T, p)
    k0, K, C = p
    # @info "check" K/T C/T K/T*exp(C/T)
    return k0*exp((-K/T)*exp(C/T))
end
myega_un = [u"s^-1", u"K", u"K"]
function obj_myega(p)
    if p[1] < 0 #|| p[2] < 100 #|| p[3] < 100 
        return NaN
    end
    ks = k_myega.(Tstops, [p.*myega_un])
    # @info "checkobj" ks p
    sum(abs2.(log.(ustrip.(u"s^-1", ks)) - klogs18))
end
p0_myega = [1e-10, 100, 2]
obj_myega(p0_myega)
opt_myega = optimize(obj_myega, p0_myega, LevenbergMarquardt(), )
opt_myega.minimizer
plot(Tcheck,  k_myega.(Tcheck, [opt_myega.minimizer.*myega_un]))
plot!(Tstops,  kstops12, yscale=:log10)
plot!(Tstops,  kstops18, yscale=:log10)

myega_params1 = @set MC1_params128.k_f = T->k_myega(T, opt_myega.minimizer.*myega_un)
myega_params2 = @set MC2_params128.k_f = T->k_myega(T, opt_myega.minimizer.*myega_un)
myega_params8 = @set MC8_params128.k_f = T->k_myega(T, opt_myega.minimizer.*myega_un)
sol1_myega = solve(ODEProblem(myega_params1), Rodas3(autodiff=true))
sol2_myega = solve(ODEProblem(myega_params2), Rodas3(autodiff=true))
sol8_myega = solve(ODEProblem(myega_params8), Rodas3(autodiff=true))
begin
pl1 =modconvtplot(sol1_myega)
plot!(fitdat1)
pl2 = modconvtplot(sol2_myega)
plot!(fitdat2)
pl3 = modconvtplot(sol8_myega)
plot!(fitdat8)
plot(pl1, pl2, pl3)
end

fit_un_myega = [u"cm^-1", u"cm^5", u"m^0.5", u"m/s", u"s^-1", u"K", u"K"]
function gen_sol_constNS_myega(fitprmlog, po; verbose=false)
    new_params = deepcopy(po)
    fitprm = exp.(fitprmlog) # Constrain all to be positive by taking exp
    verbose && @info "fitprm" fitprm
    @reset new_params.N = RampedVariable(fitprm[1]*fit_un_myega[1])
    @reset new_params.S = RampedVariable(fitprm[2]*fit_un_myega[2])
    @reset new_params.C  = fitprm[3]*fit_un_myega[3]
    @reset new_params.R0 = fitprm[4]*fit_un_myega[4]
    @reset new_params.k_f = T->k_myega(T, fitprm[5:7].*fit_un_myega[5:7])
    prob = ODEProblem(new_params)
    sol = solve(prob, Rodas3(autodiff=true))
    return sol, new_params
end
function gennsol_constNS_myega(fitprmlog, pos...)
    rets = map(pos) do po
        gen_sol_constNS_myega(fitprmlog, po)
    end
    return [ret[1] for ret in rets], [ret[2] for ret in rets]
end
function objnmyega(fitprmlog, pos, fitdats; verbose=true)
    verbose && @info "fitprm" exp.(fitprmlog)
    sols, prms = gennsol_constNS_myega(fitprmlog, pos...)
    if any([sol.retcode != ReturnCode.Terminated for sol in sols])
        verbose && @info "unfinished drying" [sol.retcode for sol in sols]
        return NaN
    end
    objs = map(sols, fitdats) do sol, fitdat
        obj_expT(sol, fitdat)
    end
    verbose && @info "obj" objs
    return sum(objs)
end
obj128myega(x) = objnmyega(x, [MC1_params12, MC2_params12, MC8_params128], [fitdat1, fitdat2, fitdat8])
p0_myegafitlog = [opt128mc.minimizer[1:4]..., log.(opt_myega.minimizer)...]
@time obj128myega(p0_myegafitlog)
opt128myega = optimize(obj128myega, p0_myegafitlog, LevenbergMarquardt(), )
@time obj128myega(opt128myega.minimizer)

sol1_myega128, myega1_params128 = gen_sol_constNS_myega(opt128myega.minimizer, MC1_params128);
sol2_myega128, myega2_params128 = gen_sol_constNS_myega(opt128myega.minimizer, MC2_params128);
sol8_myega128, myega8_params128 = gen_sol_constNS_myega(opt128myega.minimizer, MC8_params128);

begin
pl1 = blankplot_hrC(legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3, markershape=:circle)
plot!(title="Tsh = $(Tsh1(Inf*u"s")|>u"°C")")
plot!(xlim=(0, 120), ylim=(-36, -31), )
modconvtplot!(sol1_myega128, sol1_Rp128, labels=["MC" "Pikal"], markershape=:circle)
pl2 = blankplot_hrC(legend=:bottomright)
@df dat2 exptfplot!(:t, :T1, :T2, :T3, markershape=:square)
modconvtplot!(sol2_myega128, sol2_Rp128, labels=["MC" "Pikal"], markershape=:square)
plot!(xlim=(0, 120), ylim=(-38, -31), )
plot!(title="Tsh = $(Tsh2(Inf*u"s")|>u"°C")")

pl8 = blankplot_hrC(legend=:bottomright)
@df dat8 exptfplot!(:t, :T1, :T2, :T3, markershape=:utriangle)
# plot!(fitdat8)
modconvtplot!(sol8_myega128, sol8_Rp128, labels=["MC" "Pikal"], markershape=:utriangle)
plot!(xlim=(0, 20), ylim=(-40, -20))
plot!(title="Tsh = $(Tsh8(Inf*u"s")|>u"°C")")
plot(pl2, pl1, pl8, layout=(3,1), size=(600, 600), legend=:outerright)
end
savefig(plotsdir("suc128_myega_T.svg"))

# -------------------------------------
# playground temparature dependence of viscosity
function k_play(T, p)
    k0, kf, K, C = p
    rat = k0/kf
    # x = ustrip(u"K", T-K)
    # x = ustrip(u"K", T-K)
    x = (T-K)/C
    # sig = (tanh(x)+1)/2
    # return k0/(tanh(-x)+1) # This has the right behavior...
    return k0*(exp(x) + exp(-x))/((1+rat)*exp(-x)+rat*exp(x)) # but this way it can avoid numerical issues
end
play_un = [u"s^-1", u"s^-1", u"K", u"K"]
p0_play = [1e-17, 1e-10, 230, 0.8]
scatter(Tcheck,  k_play.(Tcheck, [p0_play .* play_un]), yscale=:log10)
plot!(Tstops,  kstops12, yscale=:log10)
plot!(Tstops,  kstops18, yscale=:log10)
function obj_play(p)
    if p[1] < 0 || p[2] < 0 #|| p[3] < 100 
        return NaN
    end
    ks = k_play.(Tstops, [p.*play_un])
    # @info "check" ks kstops18
    sum(abs2.(log.(ustrip.(u"s^-1", ks)) - klogs18))
end
obj_play(p0_play)
opt_play = optimize(obj_play, p0_play, LevenbergMarquardt(), )
plot(Tcheck,  k_play.(Tcheck, [opt_play.minimizer.*play_un]))
plot!(Tstops,  kstops12, yscale=:log10)
plot!(Tstops,  kstops18, yscale=:log10)

play_params1 = @set MC1_params128.k_f = T->k_play(T, opt_play.minimizer.*play_un)
play_params2 = @set MC2_params128.k_f = T->k_play(T, opt_play.minimizer.*play_un)
sol1_play = solve(ODEProblem(play_params1), Rodas3(autodiff=true))
sol2_play = solve(ODEProblem(play_params2), Rodas3(autodiff=true))
poreevolplot(sol1_play)
poreevolplot!(sol2_play)

fit_un_play = [u"cm^-1", u"cm^5", u"m^0.5", u"m/s", u"s^-1", u"s^-1", u"K", u"K"]
function gen_sol_constNS_play(fitprmlog, po; verbose=false)
    new_params = deepcopy(po)
    fitprm = exp.(fitprmlog) # Constrain all to be positive by taking exp
    verbose && @info "fitprm" fitprm
    @reset new_params.N = RampedVariable(fitprm[1]*fit_un_play[1])
    @reset new_params.S = RampedVariable(fitprm[2]*fit_un_play[2])
    @reset new_params.C  = fitprm[3]*fit_un_play[3]
    @reset new_params.R0 = fitprm[4]*fit_un_play[4]
    @reset new_params.k_f = T->k_play(T, fitprm[5:8].*fit_un_play[5:8])
    prob = ODEProblem(new_params)
    sol = solve(prob, Rodas3(autodiff=true))
    return sol, new_params
end
function gennsol_constNS_play(fitprmlog, pos...)
    rets = map(pos) do po
        gen_sol_constNS_play(fitprmlog, po)
    end
    return [ret[1] for ret in rets], [ret[2] for ret in rets]
end
function objnplay(fitprmlog, pos, fitdats; verbose=true)
    verbose && @info "fitprm" exp.(fitprmlog)
    if exp(fitprmlog[7]) < 225 || exp(fitprmlog[7]) > 250
        verbose && @info "temperature unreasonable" exp(fitprmlog[7])
        return NaN 
    end
    sols, prms = gennsol_constNS_play(fitprmlog, pos...)
    if any([sol.retcode != ReturnCode.Terminated for sol in sols])
        verbose && @info "unfinished drying" [sol.retcode for sol in sols]
        return NaN
    end
    objs = map(sols, fitdats) do sol, fitdat
        obj_expT(sol, fitdat)
    end
    verbose && @info "obj" objs
    return sum(objs)
end
p0_playfitlog = [opt128mc.minimizer[1:4]..., log.(opt_play.minimizer)...]
# p0_playfitlog = log.(p0_playfit)

obj12play(x) = objnplay(x, [MC1_params12, MC2_params12], [fitdat1, fitdat2])
obj18play(x) = objnplay(x, [MC1_params12, MC8_params128], [fitdat1, fitdat8])
obj128play(x) = objnplay(x, [MC1_params128, MC2_params128, MC8_params128], [fitdat1, fitdat2, fitdat8])
obj12play(p0_playfitlog)

opt12play = optimize(obj12play, p0_playfitlog, LevenbergMarquardt(), )
opt18play = optimize(obj18play, p0_playfitlog, LevenbergMarquardt(), )
opt128play = optimize(obj128play, p0_playfitlog, LevenbergMarquardt(), )

exp.(opt128play.minimizer) .*fit_un_play

(sol1_play18, sol2_play18, sol8_play18), playprms18 = gennsol_constNS_play(opt18play.minimizer, MC1_params128, MC2_params128, MC8_params128);
(sol1_play128, sol2_play128, sol8_play128), playprms128 = gennsol_constNS_play(opt128play.minimizer, MC1_params128, MC2_params128, MC8_params128);
# (sol1_playg, sol2_playg, sol8_playg), playprmsg = gennsol_constNS_play(p0_playfitlog, MC1_params2, MC2_params2, MC8_params128);
poreevolplot( sol1_play128)
poreevolplot!(sol2_play128)
poreevolplot!(sol8_play128)

begin
pl1 = blankplot_hrC(legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3, markershape=:circle)
# plot!(fitdat1)
plot!(title="Tsh = $(Tsh1(Inf*u"s")|>u"°C")")
plot!(xlim=(0, 120), ylim=(-36, -31), )
# modconvtplot!(sol1_playg, sol1_Rp18, labels=["MC" "Pikal"], markershape=:circle)
modconvtplot!(sol1_play18, sol1_Rp128, labels=["MC" "Pikal"], markershape=:circle)
# modconvtplot!(sol1_play128, sol1_Rp18, labels=["MC" "Pikal"], markershape=:circle)
pl2 = blankplot_hrC(legend=:bottomright)
@df dat2 exptfplot!(:t, :T1, :T2, :T3, markershape=:square)
# plot!(fitdat2)
# modconvtplot!(sol2_playg, sol2_Rp18, labels=["MC" "Pikal"], markershape=:square)
modconvtplot!(sol2_play18, sol2_Rp128, labels=["MC" "Pikal"], markershape=:square)
# modconvtplot!(sol2_play128, sol2_Rp18, labels=["MC" "Pikal"], markershape=:square)
plot!(xlim=(0, 120), ylim=(-38, -31), )
plot!(title="Tsh = $(Tsh2(Inf*u"s")|>u"°C")")

pl8 = blankplot_hrC(legend=:bottomright)
@df dat8 exptfplot!(:t, :T1, :T2, :T3, markershape=:utriangle)
# plot!(fitdat8)
# modconvtplot!(sol8_playg, sol8_Rp18, labels=["MC" "Pikal"], markershape=:utriangle)
modconvtplot!(sol8_play18, sol8_Rp128, labels=["MC" "Pikal"], markershape=:utriangle)
# modconvtplot!(sol8_play128, sol8_Rp18, labels=["MC" "Pikal"], markershape=:utriangle)
plot!(xlim=(0, 20), ylim=(-40, -20))
plot!(title="Tsh = $(Tsh8(Inf*u"s")|>u"°C")")
plot(pl2, pl1, pl8, layout=(3,1), size=(600, 600), legend=:outerright)
end
savefig(plotsdir("suc128_play_T.svg"))

begin
plot(u"°C", u"s^-1", legend=:bottomright, xlabel="T", ylabel="k\\  ", unitformat=latexsquareunitlabel)
# kineticplot!(MC1_params,    T=Tcheck, label="Interp fit to case 1")
# kineticplot!(MC2_params,    T=Tcheck, label="Interp fit to case 2")
kineticplot!(MC1_params12,  T=Tcheck, label="Interp fit to 1&2")
kineticplot!(MC1_params18,  T=Tcheck, label="Interp fit to 1&8")
kineticplot!(MC1_params128, T=Tcheck, label="Interp fit to 1&2&8")
kineticplot!(myega1_params128  , T=Tcheck, label="MYEGA fit to 1&2&8")
kineticplot!(playprms18[1]  , T=Tcheck, label="Tanh fit to 1&8")
kineticplot!(playprms128[1]  , T=Tcheck, label="Tanh fit to 1&2&8")

# kineticplot!(play1_params, T=Tcheck, label="3-param activation to 1&2")
end
savefig(plotsdir("kinetic.svg"))



l = range(0u"cm", hf0, length=50)
Rp1 = Rp1_params.Rp.(l)
Rp2 = Rp2_params.Rp.(l)
Rp12 = Rp2_params2.Rp.(l)
Rp12 = Rp2_params2.Rp.(l)

default(:markersize, 7)
begin
colors = [:red :blue :purple :brown]
rpattr = (marker=:utriangle, seriestype=:samplemarkers, linestyle=:dash)
# pl2 = poreevolplot(sol1_MC1, label="Fit 1, Cycle 1", c=colors[1], marker=:circle)
# poreevolplot!(sol2_MC1, label="Fit 1, Cycle 2", c=colors[1], marker=:square)
# poreevolplot!(sol1_MC2, label="Fit 2, Cycle 1", c=colors[2], marker=:circle)
# poreevolplot!(sol2_MC2, label="Fit 2, Cycle 2", c=colors[2], marker=:square)
# poreevolplot!(sol1_MC12, label="Fit both, Cycle 1", c=colors[3], marker=:circle)
# poreevolplot!(sol2_MC12, label="Fit both, Cycle 2", c=colors[3], marker=:square)
poreevolplot(sol1_play18, label="PBM 1&8, Cycle 1", c=colors[3], marker=:circle)
poreevolplot!(sol2_play18, label="PBM 1&8, Cycle 2", c=colors[3], marker=:square)
plot!(l, Rp1, subplot=4, label="3PP 1"; c=colors[1], rpattr...)
plot!(l, Rp2, subplot=4, label="3PP 2"; c=colors[2], rpattr...)
plot!(l, Rp12, subplot=4, label="3PP 1&2"; c=colors[3], rpattr...)
plot!(size=(800, 400))
end
savefig(plotsdir("suc12_pores.svg"))
plot(pl1, pl2, size=(1000, 400))

# -----------------
# Inspect fit results in Rp space rather than temperature space
function calc_lRp(dat)
    minteg = @df dat cumsum(:mfr[1:end-1].*(:t[2:end]-:t[1:end-1]))
    insert!(minteg, 1, 0u"g")
    errfac = minteg[end] / (fillvol*(ρ_solution-c_solid))
    ld = @. errfac*minteg/(A_p*(ρ_solution-c_solid)) |> u"cm"
end
ld1 = calc_l(dat1)
ld2 = calc_l(dat2)

function calc_lRp(sol, prm)
    l = prm.hf0 .- sol[1,:]*u"cm"
    Rp_mc = [FulbrightPores.lyo_mc_dae!(fill(0.0, 4), u, prm, t) for (u,t) in zip(sol.u, sol.t)]
    return l, Rp_mc
end
l1, Rp_mc1 = calc_lRp(sol1_play18, playprms18[1])
l8, Rp_mc8 = calc_lRp(sol8_play18, playprms18[2])
Rp_p1 = Rp1_params.Rp.(ld1)
Rp_p2 = Rp2_params.Rp.(ld2)

begin
pl1 = plot(u"cm", u"cm^2*Torr*hr/g", ylims=(0, 20), xlabel=L"h_d", ylabel=L"R_p", unitformat=(x,y)->"$x [$(latexify(y))]")
plot!(ld1, dat1.Rp, linealpha=0.8, label="TDLAS 1", c=:red)
plot!(l1, Rp_mc1,   linealpha=0.8, label="PBM 1&8, Cycle 1 ", c=:red , marker=:circle, seriestype=:samplemarkers, step=3)
plot!(ld1, Rp_p1,   linealpha=0.8, label="3PP 1", c=:red , marker=:utriangle, seriestype=:samplemarkers, step=1000)
plot!(ld1, Rp1_params18.Rp.(ld1),   linealpha=0.8, label="3PP 1&8", c=:blue, marker=:square, seriestype=:samplemarkers, step=1000)
pl2 = plot(u"cm", u"cm^2*Torr*hr/g", ylims=(0, 20), xlabel=L"h_d", ylabel=L"R_p", unitformat=(x,y)->"$x [$(latexify(y))]")
plot!(ld2, dat2.Rp, linealpha=0.8, label="TDLAS 2", c=:blue)
plot!(l2, Rp_mc2,   linealpha=0.8, label="PBM 1&8, Cycle 2 ", c=:blue, marker=:circle, seriestype=:samplemarkers, step=3)
plot!(ld2, Rp_p2,   linealpha=0.8, label="3PP 2", c=:blue, marker=:utriangle, seriestype=:samplemarkers, step=1000)
plot!(ld2, Rp1_params18.Rp.(ld2),   linealpha=0.8, label="3PP 1&8", c=:blue, marker=:square, seriestype=:samplemarkers, step=1000)
plot(pl1, pl2)
end
savefig(plotsdir("suc12_Rp.svg"))

function namefitparam_arr(fitlog, label)
    fit = exp.(fitlog)
    nt = (label = label,
        N =  fit[1]*fit_un_arr[1],
        S =  fit[2]*fit_un_arr[2],
        C =  fit[3]*fit_un_arr[3],
        R0 = fit[4]*fit_un_arr[4],
        k0 = fit[5]*fit_un_arr[5],
        E =  fit[6]*fit_un_arr[6],
    )
end
function namefitparam_arb(fitlog, label)
    fit = exp.(fitlog)
    nt = (label = label,
        N =  fit[1]*fit_un_arb[1],
        S =  fit[2]*fit_un_arb[2],
        C =  fit[3]*fit_un_arb[3],
        R0 = fit[4]*fit_un_arb[4],
        k1 = fit[5]*fit_un_arb[5],
        k2 = fit[6]*fit_un_arb[5],
        k3 = fit[7]*fit_un_arb[5],
        k4 = fit[8]*fit_un_arb[5],
        T1 = Tstops[1],
        T2 = Tstops[2],
        T3 = Tstops[3],
        T4 = Tstops[4],
    )
end
function namefitparam_play(fitlog, label)
    fit = exp.(fitlog)
    nt = (label = label,
        N =  fit[1]*fit_un_play[1],
        S =  fit[2]*fit_un_play[2],
        C =  fit[3]*fit_un_play[3],
        R0 = fit[4]*fit_un_play[4],
        k0 = fit[5]*fit_un_play[5],
        kf = fit[6]*fit_un_play[6],
        Tact = fit[7]*fit_un_play[7],
        θ  = fit[8]*fit_un_play[8],
    )
end
fits_arb = Optim.minimizer.([opt_MC1, opt_MC2, opt12mc, opt18mc, opt128mc]) 
fits_play = Optim.minimizer.([opt12play, opt18play, opt128play]) 
rep1 = Table(map(namefitparam_arb, fits_arb, ["1", "2", "1&2", "1&8", "1&2&8"]))
rep2 = Table(map(namefitparam_play, fits_play, ["1&2", "1&8", "1&2&8"]))

# safesave(datadir("exp_pro", "fit12.html"), pretty_table(HTML, report, header = ["Fit", "N", "S", "k0", "E", "C", "R0"]))
safesave(datadir("exp_pro", "fit_arb.html"), pretty_table(HTML, rep1, header = ["Fit", "N", "S", "C", "R0", "k1", "k2", "k3", "k4", "T1", "T2", "T3", "T4"]))
safesave(datadir("exp_pro", "fit_play.html"), pretty_table(HTML, rep2, header = ["Fit", "N", "S", "C", "R0", "k0", "kf", "Tact", "θ"]))
# pretty_table(report, header = ["Fit", "N", "S", "k0", "E", "C", "R0"])

# ---------------------
# Validation


pop3 = deepcopy(Rp1_params18)
@reset pop3.pch = pch3
@reset pop3.Tsh = Tsh3
sol3_Rp18 = solve(ODEProblem(pop3), Rodas3());
po_PBM3 = deepcopy(playprms18[1])
@reset po_PBM3.pch = pch3
@reset po_PBM3.Tsh = Tsh3
sol3_PBM18 = solve(ODEProblem(po_PBM3), Rodas3());

pop5 = deepcopy(Rp1_params18)
@reset pop5.pch = pch5
@reset pop5.Tsh = Tsh5
sol5_Rp18 = solve(ODEProblem(pop5), Rodas3());
po_PBM5 = deepcopy(playprms18[1])
@reset po_PBM5.pch = pch5
@reset po_PBM5.Tsh = Tsh5
sol5_PBM18 = solve(ODEProblem(po_PBM5), Rodas3());

# Training plot
begin
pl1 = blankplot_hrC(legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3, markershape=:circle)
plot!(xlim=(0, 70), ylim=(-36, -31), )
modconvtplot!(sol1_play18, sol1_Rp18, labels=["PBM 1&8" "3PP 1&8"], markershape=:circle)
plot!(title="Tsh = $(Tsh1(Inf*u"s")|>u"°C")", xlabel="")
pl8 = blankplot_hrC(legend=:bottomright)
@df dat8 exptfplot!(:t, :T1, :T2, :T3, markershape=:utriangle)
modconvtplot!(sol8_play18, sol8_Rp18, labels=["PBM 1&8" "3PP 1&8"], markershape=:utriangle)
plot!(xlim=(0, 30), ylim=(-32, -26))
plot!(title="Tsh = $(Tsh8(Inf*u"s")|>u"°C")", ylabel="")
plot(pl1, pl8, layout=(3,1), size=(600, 400), legend=:outerright)
end
savefig(plotsdir("suc128_train_T.svg"))

# Validation plot

begin
pl2 = blankplot_hrC()
@df dat2 exptfplot!(:t, :T1, :T2, :T3, markershape=:utriangle)
modconvtplot!(sol2_play18, sol2_Rp18, labels=["PBM 1&8" "3PP 1&8" ], markershape=:circle)
plot!(title="Tsh = $(Tsh2(Inf*u"s")|>u"°C")", xlabel="", ylabel="")
pl3 = blankplot_hrC()
@df dat3 exptfplot!(:t, :T1, :T2, :T3, markershape=:circle)
modconvtplot!(sol3_PBM18, sol3_Rp18, labels=["PBM 1&8" "3PP 1&8" ], markershape=:circle)
plot!(title="Tsh = $(Tsh3(Inf*u"s")|>u"°C")", xlabel="")
pl5 = blankplot_hrC()
@df dat5 exptfplot!(:t, :T1, :T2, :T3, markershape=:square)
modconvtplot!(sol5_PBM18, sol5_Rp18, labels=[ "PBM 1&8" "3PP 1&8"], markershape=:square)
plot!(title="Tsh = $(Tsh5(Inf*u"s")|>u"°C")", ylabel="")
plot(pl2, pl3, pl5, layout=(3,1))
plot!(ylim=(-37, -32), legend=:outerright)
end
savefig(plotsdir("suc18train_235val_T.svg"))


