using DrWatson
@quickactivate :FulbrightPores
using SavitzkyGolay

const OP = OptimizationProblem
const NLSP = NonlinearLeastSquaresProblem
const BFGS = Optim.BFGS
using LineSearches
alg = BFGS(linesearch=LineSearches.BackTracking())
const LP = LyoPronto
pythonplot()

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
optKv = solve(NLSP((x,p,t)->fitKv(x)[1], [10.0, 1, 1]), LevenbergMarquardt())
K_shf = fitKv(optKv.u)[2]

# --------------------------
# Cycle 1
# Labeled as case 4 in Table 8

dat1raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle1.csv"), Table)
dat1 = map(simtra_rename, dat1raw)
# dat1 = dat1[dat1.T2 .< 250u"K"]

pch1_pir_der2 = savitzky_golay(dat1.pch_pir, 391, 3, deriv=2, rate=1/u"minute").y
t_end1 = dat1.t[argmax(pch1_pir_der2[10u"hr" .< dat1.t .< 85u"hr"])] + 10u"hr"

pch1 = RampedVariable(100u"mTorr")
@df dat1 plot(:t, :pch_cm, ylim=(0, 200))
@df dat1 plot(:t, :pch_pir, ylim=(0, 200))
plot!(pch1)
tendplot!(t_end1)

# There's a slight delay before shelf temperature ramps
Tsh1 = RampedVariable(uconvert.(u"K", [-40, -40, -25]u"°C"), [1u"K/minute", 0.5u"K/minute"], [23u"minute"])
@df dat1 plot(:t, :Tsh, xlim=(0, 2))
plot!(Tsh1)

fitdat1 = PrimaryDryFit(collect(dat1.t), (
    collect(filter(x->(!isnan(x) && x<-32.5u"°C"), dat1.T1)),
    collect(filter(!isnan, dat1.T2[dat1.t .< 60u"hr"])),
    # collect(filter(!isnan, dat1.T3)),
    ),
    t_end1
)
plot(fitdat1)

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

pch2_pir_der2 = savitzky_golay(dat2.pch_pir, 391, 3, deriv=2, rate=1/u"minute").y
t_end2 = dat2.t[argmax(pch2_pir_der2[10u"hr" .< dat2.t .< 120u"hr"])] + 10u"hr"
pch2 = RampedVariable(100u"mTorr")
@df dat2 plot(:t, :pch_cm, ylim=(0, 200))
@df dat2 plot(:t, :pch_pir, ylim=(0, 200))
plot!(pch2)
tendplot!(t_end2)

fitdat2 = @df dat2 PrimaryDryFit(collect(:t), (
    # collect(filter(!isnan, :T1)),
    collect(filter(!isnan, :T2[:t .< 85u"hr"])),
    collect(filter(!isnan, :T3)),),
    t_end2)
plot(fitdat2)

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

pch8_pir_der2 = savitzky_golay(dat8.pch_pir, 391, 3, deriv=2, rate=1/u"minute").y
t_end8 = dat8.t[argmax(pch8_pir_der2[10u"hr" .< dat8.t .< 41u"hr"])] + 10u"hr"
pch8 = RampedVariable(100u"mTorr")
@df dat8 plot(:t, :pch_cm, ylim=(0, 200))
@df dat8 plot(:t, :pch_pir, ylim=(0, 200))
plot!(pch8)
tendplot!(t_end8)

Tsh8 = RampedVariable(uconvert.(u"K", [-40, -40, -10]u"°C"), [1u"K/minute", 0.5u"K/minute"], [4u"minute"])
@df dat8 plot(:t, :Tsh_i, xlim=(0, 2))
plot!(Tsh8)

fitdat8 = @df dat8 PrimaryDryFit(collect(:t), (
    collect(filter(!isnan, :T1)),
    collect(filter(!isnan, :T2[:t .< 17u"hr"])),
    collect(filter(!isnan, :T3[:t .< 17u"hr"])),),
    t_end8)
plot(fitdat8)


# -------------
# Parameter setup

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

const Rp_un = (u"cm^2*Torr*hr/g", u"cm*Torr*hr/g", u"1/cm")
function gen_sol_Rp(a123, po)
    new_params = @set po.Rp = RpFormFit((a123.*Rp_un)...)
    prob = ODEProblem(new_params; tspan=(0.0, 1000.0))
    sol = solve(prob, Rodas3())
    return sol
end
function obj_Rp(x, p)
    if any(x .< 0)
        return NaN
    end
    sol = gen_sol_Rp(x, p[1])
    return obj_expT(sol, p[2], tweight=1e-4)
end
obj_Rpf = OptimizationFunction(obj_Rp, AutoForwardDiff())
@time obj_Rpf([0.8, 90, 10], (pop8, fitdat8))

function gennsol_Rp(fitprm, pos)
    return [gen_sol_Rp(fitprm, po) for po in pos]
end
function objn_Rp(fitprm, pfs)
    if any(fitprm .< 0)
        return NaN
    end
    pos, fitdats = pfs
    sols = gennsol_Rp(fitprm, pos)
    objs = map(sols, fitdats) do sol, fitdat
        obj_expT(sol, fitdat, tweight=1e-3)
    end
    return sum(objs)
end
obj_npf = OptimizationFunction(objn_Rp, AutoForwardDiff())
# @time ForwardDiff.gradient(x->obj_Rpf(x, (pop8, fitdat8)), [0.8, 90, 1])

@time opt_Rp1 = solve(OP(obj_Rpf, [0.8, 90, 1], (pop1, fitdat1)),alg, g_tol = 1e-3)
opt2p = solve(OP(obj_Rpf, opt_Rp1.u, (pop2, fitdat2)), alg; g_tol=1e-3)
opt8p = solve(OP(obj_Rpf, opt_Rp1.u, (pop8, fitdat8)), alg; g_tol=1e-3)

sol1_Rp1 = gen_sol_Rp(opt_Rp1.u, pop1);
sol2_Rp2 = gen_sol_Rp(opt2p.u, pop2);
sol8_Rp8 = gen_sol_Rp(opt8p.u, pop8);
# plot(u"hr", u"°C", xlim=(0, 60), ylim=(-40, -30))
blankplot_hrC(ylim=(-40, -30))
modconvtplot!(sol1_Rp1)
# @df dat1 exptfplot!(:t, :T1, :T2, :T3)
plot!(fitdat1)

pass18p =  ([pop1, pop2], [fitdat1, fitdat2]);
pass128p = ([pop1, pop2, pop8], [fitdat1, fitdat2, fitdat8])
@time opt18p = solve(OP(obj_npf, opt_Rp1.u, pass18p), alg; g_tol=1e-2, show_trace=true)
@time opt128p = solve(OP(obj_npf, opt_Rp1.u, pass128p), alg; g_tol=1e-2)

sol1_Rp18 = gen_sol_Rp(opt18p.u, pop1);
sol2_Rp18 = gen_sol_Rp(opt18p.u, pop2);
sol8_Rp18 = gen_sol_Rp(opt18p.u, pop8);
sol1_Rp128 = gen_sol_Rp(opt128p.u, pop1);
sol2_Rp128 = gen_sol_Rp(opt128p.u, pop2);
sol8_Rp128 = gen_sol_Rp(opt128p.u, pop8);
Rp1_params18 = sol1_Rp18.prob.p
Rp1_params = sol1_Rp1.prob.p
Rp2_params = sol2_Rp2.prob.p
Rp8_params = sol8_Rp8.prob.p

# -------------- Set up parameters

# typical diameter: 50 μm
d_typical = 50u"μm"
# typical total volume: 
V_typical = 4/3*π*d_typical^3
N_typical = fillvol/V_typical

# Good guesses for Knudsen
N = RampedVariable(N_typical / hf0)
S = RampedVariable(V_typical^2*N_typical/hf0)
C = 2e-1u"m^0.5"
Rs = 0.05u"cm^2*Torr*hr/g"

# Good guesses for viscous
# N = RampedVariable(N_typical / hf0) 
# S = RampedVariable(1e0u"cm^5")
# C = 2e2u"m"
# Rs = 0.01u"cm^2*Torr*hr/g"

k0 = 1e-14u"1/s"

pb_mc = (
    (Rp_g, hf0, c_solid, ρ_solution),
    (K_shf, A_v, A_p),
    (pch1, Tsh1),
    (ConstPhysProp(k0), C, N, S, Rs)
)
pomc1 = ParamObjMC(pb_mc)
pomc2 = deepcopy(pomc1)
@reset pomc2.pch = pch2
@reset pomc2.Tsh = Tsh2
pomc8 = deepcopy(pomc1)
@reset pomc8.pch = pch8
@reset pomc8.Tsh = Tsh8

# const fit_un_mult = (u"cm^-1", u"cm^5", u"m^0.5", u"m/s", fill(u"s^-1", 4)...)
const fit_mul_porous = (N(0), S(0), C, Rs)
const fit_un_interp = (fill(u"s^-1", 4)...)
function gen_sol_constNS_arb(fitprm, po, k_f_gen; verbose=false)
    new_params = deepcopy(po)
    # fitprm = exp.(fitprmlog) # Constrain all to be positive by taking exp
    verbose && @info "fitprm" fitprm
    @reset new_params.N = RampedVariable(fitprm[1]*fit_un_mult[1])
    @reset new_params.S = RampedVariable(fitprm[2]*fit_un_mult[2])
    @reset new_params.C  = fitprm[3]*fit_un_mult[3]
    @reset new_params.Rs = fitprm[4]*fit_un_mult[4]
    @reset new_params.k_f = k_f_gen(fitprm[5:end])
    prob = ODEProblem(new_params; tspan=(0.0, 1000.0))
    sol = solve(prob, Rodas3())
    return sol
end
function gennsol_constNS(fitprm, pos, k_f_gen)
    sols = map(pos) do po
        gen_sol_constNS_arb(fitprm, po, k_f_gen)
    end
    return sols
end

# -----------------------------
# Fit k by interpolating at four values
Tstops = [-45.0, -38, -32, -28]u"°C" .|> u"K"
fit_mul_interp = (k0, k0, k0)
function k_f_interp(logk)
    # Lin interp in log space
    spll = LinearInterpolation(logk, Tstops, extrapolate=true)
    spl(T) = exp(spll(T))*fit_mul_interp
    # spl = LinearInterpolation(exp.(logk).*fit_un_interp, Tstops, extrapolate=true)
end
function obj_interp(x,p)
    ks = x[5:8]
    if ~(ks[1] <= ks[2] <= ks[3] <= ks[4])
        return NaN
    end
    sol, _ = gen_sol_constNS_arb(x, p[1], k_f_interp)
    return obj_expT(sol, p[2], tweight=0)
end
function obj_ninterp(fitprm, pfs; verbose=false)
    pos, fitdats = pfs
    if ~(fitprm[5] < fitprm[6] < fitprm[7] < fitprm[8])
        return NaN
    end
    sols = gennsol_constNS(fitprm, pos, k_f_interp)
    objs = map(sols, fitdats) do sol, fitdat
        obj_expT(sol, fitdat, tweight=1e-3)
    end
    return sum(objs)
end
objf_interp = OptimizationFunction(obj_interp, AutoForwardDiff())
objnf_interp = OptimizationFunction(obj_ninterp, AutoForwardDiff())
lb_MC = vcat(fill(0.0, 4), fill(-20., 4))
ub_MC = vcat(fill(Inf, 4), fill(20., 4))

@time opt_MC1 = solve(OP(obj_MCf, fitguess, (pomc1, fitdat1), lb=lb_MC, ub=ub_MC), alg;show_trace=true, g_tol = 1e-2, maxtime=10)

sol1_MC1 = gen_sol_constNS_arb(opt_MC1.u, pomc1, k_f_interp);
MC1_params1 = sol1_MC1.prob.p
begin
plot(u"hr", u"°C", xlim=(0, 120), ylim=(-36, -31), legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3, marker=:square)
modconvtplot!(sol1_MC1, sol1_Rp1, labels=["MC" "Pikal"],marker=:square)
end
savefig(plotsdir("suc1_T.svg"))


@time opt_MC8 = solve(OP(obj_MCf, opt_MC1.u, (pomc8, fitdat8), lb=lb_MC, ub=ub_MC), alg; show_trace=true, g_tol=1e-2)

sol8_MC8 = gen_sol_constNS_arb(opt_MC8.u, pomc8, k_f_interp);
MC8_params = sol8_MC8.prob.p
begin
blankplot_hrC(xlim=(0, 45), ylim=(-32, -26), legend=:bottomright)
# @df dat8 exptfplot!(:t, :T1, :T2, :T3, marker=:circle)
plot!(fitdat8)
modconvtplot!(sol8_MC8, sol8_Rp8, labels=["MC" "Pikal"],marker=:circle)
end
savefig(plotsdir("suc8_T.svg"))

g18mc = (opt_MC8.u .+ opt_MC1.u)/2
# opt18mc = optimize(obj18mc, g18mc, NelderMead(), Optim.Options(g_tol = 1e-3))
@time opt18mc = solve(OP(obj_ninterp, g18mc, ([pomc1, pomc8], [fitdat1, fitdat8]), lb=lb_MC, ub=ub_MC), alg, g_tol=1e-3, )
# opt128mc = optimize(obj128mc, g128mc, NelderMead(), Optim.Options(g_tol = 1e-3))
@time opt128mc = solve(OP(obj_ninterp, g18mc, ([pomc1, pomc2, pomc8], [fitdat1, fitdat2, fitdat8]), lb=lb_MC, ub=ub_MC), alg, g_tol=1e-3)




sol1_MC18  = gen_sol_constNS_arb(opt18mc.u, pomc1, k_f_interp);
sol2_MC18  = gen_sol_constNS_arb(opt18mc.u, pomc2, k_f_interp);
sol8_MC18  = gen_sol_constNS_arb(opt18mc.u, pomc8, k_f_interp);
MC1_params18 = sol1_MC18.prob.p
MC2_params18 = sol2_MC18.prob.p
MC8_params18 = sol8_MC18.prob.p

sol1_MC128  = gen_sol_constNS_arb(opt128mc.u, pomc1, k_f_interp);
sol2_MC128  = gen_sol_constNS_arb(opt128mc.u, pomc2, k_f_interp);
sol8_MC128  = gen_sol_constNS_arb(opt128mc.u, pomc8, k_f_interp);
MC1_params128 = sol1_MC128.prob.p
MC2_params128 = sol2_MC128.prob.p
MC8_params128 = sol8_MC128.prob.p

# begin
# pl1 = blankplot_hrC(legend=:bottomright)
# # @df dat1 exptfplot!(:t, :T1, :T2, :T3, markershape=:circle)
# plot!(fitdat1)
# plot!(xlim=(0, 100), ylim=(-36, -31), )
# modconvtplot!(sol1_MC18, sol1_Rp18, labels=["MC" "Pikal"], markeralpha=0)
# pl8 = blankplot_hrC(legend=:bottomright)
# # @df dat8 exptfplot!(:t, :T1, :T2, :T3, markershape=:square)
# plot!(fitdat8, markershap=:square)
# modconvtplot!(sol8_MC18, sol8_Rp18, labels=["MC" "Pikal"], markeralpha=0)
# plot!(xlim=(0, 45), ylim=(-32, -26), )

# plot(pl1, pl8, layout=(2,1), legend=:outerright)
# end
# savefig(plotsdir("suc18_T.svg"))

# begin
# pl1 = blankplot_hrC(legend=:bottomright)
# @df dat1 exptfplot!(:t, :T1, :T2, :T3, markershape=:circle)
# # plot!(fitdat1)
# plot!(xlim=(0, 60), ylim=(-36, -31), )
# modconvtplot!(sol1_MC128, sol1_Rp128, labels=["MC" "Pikal"], markershape=:circle)
# pl2 = blankplot_hrC(legend=:bottomright)
# @df dat2 exptfplot!(:t, :T1, :T2, :T3, markershape=:square)
# # plot!(fitdat2)
# modconvtplot!(sol2_MC128, sol2_Rp128, labels=["MC" "Pikal"], markershape=:square)
# plot!(xlim=(0, 120), ylim=(-36, -31), )

# pl8 = blankplot_hrC(legend=:bottomright)
# @df dat8 exptfplot!(:t, :T1, :T2, :T3, markershape=:circle)
# # plot!(fitdat8)
# modconvtplot!(sol8_MC128, sol8_Rp128, labels=["MC" "Pikal"], markershape=:circle)
# plot!(xlim=(0, 20), ylim=(-32, -27))
# plot(pl1, pl2, pl8, layout=(3,1), size=(600, 600), legend=:outerright)
# end
# savefig(plotsdir("suc128_T.svg"))

kstops18 = exp.(opt18mc.u[5:8]).*fit_mul_interp
kstops128 = exp.(opt128mc.u[5:8]).*fit_mul_interp
klogs18 = log.(ustrip.(u"s^-1", kstops18))
klogs128 = log.(ustrip.(u"s^-1", kstops128))
Tcheck = range(extrema(Tstops)..., length=50)
plot(Tstops, kstops18, yscale=:log10)
plot!(Tstops, kstops128, yscale=:log10)
plot!(Tstops, MC1_params.k_f, yscale=:log10)
plot!(Tstops, MC8_params.k_f, yscale=:log10)

# ------------------------
# MYEGA temparature dependence of viscosity
struct MYEGA_k
    k0
    K
    C
end
MYEGA_k(a) = MYEGA_k(a...)
(mk::MYEGA_k)(T) = mk.k0*exp((-mk.K/T)*exp(mk.C/T))
function k_myega(T, p)
    return MYEGA_k(p)(T)
end

myega_un = (u"s^-1", u"K", u"K")
myega_mult = (1e5u"s^-1", 800u"K", 600u"K")
p0_myega = [1.0, 1, 1]
function obj_kin_myega(p, _)
    if p[1] < 0 #|| p[2] < 100 #|| p[3] < 100 
        return NaN
    end
    mk = MYEGA_k(p.*myega_mult)
    ks = mk.(Tstops)
    logks = log.(ustrip.(u"s^-1", ks))
    # @info mk logks klogs18
    return sum(abs2.(logks .- klogs128))
end
objf_myega = OptimizationFunction(obj_kin_myega, AutoForwardDiff())
obj_myega(p0_myega, nothing)

opt_kin_myega = solve(OP(objf_myega, p0_myega, lb=fill(0.0, 3), ub=fill(Inf, 3)), alg)
myega_p = (opt_kin_myega.u .* myega_mult)

mk = MYEGA_k(myega_p)
mkt = MYEGA_k(myega_mult)
plot(Tcheck,  mk.(Tcheck))
plot!(Tcheck,  mkt.(Tcheck))
plot!(Tstops,  kstops18, yscale=:log10)
plot!(Tstops,  kstops128, yscale=:log10)

myega_params1 = @set MC1_params128.k_f = mk
myega_params2 = @set MC2_params128.k_f = mk
myega_params8 = @set MC8_params128.k_f = mk
sol1_myega = solve(ODEProblem(myega_params1), Rodas3())
sol2_myega = solve(ODEProblem(myega_params2, tspan=(0.0, 1000.0)), Rodas3())
sol8_myega = solve(ODEProblem(myega_params8), Rodas3())
begin
pl1 =modconvtplot(sol1_myega)
plot!(fitdat1)
pl2 = modconvtplot(sol2_myega)
plot!(fitdat2)
pl3 = modconvtplot(sol8_myega)
plot!(fitdat8)
plot(pl1, pl2, pl3)
end

fit_un_myega = (u"cm^-1", u"cm^5", u"m^0.5", u"m/s", u"s^-1", u"K", u"K")
fit_mult_myega = (myega_params1.N(0), myega_params1.S(0), myega_params1.C, myega_params1.Rs, myega_p...)
# function gen_sol_constNS_myega(fitprm, po; verbose=false)
#     new_params = deepcopy(po)
#     # fitprm = exp.(fitprmlog)
#     @reset new_params.N = RampedVariable(fitprm[1]*fit_mult_myega[1])
#     @reset new_params.S = RampedVariable(fitprm[2]*fit_mult_myega[2])
#     @reset new_params.C  = fitprm[3]*fit_mult_myega[3]
#     @reset new_params.Rs = fitprm[4]*fit_mult_myega[4]
#     @reset new_params.k_f = MYEGA_k(fitprm[5:7].*fit_mult_myega[5:7])
#     verbose && @info "fitprm" fitprm new_params
#     prob = ODEProblem(new_params, tspan=1000.0)
#     sol = solve(prob, Rodas3(autodiff=true))
#     return sol
# end
# function gennsol_constNS_myega(fitprm, pos)
#     sols = map(pos) do po
#         gen_sol_constNS_myega(fitprm, po)
#     end
#     return sols
# end
k_f_myega(vec) = MYEGA_k(exp.(vec).*mult_myega)
function objnmyega(fitprm, pfs; verbose=false)
    pos, fitdats = pfs
    # verbose && @info "fitprm" fitprm
    sols = gennsol_constNS(fitprm, pos, k_f_myega)
    if any([sol.retcode != ReturnCode.Terminated for sol in sols])
        if verbose
            for sol in sols 
                if sol.retcode != ReturnCode.Terminated
                    @info "unfinished drying" sol[1,end] fitprm 
                end
            end
        end
        # verbose && @info "unfinished drying" [sol for 
        # return NaN
    end
    objs = map(sols, fitdats) do sol, fitdat
        obj_expT(sol, fitdat, tweight=0)
    end
    verbose && @info "obj" objs fitprm 
    return sum(objs)
end
objnf_myega = OptimizationFunction(objnmyega, AutoForwardDiff())
p0_myegafit = fill(1.0, 7)
@time objnmyega(p0_myegafit, pass128, verbose=true)

pass128 = ([MC1_params128, MC2_params128, MC8_params128], [fitdat1, fitdat2, fitdat8]);
pass18 = ([MC1_params128, MC8_params128], [fitdat1, fitdat8]);

lb_myega = vcat(fill(0.0, 4), fill(-20., 3))
ub_myega = vcat(fill(Inf, 4), fill( 20., 3))

# Needs a maxiters, because the gradients have wildly different dimensions
@time opt18myega = solve(OP(objnf_myega, p0_myegafit, pass18, lb=lb_myega, ub=ub_myega), alg, g_tol=1e-2, maxiters=200)
# @time opt18myega = solve(OP(objnf_myega, p0_myegafit, pass18, lb=lb_myega, ub=ub_myega), alg, g_tol=1e-2, show_trace=true, maxiters=2)
# @time opt128myega = solve(OP(objnf_myega, p0_myegafit, pass128, lb=lb_myega, ub=ub_myega), alg, g_tol=1e-2, show_trace=true)

sol1_myega18  = gen_sol_constNS(opt18myega.u, MC1_params18, k_f_myega);
sol2_myega18  = gen_sol_constNS(opt18myega.u, MC2_params18, k_f_myega);
sol8_myega18  = gen_sol_constNS(opt18myega.u, MC8_params18, k_f_myega);
# sol1_myega128 = gen_sol_constNS_myega(opt128myega.u, MC1_params18);
# sol2_myega128 = gen_sol_constNS_myega(opt128myega.u, MC2_params18);
# sol8_myega128 = gen_sol_constNS_myega(opt128myega.u, MC8_params18);

begin
markers = [:square :circle :utriangle]
pl1 = blankplot_hrC(legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
plot!(title="Cycle 1, Tsh = $(Tsh1(Inf*u"s")|>u"°C")")
plot!(xlim=(0, 100), ylim=(-36, -31), )
modconvtplot!(sol1_myega, sol1_Rp18, labels=["MC" "Pikal"], markeralpha=0)

pl8 = blankplot_hrC(legend=:bottomright)
@df dat8 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
# plot!(fitdat8)
modconvtplot!(sol8_myega, sol8_Rp18, labels=["MC" "Pikal"], markeralpha=0)
plot!(xlim=(0, 45), ylim=(-32, -26))
plot!(title="Tsh = $(Tsh8(Inf*u"s")|>u"°C")")
plot(pl1, pl8, layout=(2,1), size=(600, 500), legend=:outerright)
end
savefig(plotsdir("suc18_myega_train_T.svg"))

# -------------------------------------
# playground temparature dependence of viscosity
struct Play_k
    k0
    kf
    K
    C
end
Play_k(a) = Play_k(a...)
function (pk::Play_k)(T)
    rat = pk.k0/pk.kf
    x = (T-pk.K)/pk.C
    # return k0/(tanh(-x)+1) # This has the right behavior...
    return k0*(exp(x) + exp(-x))/((1+rat)*exp(-x)+rat*exp(x)) # but this way it can avoid numerical issues
end
function k_play(T, p)
    return Play_k(p)(T)
end
play_un = (u"s^-1", u"s^-1", u"K", u"K")
play_mult = (1e-15u"s^-1", 1e-11u"s^-1", 236u"K", 0.8u"K")
p0_play = fill(1.0, 4)
# p0_play = [1e-17, 1e-10, 230, 0.8]
scatter(Tcheck,  k_play.(Tcheck, [p0_play .* play_mult]), yscale=:log10)
plot!(Tstops,  kstops18, yscale=:log10)
plot!(Tstops,  kstops128, yscale=:log10)

function obj_kin_play(p,_)
    if p[1] < 0 || p[2] < 0 #|| p[3] < 100 
        return NaN
    end
    ks = k_play.(Tstops, [p.*play_mult])
    # @info "check" ks kstops18
    sum(abs2.(log.(ustrip.(u"s^-1", ks)) - klogs18))
end
objf_play = OptimizationFunction(obj_play, AutoForwardDiff())
obj_kin_play(p0_play,nothing)
opt_play = solve(OP(objf_play, p0_play), alg, show_trace=true)
pk = Play_k(opt_play.u .* play_mult)
plot(Tcheck,  pk.(Tcheck))
plot!(Tstops,  kstops18, yscale=:log10)
plot!(Tstops,  kstops128, yscale=:log10)

play_params1 = @set MC1_params128.k_f = Play_k(opt_play.u.*play_mult)
play_params2 = @set MC2_params128.k_f = Play_k(opt_play.u.*play_mult)
play_params8 = @set MC8_params128.k_f = Play_k(opt_play.u.*play_mult)
sol1_play = solve(ODEProblem(play_params1), Rodas3())
sol2_play = solve(ODEProblem(play_params2), Rodas3())
sol8_play = solve(ODEProblem(play_params8), Rodas3())
poreevolplot(sol1_play)
poreevolplot!(sol2_play)
poreevolplot!(sol8_play)

begin
blankplot_hrC()
modconvtplot!(sol1_play, sol1_MC1)
plot!(fitdat1)
end
begin
blankplot_hrC()
modconvtplot!(sol8_play, sol8_MC8)
plot!(fitdat8)
end

k_f_play(vec) = Play_k(exp.(vec).*play_mult)
function objnplay(fitprm, pfs; verbose=false)
    pos, fitdats = pfs
    # verbose && @info "fitprm" fitprm
#     Tact = exp(fitprmlog[7] - 1)*fit_mult_play[7]
#     if Tact < 225u"K" || Tact > 250u"K"
#         verbose && @info "temperature unreasonable" exp(fitprmlog[7])
#         return NaN 
#     end
    sols = gennsol_constNS(fitprm, pos, k_f_play)
    if any([sol.retcode != ReturnCode.Terminated for sol in sols])
        if verbose
            for sol in sols 
                if sol.retcode != ReturnCode.Terminated
                    @info "unfinished drying" sol[1,end] fitprm 
                end
            end
        end
        # verbose && @info "unfinished drying" [sol for 
        # return NaN
    end
    objs = map(sols, fitdats) do sol, fitdat
        obj_expT(sol, fitdat, tweight=0)
    end
    verbose && @info "obj" objs fitprm 
    return sum(objs)
end
objnf_play = OptimizationFunction(objnplay, AutoForwardDiff())

lb_play = fill(0.0, 8)
lb_play[7] = log(225u"K"/play_mult[7]) + 1
ub_play = fill(Inf, 8)
ub_play[7] = log(250u"K"/play_mult[7]) + 1

opt18play = solve(OP(objnf_play, p0_playfitlog, pass18, lb=lb_play, ub=ub_play), Optim.NelderMead(), maxiters=50, time_limit=15)
# opt18play = solve(OP(objnf_play, p0_playfitlog, pass18), alg, g_tol = 1e-2, maxiters = 50, time_limit=15, show_trace=true)
# opt128play = solve(OP(objnf_play, p0_playfitlog, pass128), alg, g_tol = 1e-2, maxiters = 5, time_limit=15, show_trace=true)

sol1_play18, sol2_play18, sol8_play18 = gennsol_constNS(opt18play.u, (MC1_params128, MC2_params128, MC8_params128), k_f_play);
# (sol1_play128, sol2_play128, sol8_play128), playprms128 = gennsol_constNS_play(opt128play.u, MC1_params128, MC2_params128, MC8_params128);
# (sol1_playg, sol2_playg, sol8_playg), playprmsg = gennsol_constNS_play(p0_playfitlog, MC1_params2, MC2_params2, MC8_params128);
poreevolplot( sol1_play)
poreevolplot!(sol1_play18)
poreevolplot!(sol2_play18)
poreevolplot!(sol8_play18)

begin
pl1 = blankplot_hrC(legend=:bottomright)
markers = [:circle :utriangle :square]
@df dat1 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
plot!(title="Tsh = $(Tsh1(Inf*u"s")|>u"°C")")
plot!(xlim=(0, 120), ylim=(-36, -31), )
modconvtplot!(sol1_play18, sol1_Rp128, labels=["MC" "Pikal"], markeralpha=0)
pl2 = blankplot_hrC(legend=:bottomright)
@df dat2 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
modconvtplot!(sol2_play18, sol2_Rp128, labels=["MC" "Pikal"], markeralpha=0)
plot!(xlim=(0, 120), ylim=(-38, -31), )
plot!(title="Tsh = $(Tsh2(Inf*u"s")|>u"°C")")

pl8 = blankplot_hrC(legend=:bottomright)
@df dat8 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
modconvtplot!(sol8_play18, sol8_Rp128, labels=["MC" "Pikal"], markeralpha=0)
plot!(xlim=(0, 20), ylim=(-40, -20))
plot!(title="Tsh = $(Tsh8(Inf*u"s")|>u"°C")")
plot(pl2, pl1, pl8, layout=(3,1), size=(600, 600), legend=:outerright)
end
savefig(plotsdir("suc128_play_T.svg"))

begin
plot(u"°C", u"s^-1", legend=:bottomright, xlabel="T", ylabel="k\\  ", unitformat=latexsquareunitlabel)
# kineticplot!(MC1_params,    T=Tcheck, label="Interp fit to case 1")
kineticplot!(MC8_params,     label="Interp fit to case 8")
kineticplot!(MC1_params18,   c=1, marker=markers[2], label="Interp fit to 1&8")
kineticplot!(MC1_params128,  c=1, marker=markers[3], label="Interp fit to 1&2&8")
kineticplot!(myega_params1,   c=2, marker=markers[2], label="MYEGA fit to interpolation")
kineticplot!(play_params1  , c=3, marker=markers[2], label="Tanh fit to interpolation")
end
savefig(plotsdir("kinetic.svg"))



l = range(0u"cm", hf0, length=50)
Rp1 = Rp1_params.Rp.(l)
Rp2 = Rp2_params.Rp.(l)
Rp8 = Rp8_params.Rp.(l)
Rp18 = Rp1_params18.Rp.(l)

default(:markersize, 7)
begin
# colors = [:red :blue :purple :brown]
plot!(palette=:default)
rpattr = (seriestype=:samplemarkers, linestyle=:dash)
poreevolplot(sol1_play18, label="PBM 1&8, Cycle 1", c=1, marker=:circle)
poreevolplot!(sol2_play18, label="PBM 1&8, Cycle 2", c=2, marker=:square)
poreevolplot!(sol8_play18, label="PBM 1&8, Cycle 8", c=3, marker=:utriangle)
plot!(l, Rp1, subplot=4, label="3PP 1"; c=4, marker=:circle, rpattr...)
plot!(l, Rp2, subplot=4, label="3PP 2"; c=4, marker=:square, rpattr...)
plot!(l, Rp8, subplot=4, label="3PP 8"; c=4, marker=:utriangle, rpattr...)
# plot!(l, Rp18, subplot=4, label="3PP 1&8"; c=4, marker=:diamond, rpattr...)
plot!(size=(800, 400))
end
savefig(plotsdir("suc18_pores.svg"))

function namefitparam_arb(fitlog, label)
    fit = exp.(fitlog)
    nt = (label = label,
        # N =  fit[1]*fit_un_arb[1],
        # S =  fit[2]*fit_un_arb[2],
        # C =  fit[3]*fit_un_arb[3],
        # Rs = fit[4]*fit_un_arb[4],
        # k1 = fit[5]*fit_un_arb[5],
        # k2 = fit[6]*fit_un_arb[5],
        # k3 = fit[7]*fit_un_arb[5],
        # k4 = fit[8]*fit_un_arb[5],
        # T1 = Tstops[1],
        # T2 = Tstops[2],
        # T3 = Tstops[3],
        # T4 = Tstops[4],
    )
end
function namefitparam_play(fitlog, label)
    fit = exp.(fitlog)
    nt = (label = label,
        # N =  fit[1]*fit_un_play[1],
        # S =  fit[2]*fit_un_play[2],
        # C =  fit[3]*fit_un_play[3],
        # Rs = fit[4]*fit_un_play[4],
        # k0 = fit[5]*fit_un_play[5],
        # kf = fit[6]*fit_un_play[6],
        # Tact = fit[7]*fit_un_play[7],
        # θ  = fit[8]*fit_un_play[8],
    )
end
fits_arb = Optim.minimizer.([opt_MC1, opt_MC2, opt12mc, opt18mc, opt128mc]) 
fits_play = Optim.minimizer.([opt12play, opt18play, opt128play]) 
rep1 = Table(map(namefitparam_arb, fits_arb, ["1", "2", "1&2", "1&8", "1&2&8"]))
rep2 = Table(map(namefitparam_play, fits_play, ["1&2", "1&8", "1&2&8"]))

# safesave(datadir("exp_pro", "fit12.html"), pretty_table(HTML, report, header = ["Fit", "N", "S", "k0", "E", "C", "Rs"]))
safesave(datadir("exp_pro", "fit_arb.html"), pretty_table(HTML, rep1, header = ["Fit", "N", "S", "C", "Rs", "k1", "k2", "k3", "k4", "T1", "T2", "T3", "T4"]))
safesave(datadir("exp_pro", "fit_play.html"), pretty_table(HTML, rep2, header = ["Fit", "N", "S", "C", "Rs", "k0", "kf", "Tact", "θ"]))
# pretty_table(report, header = ["Fit", "N", "S", "k0", "E", "C", "Rs"])

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

fitdat3 = @df dat3 PrimaryDryFit(:t, (filter(!isnan, :T1), filter(!isnan, :T2), :T3[:t.<45u"hr"]))
obj3p = x->obj_expT(gen_sol_Rp(x, pop3), fitdat3)
opt3p = optimize(obj3p, opt_Rp1.u, NelderMead(), Optim.Options(g_tol = 1e-3, time_limit=15))
Rp3 = gen_sol_Rp(opt3p.u, pop3).prob.p.Rp

pop5 = deepcopy(Rp1_params18)
@reset pop5.pch = pch5
@reset pop5.Tsh = Tsh5
sol5_Rp18 = solve(ODEProblem(pop5), Rodas3());
po_PBM5 = deepcopy(playprms18[1])
@reset po_PBM5.pch = pch5
@reset po_PBM5.Tsh = Tsh5
sol5_PBM18 = solve(ODEProblem(po_PBM5), Rodas3());

fitdat5 = @df dat5 PrimaryDryFit(:t, (filter(!isnan, :T1), filter(!isnan, :T3), :T2[:t.<40u"hr"]))
plot(fitdat5)
obj5p = x->obj_expT(gen_sol_Rp(x, pop5), fitdat5)
opt5p = optimize(obj5p, opt_Rp1.u, NelderMead(), Optim.Options(g_tol = 1e-3, time_limit=15))
Rp5 = gen_sol_Rp(opt5p.u, pop5).prob.p.Rp

# Training plot
begin
pl1 = blankplot_hrC(legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
plot!(xlim=(0, 70), ylim=(-36, -31), )
modconvtplot!(sol1_play18, sol1_Rp18, labels=["PBM 1&8" "3PP 1&8"], markeralpha=0)
plot!(title="Tsh = $(Tsh1(Inf*u"s")|>u"°C")", xlabel="")
pl8 = blankplot_hrC(legend=:bottomright)
@df dat8 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
modconvtplot!(sol8_play18, sol8_Rp18, labels=["PBM 1&8" "3PP 1&8"], markeralpha=0)
plot!(xlim=(0, 30), ylim=(-32, -26))
plot!(title="Tsh = $(Tsh8(Inf*u"s")|>u"°C")", )
plot(pl1, pl8, layout=(2,1), size=(600, 400), legend=:outerright)
end
savefig(plotsdir("suc18_train_T.svg"))

# Validation plot

begin
pl2 = blankplot_hrC()
@df dat2 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
modconvtplot!(sol2_play18, sol2_Rp18, labels=["PBM 1&8" "3PP 1&8" ], markeralpha=0)
plot!(title="Tsh = $(Tsh2(Inf*u"s")|>u"°C")", xlabel="", ylabel="")
pl3 = blankplot_hrC()
@df dat3 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
modconvtplot!(sol3_PBM18, sol3_Rp18, labels=["PBM 1&8" "3PP 1&8" ], markeralpha=0)
plot!(xlim=(0, 80))
plot!(title="Tsh = $(Tsh3(Inf*u"s")|>u"°C")", xlabel="")
pl5 = blankplot_hrC()
@df dat5 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
modconvtplot!(sol5_PBM18, sol5_Rp18, labels=[ "PBM 1&8" "3PP 1&8"], markeralpha=0)
plot!(xlim=(0, 55))
plot!(title="Tsh = $(Tsh5(Inf*u"s")|>u"°C")", ylabel="")
plot(pl2, pl3, pl5, layout=(3,1))
plot!(ylim=(-37, -32), legend=:outerright)
end
savefig(plotsdir("suc18train_235val_T.svg"))

# -----------------
# Inspect fit results in Rp space rather than temperature space

begin
pl1 = blankplot_lRp(ylims=(0, 12), legend=:topleft)
plot!(title="Cycle 1, Tsh = $(Tsh1(Inf*u"s")|>u"°C")")
pbmtdlas!(sol1_play18, dat1, Rp1_params.Rp, labels=["PBM 1&8" "3PP 1" "TDLAS"])
plot!(xlabel="")
pl2 = blankplot_lRp(ylims=(0, 12), legend=:topleft,)
plot!(title="Cycle 2, Tsh = $(Tsh2(Inf*u"s")|>u"°C")")
pbmtdlas!(sol2_play18, dat2, Rp2_params.Rp, labels=["PBM 1&8" "3PP 2" "TDLAS"])
pl5 = blankplot_lRp(ylims=(0, 12), legend=:topleft,)
plot!(title="Cycle 5, Tsh = $(Tsh5(Inf*u"s")|>u"°C")")
pbmtdlas!(sol5_PBM18, filter(x->x.Rp>0u"m/s", dat5), Rp5, labels=["PBM 1&8" "3PP 5" "TDLAS"])
plot!(ylabel="")
pl8 = blankplot_lRp(ylims=(0, 12), legend=:topleft,)
plot!(title="Cycle 8, Tsh = $(Tsh8(Inf*u"s")|>u"°C")")
pbmtdlas!(sol8_play18, dat8, Rp8_params.Rp, labels=["PBM 1&8" "3PP 8" "TDLAS"])
plot!(xlabel="", ylabel="")
plot(pl1, pl8, pl2, pl5, layout=(2,2), size=(600,500))
end
savefig(plotsdir("suc18_Rp_val.svg"))

# -----------------
# Inspect with MYEGA fit instead of tanh, since seems better-conditioned


sol2_myega18, myega2_params18 = gen_sol_constNS_myega(opt18myega.u, MC2_params18);
sol3_myega18, myega3_params18 = gen_sol_constNS_myega(opt18myega.u, po_PBM3);
sol5_myega18, myega5_params18 = gen_sol_constNS_myega(opt18myega.u, po_PBM5);
begin
pl2 = blankplot_hrC(legend=:bottomright)
@df dat2 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
modconvtplot!(sol2_myega18, sol2_Rp18, labels=["MC" "Pikal"], markeralpha=0)
plot!(xlim=(0, 120), ylim=(-38, -31), )
plot!(title="Cycle 2, Tsh = $(Tsh2(Inf*u"s")|>u"°C")")
pl3 = blankplot_hrC(legend=:bottomright)
@df dat3 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
modconvtplot!(sol3_myega18, sol3_Rp18, labels=["MC" "Pikal"], markeralpha=0)
plot!(xlim=(0, 70), ylim=(-38, -31), )
plot!(title="Cycle 3, Tsh = $(Tsh3(Inf*u"s")|>u"°C")")
pl5 = blankplot_hrC(legend=:bottomright)
@df dat5 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
modconvtplot!(sol5_myega18, sol5_Rp18, labels=["MC" "Pikal"], markeralpha=0)
plot!(xlim=(0, 51), ylim=(-38, -31), )
plot!(title="Cycle 5, Tsh = $(Tsh5(Inf*u"s")|>u"°C")")
plot(pl2, pl3, pl5, layout=(3,1), size=(600, 600), legend=:outerright)
end
savefig(plotsdir("suc18_myega_val235_T.svg"))
