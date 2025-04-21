using DrWatson
@quickactivate :FulbrightPores

const OP = OptimizationProblem
const NLSP = NonlinearLeastSquaresProblem
using LineSearches
alg = Optim.BFGS(linesearch=LineSearches.BackTracking())
pythonplot()

# ------------------
# Load Kv function
# Computed in scripts/kv.jl, because it is used elsewhere too
@unpack K_shf = load(datadir("exp_pro", "simtra_kv.jld2"))

# --------------------------
# Cycle 1
# Labeled as case 4 in Table 8
dat1raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle1.csv"), Table)
dat1 = map(simtra_rename, dat1raw)

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

fitdat1 = PrimaryDryFit(dat1.t, (
    filter(x->(!isnan(x) && x<-32.5u"°C"), dat1.T1),
    filter(!isnan, dat1.T2[dat1.t .< 60u"hr"]),
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

fitdat2 = @df dat2 PrimaryDryFit(:t, (
    # collect(filter(!isnan, :T1)),
    filter(!isnan, :T2[:t .< 85u"hr"]),
    filter(!isnan, :T3),),
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
fitdat3 = @df dat3 PrimaryDryFit(:t, (filter(!isnan, :T1), filter(!isnan, :T2), :T3[:t.<45u"hr"]))

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

fitdat5 = @df dat5 PrimaryDryFit(:t, (filter(!isnan, :T1), filter(!isnan, :T3), :T2[:t.<40u"hr"]))


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

fitdat8 = @df dat8 PrimaryDryFit(:t, (
    filter(!isnan, :T1),
    filter(!isnan, :T2[:t .< 17u"hr"]),
    filter(!isnan, :T3[:t .< 17u"hr"]),),
    t_end8)
plot(fitdat8)

@info "All data loaded"

# -------------
# Parameter setup

rad_i = 1.1u"cm"
rad_o = 1.2u"cm"
A_p = π*rad_i^2
A_v = π*rad_o^2
fillvol = 5u"mL"
hf0 = fillvol/A_p

Rp_g = RpFormFit(0.8u"cm^2*Torr*hr/g", 14u"cm*Torr*hr/g", 1u"1/cm")

csolid = 0.04u"g/mL"
ρsolution = 1u"g/mL"

pb_pikal = ParamObjPikal([
    (Rp_g, hf0, csolid, ρsolution),
    (K_shf, A_v, A_p),
    (pch1, Tsh1),
])

function new_Tp(po, Tsh, pch)
    po_n = deepcopy(po)
    @reset po_n.pch = pch
    @reset po_n.Tsh = Tsh
    return po_n
end
pop1 = pb_pikal
pop2 = new_Tp(pop1, Tsh2, pch2)
pop8 = new_Tp(pop1, Tsh8, pch8)
pop5 = new_Tp(pop1, Tsh5, pch5)
pop3 = new_Tp(pop1, Tsh3, pch3)

obj_Rpf = OptimizationFunction(obj_Rp, AutoForwardDiff())
@time obj_Rpf(fill(0.0, 3), (pop8, fitdat8))

obj_npf = OptimizationFunction(objn_Rp, AutoForwardDiff())
# @time ForwardDiff.gradient(x->obj_Rpf(x, (pop8, fitdat8)), [0.8, 90, 1])

@time opt1p = solve(OP(obj_Rpf, fill(1.0, 3), (pop1, fitdat1)),alg, g_tol = 1e-3)
opt2p = solve(OP(obj_Rpf, opt1p.u, (pop2, fitdat2)), alg; g_tol=1e-3)
opt8p = solve(OP(obj_Rpf, opt1p.u, (pop8, fitdat8)), alg; g_tol=1e-3)

sol1_Rp1 = gen_sol_Rp(opt1p.u, pop1);
sol2_Rp2 = gen_sol_Rp(opt2p.u, pop2);
sol8_Rp8 = gen_sol_Rp(opt8p.u, pop8);
# plot(u"hr", u"°C", xlim=(0, 60), ylim=(-40, -30))
blankplot_hrC(ylim=(-40, -30))
modconvtplot!(sol1_Rp1)
# @df dat1 exptfplot!(:t, :T1, :T2, :T3)
plot!(fitdat1)

pass18p =  ((pop1, pop2), (fitdat1, fitdat2));
pass128p = ((pop1, pop2, pop8), (fitdat1, fitdat2, fitdat8))
@time opt18p = solve(OP(obj_npf, opt1p.u, pass18p), alg; g_tol=1e-2, show_trace=true)
@time opt128p = solve(OP(obj_npf, opt1p.u, pass128p), alg; g_tol=1e-2)

sol1_Rp18 = gen_sol_Rp(opt18p.u, pop1);
sol2_Rp18 = gen_sol_Rp(opt18p.u, pop2);
sol8_Rp18 = gen_sol_Rp(opt18p.u, pop8);
sol1_Rp128 = gen_sol_Rp(opt128p.u, pop1);
sol2_Rp128 = gen_sol_Rp(opt128p.u, pop2);
sol8_Rp128 = gen_sol_Rp(opt128p.u, pop8);
Rp1_params18 = sol1_Rp18.prob.p
Rp1_params128 = sol1_Rp128.prob.p
Rp1_params = sol1_Rp1.prob.p
Rp2_params = sol2_Rp2.prob.p
Rp8_params = sol8_Rp8.prob.p

opt3p = solve(OP(obj_Rpf, opt1p.u, (pop3, fitdat3)), alg, g_tol = 1e-3)
opt5p = solve(OP(obj_Rpf, opt1p.u, (pop5, fitdat5)), alg, g_tol = 1e-3)
# -------------- Set up parameters

# typical diameter: 50 μm
d_typical = 50u"μm"
V_typical = 4/3*π*d_typical^3
N_typical = 1/V_typical
# Good guesses for Knudsen
N = ConstPhysProp(N_typical)
S = ConstPhysProp(V_typical) # = V_typ^2*N_typ
C = 2e-1u"m^0.5"
Rs = 0.5u"cm^2*Torr*hr/g"
k0 = 1e-14u"1/s"

pb_mc = ParamObjMC([
    (hf0, csolid, ρsolution),
    (K_shf, A_v, A_p),
    (pch1, Tsh1),
    (ConstPhysProp(k0), C, N, S, Rs)
])
pomc1 = pb_mc
pomc2 = new_Tp(pomc1, Tsh2, pch2)
pomc8 = new_Tp(pomc1, Tsh8, pch8)

# -----------------------------
# Fit k by interpolating at four values
function k_f_interp(logk)
    return InterpK(logk, [-45.0, -38, -32, -28]u"°C" .|> u"K", k0)
end
# function obj_interp(x,p)
#     ks = x[5:8]
#     if ~(ks[1] <= ks[2] <= ks[3] <= ks[4])
#         return NaN
#     end
#     return objn_mc(x, (p..., k_f_interp))
# end
function obj_ninterp(fitprm, pf; verbose=false)
    if ~(fitprm[5] < fitprm[6] < fitprm[7] < fitprm[8])
        return NaN
    end
    return objn_mc(fitprm, (pf..., k_f_interp), verbose=verbose)
end
# objf_interp = OptimizationFunction(obj_interp, AutoForwardDiff())
objnf_interp = OptimizationFunction(obj_ninterp, AutoForwardDiff())
lb_MC = vcat(fill(0.0, 4), fill(-20., 4))
ub_MC = vcat(fill(Inf, 4), fill(20., 4))

interp_guess = [1.0, 1.0, 1.0, 1.0, -1, 1, 3, 12]
@time opt_MC1 = solve(OP(objnf_interp, interp_guess, (pomc1, fitdat1), lb=lb_MC, ub=ub_MC), alg;show_trace=true, g_tol = 1e-2, maxtime=10)

sol1_MC1 = gen_sol_mc_arb(opt_MC1.u, pomc1, k_f_interp);
MC1_params1 = sol1_MC1.prob.p

begin
plot(u"hr", u"°C", xlim=(0, 120), ylim=(-36, -31), legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3, marker=:square)
modconvtplot!(sol1_MC1, sol1_Rp1, labels=["MC" "Pikal"],marker=:square)
end
savefig(plotsdir("suc1_T.svg"))


@time opt_MC8 = solve(OP(objnf_interp, opt_MC1.u, (pomc8, fitdat8), lb=lb_MC, ub=ub_MC), alg; show_trace=true, g_tol=1e-2)

sol8_MC8 = gen_sol_mc_arb(opt_MC8.u, pomc8, k_f_interp);
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
@time opt18mc = solve(OP(objnf_interp, g18mc, ((pomc1, pomc8), (fitdat1, fitdat8)), lb=lb_MC, ub=ub_MC), alg, g_tol=1e-3)
# opt128mc = optimize(obj128mc, g128mc, NelderMead(), Optim.Options(g_tol = 1e-3))
@time opt128mc = solve(OP(objnf_interp, g18mc, ((pomc1, pomc2, pomc8), (fitdat1, fitdat2, fitdat8)), lb=lb_MC, ub=ub_MC), alg, g_tol=1e-3)

sol1_MC18  = gen_sol_mc_arb(opt18mc.u, pomc1, k_f_interp);
sol2_MC18  = gen_sol_mc_arb(opt18mc.u, pomc2, k_f_interp);
sol8_MC18  = gen_sol_mc_arb(opt18mc.u, pomc8, k_f_interp);
MC1_params18 = sol1_MC18.prob.p
MC2_params18 = sol2_MC18.prob.p
MC8_params18 = sol8_MC18.prob.p

sol1_MC128  = gen_sol_mc_arb(opt128mc.u, pomc1, k_f_interp);
sol2_MC128  = gen_sol_mc_arb(opt128mc.u, pomc2, k_f_interp);
sol8_MC128  = gen_sol_mc_arb(opt128mc.u, pomc8, k_f_interp);
MC1_params128 = sol1_MC128.prob.p
MC2_params128 = sol2_MC128.prob.p
MC8_params128 = sol8_MC128.prob.p

begin
pl1 = blankplot_hrC(legend=:bottomright)
# @df dat1 exptfplot!(:t, :T1, :T2, :T3, markershape=:circle)
plot!(fitdat1)
plot!(xlim=(0, 100), ylim=(-36, -31), )
modconvtplot!(sol1_MC18, sol1_Rp18, labels=["MC" "Pikal"], markeralpha=0)
pl8 = blankplot_hrC(legend=:bottomright)
# @df dat8 exptfplot!(:t, :T1, :T2, :T3, markershape=:square)
plot!(fitdat8, markershape=:square)
modconvtplot!(sol8_MC18, sol8_Rp18, labels=["MC" "Pikal"], markeralpha=0)
plot!(xlim=(0, 45), ylim=(-32, -26), )

plot(pl1, pl8, layout=(2,1), legend=:outerright)
end
savefig(plotsdir("suc18_T.svg"))

begin
pl1 = blankplot_hrC(legend=:bottomright)
# @df dat1 exptfplot!(:t, :T1, :T2, :T3, markershape=:circle)
plot!(fitdat1)
plot!(xlim=(0, 60), ylim=(-36, -31), )
modconvtplot!(sol1_MC128, sol1_Rp128, labels=["MC" "Pikal"], markeralpha=0)
pl2 = blankplot_hrC(legend=:bottomright)
# @df dat2 exptfplot!(:t, :T1, :T2, :T3, markershape=:square)
plot!(fitdat2)
modconvtplot!(sol2_MC128, sol2_Rp128, labels=["MC" "Pikal"], markeralpha=0)
plot!(xlim=(0, 120), ylim=(-36, -31), )

pl8 = blankplot_hrC(legend=:bottomright)
# @df dat8 exptfplot!(:t, :T1, :T2, :T3, markershape=:circle)
plot!(fitdat8)
modconvtplot!(sol8_MC128, sol8_Rp128, labels=["MC" "Pikal"], markeralpha=0)
plot!(xlim=(0, 20), ylim=(-32, -27))
plot(pl1, pl2, pl8, layout=(3,1), size=(600, 600), legend=:outerright)
end
savefig(plotsdir("suc128_T.svg"))

Tstops = MC1_params18.k_f.spl.t
kstops18 = get_ks(MC1_params18.k_f)
kstops128 = get_ks(MC1_params128.k_f)
klogs18 = log.(ustrip.(u"s^-1", kstops18))
klogs128 = log.(ustrip.(u"s^-1", kstops128))
Tcheck = range(extrema(Tstops)..., length=50)
plot(Tstops, kstops18, yscale=:log10)
plot!(Tstops, kstops128, yscale=:log10)
plot!(Tstops, MC1_params1.k_f.(Tstops), yscale=:log10)
plot!(Tstops, MC8_params.k_f.(Tstops), yscale=:log10)

# ------------------------
# MYEGA temparature dependence of viscosity

myega_scale = (1e10u"s^-1", 1.6e3u"K", 1e2u"K")
k_f_myega(vec) = MYEGA_k((exp.(vec).*myega_scale))
p0_myega = fill(0.1, 3)
function obj_kin_myega(p, _)
    # if p[1] < 0 #|| p[2] < 100 #|| p[3] < 100 
    #     return NaN
    # end
    mk = k_f_myega(p)
    ks = mk.(Tstops)
    logks = log.(ustrip.(u"s^-1", ks))
    return sum(abs2, logks .- klogs128)
end
objf_myega = OptimizationFunction(obj_kin_myega, AutoForwardDiff())
objf_myega(p0_myega, nothing)

# opt_kin_myega = solve(OP(objf_myega, p0_myega, lb=fill(0.0, 3), ub=fill(Inf, 3)), alg)
opt_kin_myega = solve(OP(objf_myega, p0_myega), alg)

mk = k_f_myega(opt_kin_myega.u)
mk2 = MYEGA_k(1.25e13u"s^-1", 1.6e3u"K", 5.1e2u"K")
plot(Tcheck,  mk.(Tcheck))
plot!(Tcheck,  mk2.(Tcheck))
scatter!(Tstops,  kstops18, yscale=:log10)
scatter!(Tstops,  kstops128, yscale=:log10)

myega_params1 = @set MC1_params128.k_f = mk2
myega_params2 = @set MC2_params128.k_f = mk2
myega_params8 = @set MC8_params128.k_f = mk2
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

objn_myega(fitprm, pf) = objn_mc(fitprm, (pf..., k_f_myega))
objnf_myega = OptimizationFunction(objn_myega, AutoForwardDiff())
p0_myegafit = vcat(fill(1.0, 4), fill(1.0, 3))

pass128 = ([MC1_params128, MC2_params128, MC8_params128], [fitdat1, fitdat2, fitdat8]);
pass18 = ([MC1_params128, MC8_params128], [fitdat1, fitdat8]);
@time objnf_myega(p0_myegafit, pass128)

lb_myega = vcat(fill(0.0, 4), fill(-100., 3))
ub_myega = vcat(fill(30, 4), fill( 100., 3))

# Needs a maxiters, because the gradients have wildly different dimensions
@time opt18myega = solve(OP(objnf_myega, p0_myegafit, pass18, lb=lb_myega, ub=ub_myega), alg, g_tol=1e-2, maxiters=200)
@time opt128myega = solve(OP(objnf_myega, p0_myegafit, pass128, lb=lb_myega, ub=ub_myega), alg, g_tol=1e-2, maxiters=200)

sol1_myega18  = gen_sol_mc_arb(opt18myega.u, MC1_params18, k_f_myega);
sol2_myega18  = gen_sol_mc_arb(opt18myega.u, MC2_params18, k_f_myega);
sol8_myega18  = gen_sol_mc_arb(opt18myega.u, MC8_params18, k_f_myega);
sol1_myega128 = gen_sol_mc_arb(opt128myega.u, MC1_params18, k_f_myega);
sol2_myega128 = gen_sol_mc_arb(opt128myega.u, MC2_params18, k_f_myega);
sol8_myega128 = gen_sol_mc_arb(opt128myega.u, MC8_params18, k_f_myega);

objn_mc(opt128myega.u, (pass128..., k_f_myega), verbose=true)
begin
markers = [:square :circle :utriangle]
pl1 = blankplot_hrC(legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
plot!(title="Cycle 1, Tsh = $(Tsh1(Inf*u"s")|>u"°C")")
plot!(xlim=(0, 100), ylim=(-36, -31), )
modconvtplot!(sol1_myega128, sol1_Rp128, labels=["MC" "Pikal"], markeralpha=0)
pl2 = blankplot_hrC(legend=:bottomright)
@df dat2 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
plot!(title="Cycle 2, Tsh = $(Tsh2(Inf*u"s")|>u"°C")")
plot!(xlim=(0, 100), ylim=(-36, -31), )
modconvtplot!(sol2_myega128, sol2_Rp128, labels=["MC" "Pikal"], markeralpha=0)

pl8 = blankplot_hrC(legend=:bottomright)
@df dat8 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
# plot!(fitdat8)
modconvtplot!(sol8_myega128, sol8_Rp128, labels=["MC" "Pikal"], markeralpha=0)
plot!(xlim=(0, 45), ylim=(-32, -26))
plot!(title="Tsh = $(Tsh8(Inf*u"s")|>u"°C")")
plot(pl1, pl2, pl8, layout=(3,1), size=(800, 500), legend=:outerright)
end
savefig(plotsdir("suc128_myega_train_T.svg"))

# -------------------------------------
# sigmoid temparature dependence of viscosity
const sigm_scale = (1e-15u"s^-1", 1e-9u"s^-1", 0.8u"K", 239u"K")
function k_f_sigm(vec)
    p123 = exp.(vec[1:3] .-1).*sigm_scale[1:3]
    p4 = vec[4]*u"K" + sigm_scale[4]
    Sigmoid_k(p123..., p4)
end
p0_sigm = fill(1.0, 4)
plot(Tcheck,  k_f_sigm(p0_sigm).(Tcheck), yscale=:log10)
scatter!(Tstops,  kstops18, yscale=:log10)
scatter!(Tstops,  kstops128, yscale=:log10)

function err_kin_sigm(p,_)
    ks = k_f_sigm(p).(Tstops)
    err = log.(ustrip.(u"s^-1", ks)) - klogs128
end
# objf_sigm = OptimizationFunction(obj_kin_sigm, AutoForwardDiff())
# obj_kin_sigm(p0_sigm,nothing)
opt_sigm = solve(NLSP(err_kin_sigm, p0_sigm), LevenbergMarquardt(), show_trace=Val(true))
pk = k_f_sigm(opt_sigm.u)
pkt = k_f_sigm(p0_sigm)
plot(Tcheck,  pk.(Tcheck))
plot!(Tcheck,  pkt.(Tcheck))
scatter!(Tstops,  kstops18, yscale=:log10)
scatter!(Tstops,  kstops128, yscale=:log10)

sigm_params1 = @set MC1_params128.k_f = k_f_sigm(opt_sigm.u)
sigm_params2 = @set MC2_params128.k_f = k_f_sigm(opt_sigm.u)
sigm_params8 = @set MC8_params128.k_f = k_f_sigm(opt_sigm.u)
sol1_sigm = solve(ODEProblem(sigm_params1), Rodas3())
sol2_sigm = solve(ODEProblem(sigm_params2), Rodas3())
sol8_sigm = solve(ODEProblem(sigm_params8), Rodas3())
poreevolplot(sol1_sigm)
poreevolplot!(sol2_sigm)
poreevolplot!(sol8_sigm)

begin
blankplot_hrC()
modconvtplot!(sol1_sigm, sol1_MC1)
plot!(fitdat1)
end
begin
blankplot_hrC()
modconvtplot!(sol8_sigm, sol8_MC8)
plot!(fitdat8)
end

objn_sigm(fitprm, pf) = objn_mc(fitprm, (pf..., k_f_sigm))
objnf_sigm = OptimizationFunction(objn_sigm, AutoForwardDiff())

lb_sigm = vcat(fill(0.0, 4), fill(-10, 4))
lb_sigm[8] = -30
ub_sigm = vcat(fill(Inf, 4), fill( 10, 4))
ub_sigm[8] = 30

p0_sigmfit = vcat(fill(1.0, 4), fill(1.0, 4))

opt18sigm = solve(OP(objnf_sigm, p0_sigmfit, pass18, lb=lb_sigm, ub=ub_sigm), alg, maxiters=50, time_limit=15)
# opt18sigm = solve(OP(objnf_sigm, p0_sigmfitlog, pass18), alg, g_tol = 1e-2, maxiters = 50, time_limit=15, show_trace=true)
opt128sigm = solve(OP(objnf_sigm, p0_sigmfit, pass128, lb=lb_sigm, ub=ub_sigm), alg, maxiters = 50, time_limit=15)

sol1_sigm18, sol2_sigm18, sol8_sigm18 = gennsol_mc(opt18sigm.u, (MC1_params128, MC2_params128, MC8_params128), k_f_sigm);
sol1_sigm128, sol2_sigm128, sol8_sigm128 = gennsol_mc(opt128sigm.u, (MC1_params128, MC2_params128, MC8_params128), k_f_sigm);
# poreevolplot( sol1_sigm)
poreevolplot( sol1_sigm18)
poreevolplot!(sol2_sigm18)
poreevolplot!(sol8_sigm18)

begin
pl1 = blankplot_hrC(legend=:bottomright)
markers = [:circle :utriangle :square]
@df dat1 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
plot!(title="Tsh = $(Tsh1(Inf*u"s")|>u"°C")")
plot!(xlim=(0, 120), ylim=(-36, -31), )
modconvtplot!(sol1_sigm128, sol1_Rp128, labels=["MC" "Pikal"], markeralpha=0)
pl2 = blankplot_hrC(legend=:bottomright)
@df dat2 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
modconvtplot!(sol2_sigm128, sol2_Rp128, labels=["MC" "Pikal"], markeralpha=0)
plot!(xlim=(0, 120), ylim=(-38, -31), )
plot!(title="Tsh = $(Tsh2(Inf*u"s")|>u"°C")")

pl8 = blankplot_hrC(legend=:bottomright)
@df dat8 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
modconvtplot!(sol8_sigm128, sol8_Rp128, labels=["MC" "Pikal"], markeralpha=0)
plot!(xlim=(0, 20), ylim=(-40, -20))
plot!(title="Tsh = $(Tsh8(Inf*u"s")|>u"°C")")
plot(pl2, pl1, pl8, layout=(3,1), size=(600, 600), legend=:outerright)
end
savefig(plotsdir("suc128_sigm_T.svg"))

mc_fits = Dict{Tuple,ParamObjMC}(
    (InterpK, 1, (1,)) => MC1_params1,
    (InterpK, 8, (8,)) => MC8_params,
    (InterpK, 1, (1, 8)) => MC1_params18,
    (InterpK, 8, (1, 8)) => MC8_params18,
    (InterpK, 1, (1, 2, 8)) => MC1_params128,
    (InterpK, 2, (1, 2, 8)) => MC2_params128,
    (InterpK, 8, (1, 2, 8)) => MC8_params128,
    (MYEGA_k, 1, (1, 8)) => sol1_myega18.prob.p,
    (MYEGA_k, 1, (1, 2, 8)) => sol1_myega128.prob.p,
    (MYEGA_k, 2, (1, 2, 8)) => sol2_myega128.prob.p,
    (MYEGA_k, 8, (1, 2, 8)) => sol8_myega128.prob.p,
    (Sigmoid_k, 1, (1, 8)) => sol1_sigm18.prob.p,
    (Sigmoid_k, 1, (1, 2, 8)) => sol1_sigm128.prob.p,
    (Sigmoid_k, 2, (1, 2, 8)) => sol2_sigm128.prob.p,
    (Sigmoid_k, 8, (1, 2, 8)) => sol8_sigm128.prob.p,
)

pikal_params = Dict(
    (1, (1,)) => Rp1_params,
    (2, (2,)) => Rp2_params,
    (8, (8,)) => Rp8_params,
    (1, (1,8)) => Rp1_params18,
    (1, (1,2,8)) => Rp1_params128,
)

Rp1 = Rp1_params.Rp
Rp2 = Rp2_params.Rp
Rp8 = Rp8_params.Rp
Rp18 = Rp1_params18.Rp
Rp128 = Rp1_params128.Rp
Rp3 = gen_sol_Rp(opt3p.u, pop3).prob.p.Rp
Rp5 = gen_sol_Rp(opt5p.u, pop5).prob.p.Rp

pikal_rps = Dict(
    (1,) => Rp1,
    (2,) => Rp2,
    (8,) => Rp8,
    (1,8) => Rp18,
    (1,2,8) => Rp128,
    (3,) => Rp3,
    (5,) => Rp5,
)

fitdats = Dict(
    1 => fitdat1,
    2 => fitdat2,
    3 => fitdat3,
    5 => fitdat5,
    8 => fitdat8
)

safesave(datadir("exp_pro", "simtra_sucrose_fitdata.jld2"), @strdict fitdats)
safesave(datadir("exp_pro", "simtra_sucrose_fitresults.jld2"), @strdict mc_fits pikal_params pikal_rps)

begin
    gr()
plot(u"°C", u"s^-1", legend=:bottomright, xlabel="T", ylabel="k\\  ", unitformat=latexsquareunitlabel)
# kineticplot!(MC1_params,  c  label="Interp fit to case 1")
# kineticplot!(MC8_params,     label="Interp fit to case 8")
kineticplot!(MC1_params18,   c=1, marker=markers[2], label="Interp fit to 1&8")
kineticplot!(MC1_params128,  c=1, marker=markers[3], label="Interp fit to 1&2&8")
kineticplot!(myega_params1,   c=2, marker=markers[1], label="MYEGA fit to interpolation")
kineticplot!(sol1_myega18,   c=2, marker=markers[2], label="MYEGA fit to 1&8")
kineticplot!(sol1_myega128,   c=2, marker=markers[3], label="MYEGA fit to 1&2&8")
kineticplot!(sigm_params1  , c=3, marker=markers[1], label="Tanh fit to interpolation")
kineticplot!(sol1_sigm18  , c=3, marker=markers[2], label="Tanh fit to 1&8")
kineticplot!(sol1_sigm128  , c=3, marker=markers[3], label="Tanh fit to 1&2&8")
end
savefig(plotsdir("suc_kinetic.svg"))
pythonplot()



pythonplot()
l = range(0u"cm", hf0, length=50)
default(:markersize, 7)
begin
# colors = [:red :blue :purple :brown]
plot!(palette=:default)
rpattr = (seriestype=:samplemarkers, linestyle=:dash)
poreevolplot( sol1_sigm128, label="PBM 1&2&8, Cycle 1", c=1, marker=:circle)
poreevolplot!(sol2_sigm128, label="PBM 1&2&8, Cycle 2", c=2, marker=:square)
poreevolplot!(sol8_sigm128, label="PBM 1&2&8, Cycle 8", c=3, marker=:utriangle)
plot!(l, Rp1.(l), subplot=4, label="3PP 1"; c=4, marker=:circle, rpattr...)
plot!(l, Rp2.(l), subplot=4, label="3PP 2"; c=4, marker=:square, rpattr...)
plot!(l, Rp8.(l), subplot=4, label="3PP 8"; c=4, marker=:utriangle, rpattr...)
# plot!(l, Rp18, subplot=4, label="3PP 1&8"; c=4, marker=:diamond, rpattr...)
plot!(size=(800, 400))
end
savefig(plotsdir("suc128_pores.svg"))

# # -------------------- 
# # Assume a beta distribution and see results
# gr()
# betaplot(sol1_sigm128, legend=:topright, legend_columns=1)
# # betaplot(sol8_sigm128, ylim=(1e-10, 1e35), xscale=:log10, yscale=:log10, legend=:bottomleft, legend_columns=1)
# plot!(xlabel="Pore volume [mm³]", ylabel="Probability Density")
# savefig(plotsdir("changingdist_linear.png"))

# -----------------------

function namefitparam_arb(sol, label)
    po = sol.prob.p
    nt = (label = label,
        N = po.N,
        S = po.S,
        C = po.C,
        Rs = po.Rs,
        k_f = po.k_f,
    )
end
fits_sol = [sol1_MC1, sol8_MC8, sol1_MC18, sol1_MC128, sol1_myega128, sol1_sigm128]
fits_lab = ["1", "8", "1&8", "1&2&8", "MYEGA 1&2&8", "Tanh 1&2&8"]
# fits_sigm = Optim.minimizer.([opt12sigm, opt18sigm, opt128sigm]) 
rep = Table(map(namefitparam_arb, fits_sol, fits_lab))

# safesave(datadir("exp_pro", "fit12.html"), pretty_table(HTML, report, header = ["Fit", "N", "S", "k0", "E", "C", "Rs"]))
safesave(datadir("exp_pro", "fit_arb.html"), pretty_table(HTML, rep, header = ["Fit", "N", "S", "C", "Rs", "kf"]))
# pretty_table(report, header = ["Fit", "N", "S", "k0", "E", "C", "Rs"])

# ---------------------
# Validation


sol3_Rp128 = solve(ODEProblem(pop3), Rodas3());
po_PBM3 = new_Tp(MC1_params128, Tsh3, pch3)
sol3_PBM128 = solve(ODEProblem(po_PBM3), Rodas3());


sol5_Rp128 = solve(ODEProblem(pop5), Rodas3());
po_PBM5 = new_Tp(MC1_params128, Tsh5, pch5)
sol5_PBM128 = solve(ODEProblem(po_PBM5), Rodas3());

plot(fitdat5)

# Training plot
begin
pl1 = blankplot_hrC(legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
plot!(xlim=(0, 70), ylim=(-36, -31), )
modconvtplot!(sol1_sigm128, sol1_Rp128, labels=["PBM 1&8" "3PP 1&2&8"], markeralpha=0)
plot!(title="Tsh = $(Tsh1(Inf*u"s")|>u"°C")", xlabel="")
pl2 = blankplot_hrC(legend=:bottomright)
@df dat2 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
plot!(xlim=(0, 70), ylim=(-36, -31), )
modconvtplot!(sol2_sigm128, sol2_Rp128, labels=["PBM 1&8" "3PP 1&2&8"], markeralpha=0)
plot!(title="Tsh = $(Tsh1(Inf*u"s")|>u"°C")", xlabel="")
pl8 = blankplot_hrC(legend=:bottomright)
@df dat8 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
modconvtplot!(sol8_sigm128, sol8_Rp128, labels=["PBM 1&2&8" "3PP 1&2&8"], markeralpha=0)
plot!(xlim=(0, 30), ylim=(-32, -26))
plot!(title="Tsh = $(Tsh8(Inf*u"s")|>u"°C")", )
plot(pl1, pl2, pl8, layout=(3,1), size=(600, 600), legend=:outerright)
end
savefig(plotsdir("suc128_train_T.svg"))

# Validation plot

begin
pl3 = blankplot_hrC()
@df dat3 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
modconvtplot!(sol3_PBM128, sol3_Rp128, labels=["PBM 1&2&8" "3PP 1&2&8" ], markeralpha=0)
plot!(xlim=(0, 80))
plot!(title="Tsh = $(Tsh3(Inf*u"s")|>u"°C")", xlabel="")
pl5 = blankplot_hrC()
@df dat5 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
modconvtplot!(sol5_PBM128, sol5_Rp128, labels=[ "PBM 1&2&8" "3PP 1&2&8"], markeralpha=0)
plot!(xlim=(0, 55))
plot!(title="Tsh = $(Tsh5(Inf*u"s")|>u"°C")", ylabel="")
plot(pl3, pl5, layout=(2,1))
plot!(ylim=(-37, -32), legend=:outerright)
end
savefig(plotsdir("suc128train_35val_T.svg"))

# -----------------
# Inspect fit results in Rp space rather than temperature space

begin
pl1 = blankplot_lRp(ylims=(0, 12), legend=:topleft)
plot!(title="Cycle 1, Tsh = $(Tsh1(Inf*u"s")|>u"°C")")
pbmtdlas!(sol1_sigm128, dat1, Rp1_params.Rp, labels=["PBM 1&8" "3PP 1" "TDLAS"])
plot!(xlabel="")
pl2 = blankplot_lRp(ylims=(0, 12), legend=:topleft,)
plot!(title="Cycle 2, Tsh = $(Tsh2(Inf*u"s")|>u"°C")")
pbmtdlas!(sol2_sigm128, dat2, Rp2_params.Rp, labels=["PBM 1&8" "3PP 2" "TDLAS"])
pl5 = blankplot_lRp(ylims=(0, 12), legend=:topleft,)
plot!(title="Cycle 5, Tsh = $(Tsh5(Inf*u"s")|>u"°C")")
pbmtdlas!(sol5_PBM128, filter(x->x.Rp>0u"m/s", dat5), Rp5, labels=["PBM 1&8" "3PP 5" "TDLAS"])
plot!(ylabel="")
pl8 = blankplot_lRp(ylims=(0, 12), legend=:topleft,)
plot!(title="Cycle 8, Tsh = $(Tsh8(Inf*u"s")|>u"°C")")
pbmtdlas!(sol8_sigm128, dat8, Rp8_params.Rp, labels=["PBM 1&8" "3PP 8" "TDLAS"])
plot!(xlabel="", ylabel="")
plot(pl1, pl8, pl2, pl5, layout=(2,2), size=(600,500))
end
savefig(plotsdir("suc128_Rp_val.svg"))

# -----------------
# Inspect with MYEGA fit instead of tanh, since seems better-conditioned


sol2_myega128 = gen_sol_mc_arb(opt128myega.u, MC2_params18, k_f_myega);
sol3_myega128 = gen_sol_mc_arb(opt128myega.u, po_PBM3, k_f_myega);
sol5_myega128 = gen_sol_mc_arb(opt128myega.u, po_PBM5, k_f_myega);
begin
pl3 = blankplot_hrC(legend=:bottomright)
@df dat3 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
modconvtplot!(sol3_myega128, sol3_Rp128, labels=["MC" "Pikal"], markeralpha=0)
plot!(xlim=(0, 70), ylim=(-38, -31), )
plot!(title="Cycle 3, Tsh = $(Tsh3(Inf*u"s")|>u"°C")")
pl5 = blankplot_hrC(legend=:bottomright)
@df dat5 exptfplot!(:t, :T1, :T2, :T3, markershape=markers)
modconvtplot!(sol5_myega128, sol5_Rp128, labels=["MC" "Pikal"], markeralpha=0)
plot!(xlim=(0, 51), ylim=(-38, -31), )
plot!(title="Cycle 5, Tsh = $(Tsh5(Inf*u"s")|>u"°C")")
plot(pl3, pl5, layout=(2,1), size=(600, 400), legend=:outerright)
end
savefig(plotsdir("suc128_myega_val35_T.svg"))

