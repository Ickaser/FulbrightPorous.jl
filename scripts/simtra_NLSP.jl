using DrWatson
@quickactivate :FulbrightPores
import LeastSquaresOptim

const NLSP = NonlinearLeastSquaresProblem
const OP = OptimizationProblem
LSOLM = LeastSquaresOptimJL(:lm; autodiff=:forward)
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

# tcond = (termination_condition=NonlinearSolve.AbsTerminationMode(),)
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

pb_pikal = ParamObjPikal([
    (Rp_g, hf0, c_solid, ρ_solution),
    (K_shf, A_v, A_p),
    (pch1, Tsh1),
])

pop1 = pb_pikal
pop2 = deepcopy(pop1)
@reset pop2.pch = pch2
@reset pop2.Tsh = Tsh2
pop8 = deepcopy(pop1)
@reset pop8.pch = pch8
@reset pop8.Tsh = Tsh8


function err_Rp(x, p; verbose=false)
    if any(x .< 0)
        return NaN
    end
    sol = gen_sol_Rp(x, p[1])
    return err_expT(sol, p[2], tweight=1e-4, verbose=verbose)
end
# err_Rp(x,p; verbose=false) = LyoPronto.err_expT(gen_sol_Rp(x, p[1])[1], p[2]; verbose=verbose)
@time err_Rp(fill(0.0, 3), (pop1, fitdat1))

opt1p = solve(NLSP(err_Rp, fill(0.0, 3), p=(pop1, fitdat1)), LM(); maxiters = 100, abstol=1e-2, show_trace=Val(true), tcond...)
@time opt_Rp1 = solve(OP(obj_Rpf, fill(0.0, 3), (pop1, fitdat1)),Optim.BFGS(), g_tol=1e-3)
# opt_Rp1 = solve(NLSP(obj_Rp, [0.8, 14, 1], p=(pop1, fitdat1)), LM(); reltol=1e-1, show_trace=Val(true))

sol1_Rp1 = gen_sol_Rp(opt_Rp1.u, pop1);
Rp1_params = sol1_Rp1.prob.p
plot(u"hr", u"°C", xlim=(0, 60), ylim=(-40, -30))
modconvtplot!(sol1_Rp1)
@df dat1 exptfplot!(:t, :T1, :T2, :T3)

# -------------- Provide guess values for new parameters

# typical diameter: 50 μm
d_typical = 50u"μm"
# typical total volume: 
V_typical = 4/3*π*d_typical^3
N_typical = 1/V_typical
# Good guesses for Knudsen
N = ConstPhysProp(N_typical)
S = ConstPhysProp(V_typical) # = V_typ^2*N_typ
C = 2e-1u"m^0.5"
Rs = 0.5u"cm^2*Torr*hr/g"
k0 = 1e-14u"1/s"

pb_mc = (
    (hf0, c_solid, ρ_solution),
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

# tprob = ODEProblem(pomc1)
# tsol = solve(tprob, Rodas3(autodiff=true))
# poreevolplot(tsol)
function k_f_interp_suc(logk)
    return InterpK(logk, [-45.0, -38, -32, -28]u"°C" .|> u"K", k0)
end

function err_MC_interp(x,p)
    if ~(x[5] < x[6] < x[7] < x[8])
        # count number of data points
        return fill(NaN, sum(p[2].Tf_iend) + (ismissing(p[2].Tvws) ? 0 : sum(p[2].Tvw_iend)) + !ismissing(p[2].t_end))
    end
    sol = gen_sol_mc_arb(x, p[1], k_f_interp_suc)
    return LyoPronto.err_expT(sol, p[2])
end
# opt2p = solve(NLSP(err_Rp, opt_Rp1.u, (pop2, fitdat2)), TR(); maxiters=200)
# opt_MC2 = solve(NLSP(err_MC, opt_MC1.u, (pomc2, fitdat2)), TR(); maxiters=200)

# opt8p = solve(NLSP(obj_Rp, [0.8, 14, 1], p=(pop8,fitdat8)), LM(); reltol=1e-1)
# opt_MC8 = solve(NLSP(err_MC, opt_MC2.u, p=(pomc8, fitdat8)), LM(); maxiters=300)

# ----------- Plot separately
begin
plot(u"hr", u"°C", xlim=(0, 70), ylim=(-36, -31), legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3, marker=:square)
modconvtplot!(sol1_MC1, sol1_Rp1, labels=["MC" "Pikal"],marker=:square)
end
begin
blankplot_hrC(xlim=(0, 20), ylim=(-36, -28), legend=:bottomright)
@df dat8 exptfplot!(:t, :T1, :T2, :T3, marker=:circle)
modconvtplot!(sol8_MC8, sol8_Rp8, labels=["MC" "Pikal"],marker=:circle)
end