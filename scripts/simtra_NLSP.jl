using DrWatson
@quickactivate :FulbrightPorous
import LeastSquaresOptim

const NLSP = NonlinearLeastSquaresProblem
const OP = OptimizationProblem
LSOLM = LeastSquaresOptimJL(:lm; autodiff=:forward)
# pythonplot()


begin
# ------------------
# Load Kv function
# Computed in scripts/kv.jl, because it is used elsewhere too
Kshf = load(datadir("exp_pro", "simtra_kv.jld2"))["K_shf"]

# --------------------------
# Cycle 1 with original filenames, but labeled as case 4 in Table 8: use that numbering
dat4raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle4.csv"), Table)
dat4 = map(simtra_rename, dat4raw)

t_end4 = identify_pd_end(dat4.t, dat4.pch_pir, Val(:onoff), window_width=391,  tmin=10u"hr",  tmax=85u"hr")

pch4 = RampedVariable(100u"mTorr")
# @df dat4 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat4 plot(:t, :pch_pir, ylim=(0, 200))
# plot!(pch4)
# tendplot!(t_end4)

# There's a slight delay before shelf temperature ramps
Tsh4 = RampedVariable(uconvert.(u"K", [-40, -40, -25]u"°C"), [1u"K/minute", 0.5u"K/minute"], [23u"minute"])
# @df dat4 plot(:t, :Tsh, xlim=(0, 2))
# plot!(Tsh4)

fitdat4 = ExpFitData(dat4.t, 
    TfData(filter(x->(!isnan(x) && x<-32.5u"°C"), dat4.T1)),
    TfData(filter(!isnan, dat4.T2[dat4.t .< 60u"hr"])),
    # TfData(collect(filter(!isnan, dat4.T3))),
    EndTimeData(t_end4)
)
# plot(fitdat4)

# ------------------------------
# Cycle 6 according to publication naming

dat6raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle6.csv"), Table)
dat6 = map(simtra_rename, dat6raw)
filter!(x->x.t .> 0u"hr", dat6)

# @df dat6 exptfplot!(:t, :T1, :T2, :T3, c=:red)
# plot!(ylim=(-40, -30))

t_end6 = identify_pd_end(dat6.t, dat6.pch_pir, Val(:onoff), window_width=391,  tmin=10u"hr",  tmax=120u"hr")
pch6 = RampedVariable(100u"mTorr")
# @df dat6 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat6 plot(:t, :pch_pir, ylim=(0, 200))
# plot!(pch6)
# tendplot!(t_end6)

fitdat6 = @df dat6 ExpFitData(:t,
    # TfData(collect(filter(!isnan, :T1))),
    TfData(filter(!isnan, :T2[:t .< 85u"hr"])),
    TfData(filter(!isnan, :T3)),
    EndTimeData(t_end6))
# plot(fitdat6)

Tsh6 = RampedVariable(uconvert.(u"K", [-45, -45, -30]u"°C"), [1u"K/minute", 0.5u"K/minute"], [23u"minute"])
# @df dat6 plot(:t, :Tsh, xlim=(0, 2))
# plot!(Tsh6)

# -------------------------------  
# Cycle 5
dat5raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle5.csv"), Table)
dat5 = map(simtra_rename, dat5raw)

pch5 = RampedVariable(60u"mTorr")
t_end5 = identify_pd_end(dat5.t, dat5.pch_pir, Val(:onoff), window_width=391,  tmin=10u"hr",  tmax=68u"hr")
# @df dat5 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat5 plot!(:t, :pch_pir, ylim=(0, 200))
# plot!(pch5)
# tendplot!(t_end5)

Tsh5 = RampedVariable(uconvert.(u"K", [-40, -25]u"°C"), 0.5u"K/minute")
# @df dat5 plot(:t, :Tsh, xlim=(0, 2))
# plot!(Tsh5)
fitdat5 = @df dat5 ExpFitData(:t,
    TfData(filter(!isnan, :T1)),
    TfData(filter(!isnan, :T2)),
    TfData(:T3[:t.<45u"hr"]),
    EndTimeData(t_end5))
# plot(fitdat5)

# -----------------
# Cycle 5
dat2raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle2.csv"), Table)
dat2 = map(simtra_rename, dat2raw)

pch2 = RampedVariable(50u"mTorr")
t_end2 = identify_pd_end(dat2.t, dat2.pch_pir, Val(:onoff), window_width=391,  tmin=10u"hr",  tmax=50u"hr")
# @df dat2 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat2 plot!(:t, :pch_pir, ylim=(0, 200))
# plot!(pch2)
# tendplot!(t_end2)

Tsh2 = RampedVariable(uconvert.(u"K", [-40, -40, -20]u"°C"), [1u"K/minute", 0.5u"K/minute"], [12u"minute"])
# @df dat2 plot(:t, :Tsh_i, xlim=(0, 2))
# plot!(Tsh2)

fitdat2 = @df dat2 ExpFitData(:t, 
    TfData(filter(!isnan, :T1)),
    TfData(filter(!isnan, :T3)),
    TfData(:T2[:t.<40u"hr"]),
    EndTimeData(t_end2))
# plot(fitdat2)

# -----------------
# Cycle 1, high temperature
dat1raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle1.csv"), Table)
dat1 = map(simtra_rename, dat1raw)
# @df dat1 exptfplot(:t, :T1, :T2, :T3)
# @df dat1 plot!(:t, :Tsh_i, xlim=(0, 20), ylim=(-40, -20))

pch1 = RampedVariable(100u"mTorr")
t_end1 = identify_pd_end(dat1.t, dat1.pch_pir, Val(:onoff), window_width=391,  tmin=10u"hr",  tmax=41u"hr")
# @df dat1 plot(:t, :pch_cm, ylim=(0, 200))
# @df dat1 plot(:t, :pch_pir, ylim=(0, 200))
# plot!(pch1)
# tendplot!(t_end1)

Tsh1 = RampedVariable(uconvert.(u"K", [-40, -40, -10]u"°C"), [1u"K/minute", 0.5u"K/minute"], [4u"minute"])
# @df dat1 plot(:t, :Tsh_i, xlim=(0, 2))
# plot!(Tsh1)

fitdat1 = @df dat1 ExpFitData(:t,
    TfData(filter(!isnan, :T1)),
    TfData(filter(!isnan, :T2[:t .< 17u"hr"])),
    TfData(filter(!isnan, :T3[:t .< 17u"hr"])),
    EndTimeData(t_end1))
# plot(fitdat1)
end

# tcond = (termination_condition=NonlinearSolve.AbsTerminationMode(),)
# -------------
# Parameter setup

begin
rad_i = 1.1u"cm"
rad_o = 1.2u"cm"
Ap = π*rad_i^2
Av = π*rad_o^2
fillvol = 5u"mL"
hf0 = fillvol/Ap

Rp_g = RpFormFit(0.8u"cm^2*Torr*hr/g", 14u"cm*Torr*hr/g", 1u"1/cm")

csolid = 0.04u"g/mL"
ρsolution = 1u"g/mL"

# pb_pikal = ParamObjPikal([
#     (Rp_g, hf0, csolid, ρsolution),
#     (K_shf, A_v, A_p),
#     (pch1, Tsh1),
# ])

# pop1 = pb_pikal
# pop2 = copy(pop1)
# @reset pop2.pch = pch2
# @reset pop2.Tsh = Tsh2
# pop8 = copy(pop1)
# @reset pop8.pch = pch8
# @reset pop8.Tsh = Tsh8


# @time opt_Rp1 = solve(OP(obj_Rpf, fill(0.0, 3), (pop1, fitdat1)),Optim.BFGS(), g_tol=1e-3)
# # opt_Rp1 = solve(NLSP(obj_Rp, [0.8, 14, 1], p=(pop1, fitdat1)), LM(); reltol=1e-1, show_trace=Val(true))

# sol1_Rp1 = gen_sol_Rp(opt_Rp1.u, pop1);
# Rp1_params = sol1_Rp1.prob.p
# plot(u"hr", u"°C", xlim=(0, 60), ylim=(-40, -30))
# modconvtplot!(sol1_Rp1)
# @df dat1 exptfplot!(:t, :T1, :T2, :T3)

# # -------------- Provide guess values for new parameters

# # typical diameter: 50 μm
d_typical = 50u"μm"
V_typical = 4/3*π*d_typical^3
N_typical = 1/V_typical
# Dummy values
N = ConstPhysProp(N_typical)
S = ConstPhysProp(V_typical) # = V_typ^2*N_typ
C = 2e-1u"m^0.5"
Rs = 0.5u"cm^2*Torr*hr/g"
k0 = 1e-14u"1/s"

pop1 = ParamObjPikal([
    (Rp_g, hf0, csolid, ρsolution),
    (Kshf, Av, Ap),
    (pch1, Tsh1),
])
pomc1 = ParamObjMC([
    (hf0, csolid, ρsolution),
    (Kshf, Av, Ap),
    (pch1, Tsh1),
    (ConstPhysProp(k0), C, N, S, Rs)
])

function new_Tp(po, Tsh, pch)
    po_n = deepcopy(po)
    @reset po_n.pch = pch
    @reset po_n.Tsh = Tsh
    return po_n
end
pop2 = new_Tp(pop1, Tsh2, pch2)
pop3 = new_Tp(pop1, Tsh3, pch3)
pop5 = new_Tp(pop1, Tsh5, pch5)
pop8 = new_Tp(pop1, Tsh8, pch8)
pomc2 = new_Tp(pomc1, Tsh2, pch2)
pomc3 = new_Tp(pomc1, Tsh3, pch3)
pomc5 = new_Tp(pomc1, Tsh5, pch5)
pomc8 = new_Tp(pomc1, Tsh8, pch8)
dat1 = Table(dat1, (ld=FulbrightPorous.calc_l(dat1, pop1.hf0)))
dat2 = Table(dat2, (ld=FulbrightPorous.calc_l(dat2, pop2.hf0)))
dat3 = Table(dat3, (ld=FulbrightPorous.calc_l(dat3, pop3.hf0)))
dat5 = Table(dat5, (ld=FulbrightPorous.calc_l(dat5, pop5.hf0)))
dat8 = Table(dat8, (ld=FulbrightPorous.calc_l(dat8, pop8.hf0)))
end

cases = Table([ (name="1", pop=pop1, pomc=pomc1, fitdat=fitdat1, dat=dat1),
                (name="2", pop=pop2, pomc=pomc2, fitdat=fitdat2, dat=dat2),
                (name="3", pop=pop3, pomc=pomc3, fitdat=fitdat3, dat=dat3),
                (name="5", pop=pop5, pomc=pomc5, fitdat=fitdat5, dat=dat5),
                (name="8", pop=pop8, pomc=pomc8, fitdat=fitdat8, dat=dat8),
                ])

trans_Rp = Rp_transform_basic(0.8u"cm^2*Torr*hr/g", 14.0u"cm*Torr*hr/g", 1e0u"1/cm")
obj_Rpf = OptimizationFunction((x,y)->obj_pd(x,y,tweight=1e-2), AutoForwardDiff(chunksize=3))
sqe_Rpf = NonlinearFunction{false}((x,y)->nls_pd(x,y,tweight=1e-2))
# obj_Rpf = OptimizationFunction(obj_pd, AutoForwardDiff(chunksize=3))
@time obj_Rpf([0.0, 1.0, 0.0], (trans_Rp, pop1, fitdat1))
@time sqe_Rpf([0.0, 1.0, 0.0], (trans_Rp, pop1, fitdat1))

@time opt1p = solve(NLSP(sqe_Rpf, fill(0.0, 3), p=(trans_Rp, pop1, fitdat1)), LevenbergMarquardt(); maxiters = 100, abstol=1e-2, show_trace=Val(true))
obj_npf = OptimizationFunction((x,y)->objn_pd(x,y,tweight=1e-2), AutoForwardDiff(chunksize=3))

single_Rps = map(cases) do case
    (;name, pop, fitdat, dat) = case
    u0 = [0., 0., 0.]
    @time opt = solve(OP(obj_Rpf, u0, (trans_Rp, pop, fitdat)), optalg)#, g_tol=1e-6)
    Rp = RpFormFit(transform(trans_Rp, opt.u).Rp...)
    sol = gen_sol_pd(opt.u, trans_Rp, pop, fitdat)
    @info "analyis" name opt Rp
    hRps = map(x->calc_hRp_T(case.pop, case.fitdat, i=x), 1:length(case.fitdat.Tfs))
    return (opt=opt, sol_p=sol, Rp, hRps)
end;

cases_Rps = Table(cases, single_Rps)

pls = map(cases_Rps) do case 
    plT = blankplot_hrC()
    plot!(case.fitdat, title=case.name)
    modconvtplot!(case.sol_p)
    plot!(legend = :bottomright)
    # plot!(case.pop.Tsh, tmax=case.fitdat.t_end, c=:black, label=L"T_\mathrm{sh}")
    return plT
end
plot(pls..., layout=(3,2), size=(800, 1000))

l = range(0u"cm", hf0, length=100)
begin
pls = map(cases_Rps) do case
    plr = blankplot_lRp()
    plot!(case.dat.ld, case.dat.Rp, label="TDLAS")
    for i in 1:length(case.fitdat.Tfs)
        hdT, RpT = case.hRps[i]
        plot!(hdT, RpT, c=i+1, label="Rp from TP$i")
    end
    plot!(l, case.Rp.(l), c=5, label="3PP fit directly to T")
    plot!(ylim=(0, 10), legend = :bottomright)
    plot!(title="\$p_{ch}=\\mathrm{$(case.pop.pch(Inf*u"s")|>u"mTorr")}\$,\n \$T_{sh}=$(float(case.pop.Tsh(Inf*u"s")|>u"°C"))\$  ")
    plot!(top_margin=20Plots.px)
    plot!(bottom_margin=15Plots.px)
    case.name ∈ ["1", "5"] && plot!(left_margin=20Plots.px)
    plot!(legend=:none)
    return plr
end
explot = plot(framestyle=:none, legendfontsize=12, legend_position=:inside)
plot!([Inf], [Inf], label="TDLAS", c=1)
plot!([Inf], [Inf], label="Rp from TP1", c=2)
plot!([Inf], [Inf], label="Rp from TP2", c=3)
plot!([Inf], [Inf], label="Rp from TP3", c=4)
plot!([Inf], [Inf], label="3PP fit directly to T", c=5)
plot(pls..., explot, layout=(2,3), size=(900, 500))
end
savefig(plotsdir("simtra_suc_Rp.png"))
savefig(plotsdir("simtra_suc_Rp.svg"))
savefig(plotsdir("simtra_suc_Rp.pdf"))

function make_Rp_pass(casenums)
    subdat = filter(x->parse(Int, x.name)∈casenums, cases_Rps)
    (subdat.pop, subdat.fitdat)
end
function make_mc_pass(casenums)
    subdat = filter(x->parse(Int, x.name)∈casenums, cases)
    (subdat.pomc, subdat.fitdat)
end

rpfit_cases = Table(cases = [(1,), (2,), (3,), (5,), (8,), (1,8), (1,2,8)],
        p0 = fill([0.5, -0.3, 0], 7),
        trans = fill(trans_Rp, 7))
rpfits = map(rpfit_cases) do case
    @time solve(OP(obj_npf, case.p0, (case.trans, make_Rp_pass(case.cases)...)), optalg;)
end

# -----------------------------
# Fit k by interpolating at four values
# interp_nt = (kl1 = TV.Identity(),
#             kl2 = TV.Identity(),
#             kl3 = TV.Identity(),
#             kl4 = TV.Identity())
const TV = TransformVariables
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