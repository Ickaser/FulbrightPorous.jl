@quickactivate :FulbrightPores

const OP = OptimizationProblem
using LineSearches
optalg = Optim.BFGS(linesearch=LineSearches.BackTracking())

# ------------------
# Load Kv function
# Computed in scripts/kv.jl, because it is used elsewhere too
@unpack K_shf = load(datadir("exp_pro", "simtra_kv.jld2"))

begin
# ---------------
# Load cycle data
dat1raw = CSV.read(datadir("exp_raw", "Simtra_Man", "Man_Cycle1.csv"), Table)
dat1 = map(simtra_rename, dat1raw)

pch1_pir_der2 = savitzky_golay(dat1.pch_pir, 391, 3, deriv=2, rate=1/u"minute").y
t_end1 = dat1.t[argmax(pch1_pir_der2[10u"hr" .< dat1.t .< 50u"hr" ])] + 10u"hr"

pch1 = RampedVariable(100u"mTorr")
@df dat1 plot(:t, :pch_cm, ylim=(0, 200))
@df dat1 plot!(:t, :pch_pir, ylim=(0, 200))
plot!(pch1)
tendplot!(t_end1)

# There's a slight delay before shelf temperature ramps
Tsh1 = RampedVariable(uconvert.(u"K", [-15, -15, -10]u"°C"), [1u"K/minute", 0.5u"K/minute"], [11u"minute"])
@df dat1 plot(:t, :Tsh, xlim=(0, 2))
plot!(Tsh1)

@df dat1[10:end] exptfplot(:t, :T1, :T2, :T3)

fitdat1 = @df dat1[12:end] ExpFitData(:t,
    TfData(filter(!isnan, :T1)),
    TfData(:T2[:t .< 22u"hr"]),
    TfData(:T3[:t .< 20u"hr"]),
    EndTimeData(t_end1)
)
plot(fitdat1)

# ---------------
# Cycle 2
dat2raw = CSV.read(datadir("exp_raw", "Simtra_Man", "Man_Cycle2.csv"), Table)
dat2 = map(simtra_rename, dat2raw)

pch2_pir_der2 = savitzky_golay(dat2.pch_pir, 191, 3, deriv=2, rate=1/u"minute").y
t_end2 = dat2.t[argmax(pch2_pir_der2[10u"hr" .< dat2.t .< 20u"hr" ])] + 10u"hr"

pch2 = RampedVariable(160u"mTorr")
@df dat2 plot(:t, :pch_cm, ylim=(0, 200))
@df dat2 plot!(:t, :pch_pir, ylim=(0, 400))
plot!(pch2)
tendplot!(t_end2)

# There's a slight delay before shelf temperature ramps
Tsh2 = RampedVariable(uconvert.(u"K", [-15, -15, 20]u"°C"), [1u"K/minute", 0.5u"K/minute"], [10u"minute"])
@df dat2 plot(:t, :Tsh, xlim=(0, 2))
plot!(Tsh2)

@df dat2 exptfplot(:t, :T1, :T2, :T3)

fitdat2 = @df dat2[10:end] ExpFitData(:t,
    TfData(filter(!isnan, :T1)),
    TfData(filter(!isnan, :T2)),
    TfData(:T3[:t .< 7.5u"hr"]),
    EndTimeData(t_end2)
)
plot(fitdat2)

# ---------------
# Cycle 3
dat3raw = CSV.read(datadir("exp_raw", "Simtra_Man", "Man_Cycle3.csv"), Table)
dat3_rows = Table(filter(x->!ismissing(x.t), map(simtra_rename, dat3raw)))
dat3 = filter(x->x.pch_cm < 200u"mTorr", dat3_rows) # Data starts before PD
dat3.t .-= dat3.t[begin]


pch3_pir_der2 = savitzky_golay(dat3.pch_pir, 191, 3, deriv=2, rate=1/u"minute").y
t_end3 = dat3.t[argmax(pch3_pir_der2[10u"hr" .< dat3.t .< 25u"hr" ])] + 10u"hr"

pch3 = RampedVariable(180u"mTorr")
@df dat3 plot(:t, :pch_cm, ylim=(0, 200))
@df dat3 plot!(:t, :pch_pir, ylim=(0, 400))
plot!(pch3)
tendplot!(t_end3)

# There's a slight delay before shelf temperature ramps
Tsh3 = RampedVariable(uconvert.(u"K", [-15, -15, 0]u"°C"), [1u"K/minute", 0.5u"K/minute"], [15u"minute"])
@df dat3 plot(:t, :Tsh, xlim=(0, 4))
plot!(Tsh3)

@df dat3 exptfplot(:t, :T1, :T2, :T3)

fitdat3 = @df dat3[10:end] ExpFitData(:t,
    TfData(filter(!isnan, :T1)),
    TfData(:T2[:t .< 13u"hr"]),
    TfData(filter(!isnan, :T3)),
    EndTimeData(t_end3)
)
plot(fitdat3)

# ---------------
# Cycle 6
dat6raw = CSV.read(datadir("exp_raw", "Simtra_Man", "Man_Cycle6.csv"), Table)
dat6 = map(simtra_rename, dat6raw)

pch6_pir_der2 = savitzky_golay(dat6.pch_pir, 191, 3, deriv=2, rate=1/u"minute").y
t_end6 = dat6.t[argmax(pch6_pir_der2[10u"hr" .< dat6.t .< 20u"hr" ])] + 10u"hr"

pch6 = RampedVariable(100u"mTorr")
@df dat6 plot(:t, :pch_cm, ylim=(0, 200))
@df dat6 plot!(:t, :pch_pir, ylim=(0, 200))
plot!(pch6)
tendplot!(t_end6)

# There's a slight delay before shelf temperature ramps
Tsh6 = RampedVariable(uconvert.(u"K", [-15, -15, 40]u"°C"), [1u"K/minute", 0.5u"K/minute"], [4u"minute"])
@df dat6 plot(:t, :Tsh, xlim=(0, 2))
plot!(Tsh6)

@df dat6 exptfplot(:t, :T1, :T2, :T3)

fitdat6 = @df dat6[20:end] ExpFitData(:t,
    TfData(filter(!isnan, :T1)),
    TfData(:T2[:t .< 7u"hr"]),
    TfData(filter(!isnan, :T3)),
    EndTimeData(t_end6)
)
plot(fitdat6)

# ---------------
# Cycle 7
dat7raw = CSV.read(datadir("exp_raw", "Simtra_Man", "Man_Cycle7.csv"), Table)
dat7 = map(simtra_rename, dat7raw)

pch7_pir_der2 = savitzky_golay(dat7.pch_pir, 191, 3, deriv=2, rate=1/u"minute").y
t_end7 = dat7.t[argmax(pch7_pir_der2[10u"hr" .< dat7.t .< 20u"hr" ])] + 10u"hr"

pch7 = RampedVariable(110u"mTorr")
@df dat7 plot(:t, :pch_cm, ylim=(0, 200))
@df dat7 plot!(:t, :pch_pir, ylim=(0, 200))
plot!(pch7)
tendplot!(t_end7)

# There's a slight delay before shelf temperature ramps
Tsh7 = RampedVariable(uconvert.(u"K", [-15, -15, 20]u"°C"), [1u"K/minute", 0.5u"K/minute"], [4u"minute"])
@df dat7 plot(:t, :Tsh, xlim=(0, 2))
plot!(Tsh7)

@df dat7 exptfplot(:t, :T1, :T2, :T3)

fitdat7 = @df dat7[15:end] ExpFitData(:t,
    TfData(filter(!isnan, :T1)),
    TfData(filter(!isnan, :T2)),
    TfData(:T3[:t .< 9u"hr"]),
    EndTimeData(t_end7)
)
plot(fitdat7)
end

# -----------------
# Parameter setup

begin
rad_i = 1.1u"cm"
rad_o = 1.2u"cm"
A_p = π*rad_i^2
A_v = π*rad_o^2
fillvol = 5u"mL"
hf0 = fillvol/A_p

Rp_g = RpFormFit(0.8u"cm^2*Torr*hr/g", 14u"cm*Torr*hr/g", 1u"1/cm")

csolid = 0.04u"g/mL"
ρsolution = 1u"g/mL"

d_typical = 50u"μm"
V_typical = 4/3*π*d_typical^3
N_typical = 1/V_typical
# Good guesses for Knudsen
N = ConstPhysProp(N_typical)
S = ConstPhysProp(V_typical) # = V_typ^2*N_typ
C = 2e-1u"m^0.5"
Rs = 0.5u"cm^2*Torr*hr/g"
k0 = 1e-14u"1/s"

pop1 = ParamObjPikal([
    (Rp_g, hf0, csolid, ρsolution),
    (K_shf, A_v, A_p),
    (pch1, Tsh1),
])
pomc1 = ParamObjMC([
    (hf0, csolid, ρsolution),
    (K_shf, A_v, A_p),
    (pch1, Tsh1),
    (ConstPhysProp(k0), C, N, S, Rs)
])
end

function new_Tp(po, Tsh, pch)
    po_n = copy(po)
    @reset po_n.pch = pch
    @reset po_n.Tsh = Tsh
    return po_n
end

pop2 = new_Tp(pop1, Tsh2, pch2)
pop3 = new_Tp(pop1, Tsh3, pch3)
pop6 = new_Tp(pop1, Tsh6, pch6)
pop7 = new_Tp(pop1, Tsh7, pch7)
pomc2 = new_Tp(pomc1, Tsh2, pch2)
pomc3 = new_Tp(pomc1, Tsh3, pch3)
pomc6 = new_Tp(pomc1, Tsh6, pch6)
pomc7 = new_Tp(pomc1, Tsh7, pch7)

begin
pl_all = blankplot_hrC()
plot!(fitdat1, c=1, labsuffix=" 1")
plot!(fitdat2, c=2, labsuffix=" 2")
plot!(fitdat3, c=3, labsuffix=" 3")
plot!(fitdat6, c=4, labsuffix=" 6")
plot!(fitdat7, c=5, labsuffix=" 7")
end

# On this basis, let's try using 1 and 3 to train and 2,6,7 to validate
# For a crystalline, using only the lowest cases to measure should not be an issue

# ----------------------
# Fit Rp
obj_Rpf = OptimizationFunction(obj_pd, AutoForwardDiff())
@time obj_Rpf(fill(0.0, 3), (pop1, fitdat1))

obj_npf = OptimizationFunction(objn_pd, AutoForwardDiff())

cases = Table([ (name="1", pop=pop1, pomc=pomc1, fitdat=fitdat1, dat=dat1),
                (name="2", pop=pop2, pomc=pomc2, fitdat=fitdat2, dat=dat2),
                (name="3", pop=pop3, pomc=pomc3, fitdat=fitdat3, dat=dat3),
                (name="6", pop=pop6, pomc=pomc6, fitdat=fitdat6, dat=dat6),
                (name="7", pop=pop7, pomc=pomc7, fitdat=fitdat7, dat=dat7)])

single_Rps = map(cases) do case
    (;name, pop, fitdat, dat) = case
    opt = solve(OP(obj_Rpf, [0.0, 0, -5], (pop, fitdat)), optalg, g_tol=1e-3)
    sol = gen_sol_Rp(opt.u, pop)
    @info "analyis" name opt sol
    Rp = sol.prob.p.Rp
    return (opt=opt, sol_p=sol, Rp=Rp)
end

cases_Rps = Table(cases, single_Rps)

pls = map(cases_Rps) do case 
    plT = blankplot_hrC()
    plot!(case.fitdat)
    modconvtplot!(case.sol_p)
    return plT
end
plot(pls..., layout=(5,1), size=(500, 1000))

lr = range(0u"cm", hf0, length=50)
begin
sts = Table((c=2:6, marker=[:o, :square, :ltriangle, :utriangle, :diamond]))
pl_Rp = plot(u"cm", u"cm^2*hr*Torr/g", ylabel="Rp", xlabel="l", color_palette = palette(:magma, 7))
for (case, st) in zip(sort(cases_Rps, by=x->maximum(x.fitdat.Tfs[2])), sts)
# pls = map(cases_Rps, sts) do case, st
    # pl_Rp = plot(u"cm", u"cm^2*hr*Torr/g", ylabel="Rp", xlabel="l")
    Tsh = case.pop.Tsh(Inf*u"s") |> u"°C"
    plot!(lr, case.Rp.(lr), c=st.c, label="$(case.name) fit; Tsh=$(Tsh)")
    ld = FulbrightPores.calc_l(case.dat, case.pop.hf0)
    plot!(ld, case.dat.Rp, c=st.c, marker=st.marker, markersize=7, la=0.5, lw=1, seriestype=:samplemarkers, step=100, label="$(case.name) TDLAS")
    plot!(ylim=(0, 15), legend=:topleft)
end
# plot(pls..., layout=(3,2), size=(800,800))
pl_Rp
end
savefig(plotsdir("simtra_mannitol_plainRp.svg"))
savefig(plotsdir("simtra_mannitol_plainRp.pdf"))

cases_Rps.Rp

# -------------------------------
# Fit to MC model
# Fit to 1 and 3, validate on 2, 6, 7

function make_mcpass(casenums)
    subdat = filter(x->parse(Int, x.name)∈casenums, cases)
    (subdat.pomc, subdat.fitdat)
end

pass13 = make_mcpass([1, 2])
pass123 = make_mcpass([1, 2, 3])
pass16 = make_mcpass([1, 6])
pass136 = make_mcpass([1, 3, 6])
pass17 = make_mcpass([1, 7])

myega_scale = (1e10u"s^-1", 1.6e3u"K", 1e2u"K")
k_f_myega(vec) = MYEGA_k((exp.(vec).*myega_scale))
objn_myega(fitprm, pf) = objn_mc(fitprm, (pf..., k_f_myega), tweight=1e-1)
objnf_myega = OptimizationFunction(objn_myega, AutoForwardDiff())
lb_myega = vcat(fill(0.0, 4), fill(-100, 3))
ub_myega = vcat(fill(Inf, 4), fill( 100, 3))

sigm_scale = (1e-15u"s^-1", 1e-9u"s^-1", 0.8u"K", 239u"K")
function k_f_sigm(vec)
    p123 = exp.(vec[1:3] .-1).*sigm_scale[1:3]
    p4 = vec[4]*u"K" + sigm_scale[4]
    Sigmoid_k(p123..., p4)
end
objn_sigm(fitprm, pf) = objn_mc(fitprm, (pf..., k_f_sigm), tweight=1e-1)
objnf_sigm = OptimizationFunction(objn_sigm, AutoForwardDiff())
lb_sigm = vcat(fill(0.0, 4), fill(-30, 4))
# lb_sigm[8] = -30
ub_sigm = vcat(fill(Inf, 4), fill( 30, 4))
# ub_sigm[8] = 30

# opt13sigm = solve(OP(objnf_sigm, fill(1.0, 8), pass13), optalg, maxiters=50, time_limit=15)
opt13sigm = solve(OP(objnf_sigm, fill(1.0, 8), pass13, lb=lb_sigm, ub=ub_sigm), optalg, maxiters=50, time_limit=15)
opt16sigm = solve(OP(objnf_sigm, fill(1.0, 8), pass16, lb=lb_sigm, ub=ub_sigm), optalg, maxiters=50, time_limit=15)
opt123sigm = solve(OP(objnf_sigm, fill(0.1, 8), pass123, lb=lb_sigm, ub=ub_sigm), optalg, maxiters=50, time_limit=15)
opt136sigm = solve(OP(objnf_sigm, fill(1.0, 8), pass136, lb=lb_sigm, ub=ub_sigm), optalg, maxiters=50, time_limit=15)
# @profview solve(OP(objnf_sigm, fill(1.0, 8), pass136, lb=lb_sigm, ub=ub_sigm), optalg, maxiters=50, time_limit=15)
opt123myega = solve(OP(objnf_myega, fill(0.1, 7), pass123, lb=lb_myega, ub=ub_myega), optalg, maxiters=50, time_limit=15)
opt136myega = solve(OP(objnf_myega, fill(1.0, 7), pass136, lb=lb_myega, ub=ub_myega), optalg, maxiters=50, time_limit=15)

fits = Table((opt=[opt13sigm, opt16sigm, opt123sigm, opt136sigm, opt123myega, opt136myega], 
              k_type=[k_f_sigm, k_f_sigm, k_f_sigm, k_f_sigm, k_f_myega, k_f_myega],
              cases=[(1,3),  (1,6),  (1,2,3),  (1,3,6),  (1,2,3),  (1,3,6)]))
fits = Table(fits, (kf= map(x->x.k_type(x.opt[5:end]), fits),))

# cases
function gen_pbmsols(opt, k_f_gen)
    map(cases) do case
        sol = gen_sol_mc_arb(opt.u, case.pomc, k_f_gen)
        return (sol_mc=sol,k_f=sol.prob.p.k_f)
    end
end

cases_both = Table(cases, single_Rps, gen_pbmsols(opt123myega, k_f_myega))
# cases_both = Table(cases, single_Rps, gen_pbmsols(opt123sigm, k_f_sigm))

begin
pythonplot()
pls_T = map(cases_both) do case
    pl = blankplot_hrC(legend=:outerright, title=case.name)
    plot!(case.fitdat)
    modconvtplot!(case.sol_mc, case.sol_p, linealpha=0.7, markeralpha=0, labels=["PBM: M 1&2&3" "3PP $(case.name)"])
    return pl
end
# plot!(pls_T[5], legend=:bottomright)
plot(pls_T..., layout=(3,2), size=(800, 600))
end
savefig(plotsdir("simtra_mannitol_T.svg"))
savefig(plotsdir("simtra_mannitol_T.pdf"))

begin
gr()
plk = plot(u"°C", u"s^-1", ylabel="k", xlabel="T")
for fit in fits
    symb = fit.k_type == k_f_sigm ? :circle : :square
    c = findfirst([fit.cases] .== [(1,2,3), (1,3,6), (1,3), (1,6)]) 
    kineticplot!(fit.kf, T=range(230u"K", 265u"K", length=101), label="$(fit.cases), $(fit.k_type)", c=c, marker=symb)
end
plot!(legend=:outerright, ylim=(1e-25, 1e-8))
end
savefig(plotsdir("man_kinetic.svg"))

# begin
# pl = blankplot_hrC()
# modconvtplot!(sol1_mc13, sol3_mc13)
# plot!(fitdat1, c=1, labsuffix=" 1")
# plot!(fitdat3, c=3, labsuffix=" 3")
# end

begin
pls = map(cases_both) do case
    pl = blankplot_lRp(ylim=(0, 15), legend=:bottomright, title=case.name)
    pbmtdlas!(case.sol_mc, case.dat, case.Rp, labels=["PBM 1&2&3" "3PP $(case.name)" "TDLAS $(case.name)"])
end
plot(pls..., layout=(3,2), size=(800, 600))
end
savefig(plotsdir("simtra_mannitol_Rp.svg"))
savefig(plotsdir("simtra_mannitol_Rp.pdf"))

kineticplot(sols_sigm136[1].sol_mc, T=range(230u"K", 250u"K", length=101))
# for case in cases_both

begin

poreevolplot( cases_both[1].sol_mc, label="1", c=2, palette=palette(:magma, 7))
poreevolplot!(cases_both[2].sol_mc, label="2", c=5)
poreevolplot!(cases_both[3].sol_mc, label="3", c=3)
poreevolplot!(cases_both[4].sol_mc, label="6", c=6)
poreevolplot!(cases_both[5].sol_mc, label="7", c=4)
end
savefig(plotsdir("man_pores.svg"))
savefig(plotsdir("man_pores.pdf"))
