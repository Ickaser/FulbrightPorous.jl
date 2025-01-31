@quickactivate :FulbrightPores

const OP = OptimizationProblem
using LineSearches
alg = Optim.BFGS(linesearch=LineSearches.BackTracking())
pythonplot()

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

fitdat1 = @df dat1[12:end] PrimaryDryFit(:t, (
    filter(!isnan, :T1),
    :T2[:t .< 22u"hr"],
    :T3[:t .< 20u"hr"],
    ),
    t_end1
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

fitdat2 = @df dat2[10:end] PrimaryDryFit(:t, (
    filter(!isnan, :T1),
    filter(!isnan, :T2),
    :T3[:t .< 7.5u"hr"],
    ),
    t_end2
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

fitdat3 = @df dat3[10:end] PrimaryDryFit(:t, (
    filter(!isnan, :T1),
    :T2[:t .< 13u"hr"],
    filter(!isnan, :T3),
    ),
    t_end3
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

fitdat6 = @df dat6[20:end] PrimaryDryFit(:t, (
    filter(!isnan, :T1),
    :T2[:t .< 7u"hr"],
    filter(!isnan, :T3),
    ),
    t_end6
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

fitdat7 = @df dat7[15:end] PrimaryDryFit(:t, (
    filter(!isnan, :T1),
    filter(!isnan, :T2),
    :T3[:t .< 9u"hr"],
    ),
    t_end7
)
plot(fitdat7)
end

# -----------------
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
    (Rp_g, hf0, c_solid, ρ_solution),
    (K_shf, A_v, A_p),
    (pch1, Tsh1),
])
pomc1 = (
    (hf0, c_solid, ρ_solution),
    (K_shf, A_v, A_p),
    (pch1, Tsh1),
    (ConstPhysProp(k0), C, N, S, Rs)
)

function new_Tp(po, Tsh, pch)
    po_n = deepcopy(po)
    @reset po_n.pch = pch
    @reset po_n.Tsh = Tsh
    return po_n
end

pop2 = new_Tp(pop1, Tsh2, pch2)
pop3 = new_Tp(pop1, Tsh3, pch3)
pop6 = new_Tp(pop1, Tsh6, pch6)
pop7 = new_Tp(pop1, Tsh7, pch7)

begin
pl_all = blankplot_hrC()
plot!(fitdat1, c=1, label="1")
plot!(fitdat2, c=2, label="2")
plot!(fitdat3, c=3, label="3")
plot!(fitdat6, c=4, label="6")
plot!(fitdat7, c=5, label="7")
end

# On this basis, let's try using 1 and 3 to train and 2,6,7 to validate
# For a crystalline, using only the lowest cases to measure should not be an issue

# ----------------------
# Fit Rp
obj_Rpf = OptimizationFunction(obj_Rp, AutoForwardDiff())
@time obj_Rpf(fill(0.0, 3), (pop1, fitdat1))

obj_npf = OptimizationFunction(objn_Rp, AutoForwardDiff())
# @time ForwardDiff.gradient(x->obj_Rpf(x, (pop8, fitdat8)), [0.8, 90, 1])

cases = Table([(name="1", po=pop1, fitdat=fitdat1, dat=dat1),
         (name="2", po=pop2, fitdat=fitdat2, dat=dat2),
         (name="3", po=pop3, fitdat=fitdat3, dat=dat3),
         (name="6", po=pop6, fitdat=fitdat6, dat=dat6),
         (name="7", po=pop7, fitdat=fitdat7, dat=dat7)])

single_Rps = map(cases) do case
    (;name, po, fitdat, dat) = case
    opt = solve(OP(obj_Rpf, [0.0, 0, -5], (po, fitdat)), alg, g_tol=1e-3)
    sol = gen_sol_Rp(opt.u, po)
    @info "analyis" name opt sol
    Rp = sol.prob.p.Rp
    return (opt=opt, sol=sol, Rp=Rp)
end

cases_res = Table(cases, single_Rps)

pls = map(cases_res) do case 
    plT = blankplot_hrC()
    plot!(case.fitdat)
    modconvtplot!(case.sol)
    return plT
end
plot(pls..., layout=(5,1), size=(500, 1000))

lr = range(0u"cm", hf0, length=50)
begin
sts = Table((c=[1,2,3,4,5], marker=[:o, :square, :ltriangle, :utriangle, :diamond]))
# pl_Rp = plot(u"cm", u"cm^2*hr*Torr/g", ylabel="Rp", xlabel="l")
# for (case, st) in zip(cases_res, sts)
pls = map(cases_res, sts) do case, st
    pl_Rp = plot(u"cm", u"cm^2*hr*Torr/g", ylabel="Rp", xlabel="l")
    Tsh = case.po.Tsh(Inf*u"s") |> u"°C"
    plot!(lr, case.Rp.(lr), c=st.c, label="$(case.name) fit $(Tsh)")
    ld = FulbrightPores.calc_l(case.dat, case.po.hf0)
    plot!(ld, case.dat.Rp, c=st.c, marker=st.marker, lw=1, seriestype=:samplemarkers, step=100, label="TDLAS $(case.name)")
    plot!(ylim=(0, 20), legend=:topleft)
end
plot(pls..., layout=(3,2), size=(800,800))
end

cases_res.Rp

# -------------------------------
# Fit to MC model
# 
