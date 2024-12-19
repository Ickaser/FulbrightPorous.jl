using DrWatson
@quickactivate :FulbrightPores
using DataInterpolations
using UnitfulLatexify, Latexify
using PrettyTables

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
optKv = optimize(x->fitKv(x)[1], [10.0, 1, 1])
K_shf = fitKv(optKv.minimizer)[2]

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

fitdat2 = PrimaryDryFit(collect(dat2.t), (
    # collect(filter(!isnan, dat2.T1)),
    collect(filter(!isnan, dat2.T2[dat2.t .< 85u"hr"])),
    collect(filter(!isnan, dat2.T3)),))
plot(fitdat2)

pch2 = RampedVariable(100u"mTorr")
# @df dat2 plot(:t, :pch_cm, ylim=(0, 200))
# plot!(pch2)
Tsh2 = Tsh1 - 5u"K"
# @df dat2 plot(:t, :Tsh, xlim=(0, 2))
# plot!(Tsh2)

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
# Cycle 8
dat8raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle8.csv"), Table)
dat8 = map(simtra_rename, dat8raw)
@df dat8 exptfplot(:t, :T1, :T2, :T3)

pch8 = RampedVariable(100u"mTorr")
@df dat8 plot(:t, :pch_cm, ylim=(0, 200))
plot!(pch8)
Tsh8 = RampedVariable(uconvert.(u"K", [-40, -40, -10]u"°C"), [1u"K/minute", 0.5u"K/minute"], [4u"minute"])
@df dat8 plot(:t, :Tsh_i, xlim=(0, 2))
plot!(Tsh8)

# -------------
# Model setup, Pikal

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

# model setup, additional parameters for PBM

# typical diameter: 50 μm
d_typical = 50u"μm"
# typical total volume: 
V_typical = 4/3*π*d_typical^3
N_typical = fillvol/V_typical

N = RampedVariable(N_typical / hf0)
S = RampedVariable(V_typical^2*N_typical/hf0)
C = 2e-1u"m^0.5"
R0 = 0.01u"cm^2*Torr*hr/g"
function k_act(T, p)
    k0, K, C = p
    x = (T-K)/C
    return k0*(exp(x) + exp(-x))/exp(-x) 
end
k0 = 1e-16u"1/s"
Tact = -37u"°C"
θ = 0.4u"K"
kf = T->k_act(T, [k0, Tact, θ])
pb_PBM = (
    (Rp_g, hf0, c_solid, ρ_solution),
    (K_shf, A_v, A_p),
    (pch1, Tsh1),
    (kf, C, N, S, R0)
)
poPBM1 = ParamObjMC(pb_PBM)
poPBM2 = deepcopy(poPBM1)
@reset poPBM2.pch = pch2
@reset poPBM2.Tsh = Tsh2

# ---------------- Functions needed for fitting

# One at a time, 3PP
function gen_sol_Rp(a123, po)
    Rp_un = [u"cm^2*Torr*hr/g", u"cm*Torr*hr/g", u"1/cm"]
    new_params = @set po.Rp = RpFormFit((a123.*Rp_un)...)
    prob = ODEProblem(new_params)
    sol = solve(prob, Rodas3())
    return sol, new_params
end
# Two at at a time, 3PP
function gen2sol_Rp(fitprm, po1, po2)
    sol1, prm1 = gen_sol_Rp(fitprm, po1)
    sol2, prm2 = gen_sol_Rp(fitprm, po2)
    return sol1, sol2, prm1, prm2
end
# One at a time, PBM
fit_un_PBM = [u"cm^-1", u"cm^5", u"s^-1", u"K", u"K", u"m^0.5", u"m/s"]
function gen_sol_PBM(fitprmlog, po; verbose=false)
    new_params = deepcopy(po)
    fitprm = exp.(fitprmlog) # Constrain all to be positive by taking exp
    verbose && @info "fitprm" fitprm
    @reset new_params.N = RampedVariable(fitprm[1]*fit_un_PBM[1])
    @reset new_params.S = RampedVariable(fitprm[2]*fit_un_PBM[2])
    @reset new_params.k_f = T->k_act(T, fitprm[3:5].*fit_un_PBM[3:5])
    @reset new_params.C  = fitprm[6]*fit_un_PBM[6]
    @reset new_params.R0 = fitprm[7]*fit_un_PBM[7]
    prob = ODEProblem(new_params)
    sol = solve(prob, Rodas3())
    return sol, new_params
end
# Two at a time, PBM
function gen2sol_PBM(fitprmlog, po1, po2)
    sol1, prm1 = gen_sol_PBM(fitprmlog, po1)
    sol2, prm2 = gen_sol_PBM(fitprmlog, po2)
    return sol1, sol2, prm1, prm2
end

# Fit Pikal model separately, then together
obj_Rp1 = x->obj_expT(gen_sol_Rp(x, pop1)[1], fitdat1)
obj_Rp2 = x->obj_expT(gen_sol_Rp(x, pop2)[1], fitdat2)
function obj_Rp12(fitrpm)
    sol1, sol2, prm1, prm2 = gen2sol_Rp(fitrpm, pop1, pop2)
    obj1 = obj_expT(sol1, fitdat1)
    obj2 = obj_expT(sol2, fitdat2)
    return obj1 + obj2
end

p0 = [0.8, 14, 1]
opt_Rp1 = optimize(obj_Rp1, p0, LevenbergMarquardt())
opt_Rp1 = optimize(obj_Rp1, p0, Optim.Options(g_tol=1e-3))
sol1_Rp1, Rp1_params = gen_sol_Rp(Optim.minimizer(opt_Rp1), pop1);
opt_Rp2 = optimize(obj_Rp2, p0, Optim.Options(g_tol=1e-3))
sol2_Rp2, Rp2_params = gen_sol_Rp(Optim.minimizer(opt_Rp2), pop2);
opt_Rp12 = optimize(obj_Rp12, p0, Optim.Options(g_tol=1e-3))
sol1_Rp12, sol2_Rp12, Rp12_params1, Rp12_params2 = gen2sol_Rp(Optim.minimizer(opt_Rp12), pop1, pop2);

# Plot results
begin
pl1 = blankplot_hrC()
plot!(fitdat1)
modconvtplot!(sol1_Rp1, sol1_Rp12, labels=["3PP 1, 1" "3PP 1&2, 1"], linealpha=0.8)
pl2 = blankplot_hrC()
plot!(fitdat2)
modconvtplot!(sol2_Rp2, sol2_Rp12, labels=["3PP 2, 2" "3PP 1&2, 2"], linealpha=0.8)
plot!(pl2, ylabel="")
plot(pl1, pl2, legend=:bottomright, ylim=(-38, -32), layout=(1,2), link=:all)
end

# Fit PBM model to both cases at once 
function obj12PBM(fitprmlog; verbose=true)
    verbose && @info "fitprm" exp.(fitprmlog)
    if exp(fitprmlog[4]) < 220 || exp(fitprmlog[4]) > 270
        verbose && @info "temperature unreasonable" exp(fitprmlog[4])
        return NaN 
    end
    sol1, sol2, prm1, prm2 = gen2sol_PBM(fitprmlog, poPBM1, poPBM2)
    # if sol1.retcode != ReturnCode.Terminated || sol2.retcode != ReturnCode.Terminated
    #     verbose && @info "unfinished drying" sol1.retcode sol2.retcode
    #     return NaN
    # end
    obj1 = obj_expT(sol1, fitdat1)
    obj2 = obj_expT(sol2, fitdat2)
    verbose && @info "obj" obj1 obj2
    return obj1 + obj2
end

g12_PBM = log.(ustrip.(fit_un_PBM, [N(0), S(0), k0, Tact, θ, C, R0]))
obj12PBM(g12_PBM)
# opt12_PBM = optimize(obj12PBM, g12_PBM, NelderMead(), Optim.Options(g_tol = 1e-3, time_limit=30))
opt12_PBM = optimize(obj12PBM, g12_PBM, LBFGS(), Optim.Options(g_tol = 1e-3, time_limit=30))

sol1_PBM12, sol2_PBM12, PBM1_params12, PBM2_params12 = gen2sol_PBM(Optim.minimizer(opt12_PBM), poPBM1, poPBM2);

begin
pl1 = blankplot_hrC(legend=:bottomright,)
@df dat1 exptfplot!(:t, :T1, :T2, :T3, markershape=:circle)
modconvtplot!(sol1_PBM12, sol1_Rp12, labels=["PBM 1&2" "3PP 1&2"], markershape=:circle)
pl2 = blankplot_hrC(legend=:topleft)
@df dat2 exptfplot!(:t, :T1, :T2, :T3, markershape=:square)
modconvtplot!(sol2_PBM12, sol2_Rp12, labels=["PBM 1&2" "3PP 1&2"], markershape=:square)
pl12 = plot(pl1, pl2,  ylim=(-36, -31), layout=(1,2), link=:all)
end
savefig(plotsdir("suc12_T.svg"))

l = range(0u"cm", hf0, length=50)
Rp1 = Rp1_params.Rp.(l)
Rp2 = Rp2_params.Rp.(l)
Rp12 = Rp12_params1.Rp.(l)

default(:markersize, 7)
begin
colors = [:red :blue :purple :brown]
rpattr = (marker=:utriangle, seriestype=:samplemarkers, linestyle=:dash)
poreevolplot(sol1_PBM12, label="PBM 1&2, Cycle 1", c=colors[3], marker=:circle)
poreevolplot!(sol2_PBM12, label="PBM 1&2, Cycle 2", c=colors[3], marker=:square)
plot!(l, Rp1, subplot=4, label="3PP 1"; c=colors[1], rpattr...)
plot!(l, Rp2, subplot=4, label="3PP 2"; c=colors[2], rpattr...)
plot!(l, Rp12, subplot=4, label="3PP 1&2"; c=colors[3], rpattr...)
plot!(size=(800, 400))
end
savefig(plotsdir("suc12_pores.svg"))

# -----------------
function namefitparam_PBM(fitlog, label)
    fit = exp.(fitlog)
    nt = (label = label,
        N =  fit[1]*fit_un_PBM[1],
        S =  fit[2]*fit_un_PBM[2],
        C =  fit[6]*fit_un_PBM[6],
        R0 = fit[7]*fit_un_PBM[7],
        k0 = fit[3]*fit_un_PBM[3],
        Tact = fit[4]*fit_un_PBM[4],
        θ  = fit[5]*fit_un_PBM[5],
    )
end
fits_PBM = Optim.minimizer.([opt12_PBM]) 
report = Table(map(namefitparam_PBM, fits_PBM, ["1&2"]))

safesave(datadir("exp_pro", "fit12.html"), pretty_table(HTML, report, header = ["Fit", "N", "S", "C", "R0", "k0", "Tact", "θ"]))
pretty_table(report, header = ["Fit", "N", "S", "C", "R0", "k0", "Tact", "θ"])

# ---------------------- Compare fit to other experiments

pop3 = deepcopy(Rp12_params1)
@reset pop3.pch = pch3
@reset pop3.Tsh = Tsh3
sol3_Rp12 = solve(ODEProblem(pop3), Rodas3());
po_PBM3 = deepcopy(PBM1_params12)
@reset po_PBM3.pch = pch3
@reset po_PBM3.Tsh = Tsh3
sol3_PBM12 = solve(ODEProblem(po_PBM3), Rodas3());

pop5 = deepcopy(Rp12_params1)
@reset pop5.pch = pch5
@reset pop5.Tsh = Tsh5
sol5_Rp12 = solve(ODEProblem(pop5), Rodas3());
po_PBM5 = deepcopy(PBM1_params12)
@reset po_PBM5.pch = pch5
@reset po_PBM5.Tsh = Tsh5
sol5_PBM12 = solve(ODEProblem(po_PBM5), Rodas3());

pop8 = deepcopy(Rp12_params1)
@reset pop8.pch = pch8
@reset pop8.Tsh = Tsh8
sol8_Rp12 = solve(ODEProblem(pop8), Rodas3());
po_PBM8 = deepcopy(PBM1_params12)
@reset po_PBM8.pch = pch8
@reset po_PBM8.Tsh = Tsh8
sol8_PBM12 = solve(ODEProblem(po_PBM8), Rodas3());

begin
pl3 = blankplot_hrC(legend=:bottomright)
@df dat3 exptfplot!(:t, :T1, :T2, :T3, markershape=:circle)
modconvtplot!(sol3_PBM12, sol3_Rp12, labels=["PBM 1&2, 3", "3PP 1&2, 3"], markershape=:circle)
plot!(ylim=(-36, -31))
pl5 = blankplot_hrC(legend=:bottomright)
@df dat5 exptfplot!(:t, :T1, :T2, :T3, markershape=:square)
modconvtplot!(sol5_PBM12, sol5_Rp12, labels=["PBM 1&2, 5", "3PP 1&2, 5"], markershape=:square)
plot!(xlim=(0, 55), ylim=(-36, -31))
pl8 = blankplot_hrC(legend=:bottomright)
@df dat8 exptfplot!(:t, :T1, :T2, :T3, markershape=:utriangle)
modconvtplot!(sol8_PBM12, labels="PBM 1&2, 8", markershape=:utriangle)
plot!(xlim=(0, 25), ylim=(-33, -27))
plot(pl12, pl3, pl5, pl8, legend=:bottomright, link=:none, size=(1000, 600))
end
savefig(plotsdir("suc358_T.svg"))

