using DrWatson
@quickactivate :FulbrightPores
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
K_shf = fitKv(Optim.minimizer(optKv))[2]

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


po_pikal = ParamObjPikal(pb_pikal)
# u0_pikal = ustrip.([u"cm", u"K"], [hf0, Tsh1(0u"s")])
# tspan = (0.0, 200.0)
# prob = ODEProblem(lyo_1d_dae_f, u0_pikal, tspan, po_pikal; callback=end_drying_callback)
prob = ODEProblem(po_pikal)
tsol = solve(prob, Rodas3(autodiff=false))

plot(u"hr", u"°C", xlim=(0, 60), ylim=(-40, -30))
modconvtplot!(tsol)
@df dat1 exptfplot!(:t, :T1, :T2, :T3)

function gen_sol_Rp(a123, po)
    Rp_un = [u"cm^2*Torr*hr/g", u"cm*Torr*hr/g", u"1/cm"]
    new_params = @set po.Rp = RpFormFit((a123.*Rp_un)...)
    prob = ODEProblem(new_params)
    sol = solve(prob, Rodas3(autodiff=false))
    return sol, new_params
end
obj_Rp1 = x->obj_expT(gen_sol_Rp(x, po_pikal)[1], fitdat1)
@time obj_Rp1([0.8, 14, 1])

opt_Rp1 = optimize(obj_Rp1, [0.8, 14, 1], Optim.Options(g_tol=1e-3))

sol_Rp1, Rp1_params = gen_sol_Rp(Optim.minimizer(opt_Rp1), po_pikal);
plot(u"hr", u"°C", xlim=(0, 60), ylim=(-40, -30))
modconvtplot!(sol_Rp1)
@df dat1 exptfplot!(:t, :T1, :T2, :T3)

# -------------- Provide guess values for new parameters

# C = 2e9u"m"
# typical diameter: 50 μm
d_typical = 50u"μm"
# typical total volume: 
V_typical = 4/3*π*d_typical^3
N_typical = fillvol/V_typical
 
begin
# Good guesses for Knudsen
N = RampedVariable(N_typical / hf0)
S = RampedVariable(V_typical^2*N_typical/hf0)
C = 2e-2u"m"
R0 = 0.05u"cm^2*Torr*hr/g"

# Good guesses for viscous
# N = RampedVariable(N_typical / hf0) 
# S = RampedVariable(1e0u"cm^5")
# C = 2e2u"m"
# R0 = 0.01u"cm^2*Torr*hr/g"

k0 = 1e-11u"1/s"
E = u"R"*Tsh1(0u"s")

pb_mc = (
    (Rp_g, hf0, c_solid, ρ_solution),
    (K_shf, A_v, A_p),
    (pch1, Tsh1),
    (k0, E, C, N, S, R0)
)
po_mc = ParamObjMC(pb_mc)

tprob = ODEProblem(po_mc)
tsol = solve(tprob, Rodas3(autodiff=true))
poreevolplot(tsol)
end

fit_un = [u"cm^-1", u"cm^5", u"s^-1", u"J/mol", u"m", u"m/s"]
fitguess = ustrip.(fit_un, [N(0), S(0), k0, E, C, R0])
fitguesslog = log.(ustrip.(fit_un, [N(0), S(0), k0, E, C, R0]))


function gen_sol_constNS(fitprmlog, po; verbose=false)
    new_params = deepcopy(po)
    fitprm = exp.(fitprmlog) # Constrain all to be positive by taking exp
    verbose && @info "fitprm" fitprm
    @reset new_params.N = RampedVariable(fitprm[1]*fit_un[1])
    @reset new_params.S = RampedVariable(fitprm[2]*fit_un[2])
    @reset new_params.k0 = fitprm[3]*fit_un[3] 
    @reset new_params.E  = fitprm[4]*fit_un[4]
    @reset new_params.C  = fitprm[5]*fit_un[5]
    @reset new_params.R0 = fitprm[6]*fit_un[6]
    prob = ODEProblem(new_params)
    sol = solve(prob, Rodas3(autodiff=true))
    return sol, new_params
end

tsol, tparams = gen_sol_constNS(fitguesslog, po_mc); 
modconvtplot(tsol)
plot!(fitdat1)

obj_MC1 = x->(obj_expT(gen_sol_constNS(x, po_mc; verbose=true)[1], fitdat1))
@time obj_MC1(fitguesslog)

opt_MC1 = optimize(obj_MC1, fitguesslog, NelderMead(), Optim.Options(g_tol=1e-2))
# opt_MC1 = optimize(obj_MC1, fitguesslog, LBFGS(), Optim.Options(g_tol = 1e-2))

exp.(Optim.minimizer(opt_MC1)) .*fit_un

sol_MC1, MC1_params = gen_sol_constNS(Optim.minimizer(opt_MC1), po_mc);
begin
plot(u"hr", u"°C", xlim=(0, 60), ylim=(-40, -30), legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3)
modconvtplot!(sol_MC1, sol_Rp1, labels=["MC" "Pikal"])
end

poreevolplot(sol_MC1)


# -------------
# Everything from here on out needs to be reviewd
# Especially make sure simulations are all getting right parameters

# --------------------------
# Introduce Cycle 2
# Labeled as case 6 in table 8


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
@df dat2 plot(:t, :pch_cm, ylim=(0, 200))
plot!(pch2)
Tsh2 = Tsh1 - 5u"K"
@df dat2 plot(:t, :Tsh, xlim=(0, 2))
plot!(Tsh2)
pop2 = deepcopy(Rp1_params)
@reset pop2.pch = pch2
@reset pop2.Tsh = Tsh2
pomc2 = deepcopy(MC1_params)
@reset pomc2.pch = pch2
@reset pomc2.Tsh = Tsh2

# ------ Compare existing fit to Cycle 2 data

prob2p = ODEProblem(pop2)
sol2p = solve(prob2p, Rodas3(autodiff=false))

prob2mc = ODEProblem(pomc2)
sol2mc = solve(prob2mc, Rodas3(autodiff=false))


begin
plot(u"hr", u"°C", xlim=(0, 150), ylim=(-36, -31), legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3, marker=:square)
@df dat2 exptfplot!(:t, :T1, :T2, :T3, marker=:circle)
modconvtplot!(sol_MC1, sol_Rp1, labels=["MC" "Pikal"],marker=:square)
modconvtplot!(sol2mc, sol2p, labels=["MC" "Pikal"],marker=:circle)
end


# ------ Fit to Cycle 2 data separately

# GitHub Copilot guessed the next few lines with very few corrections!
obj2p = x->obj_expT(gen_sol_Rp(x, pop2)[1], fitdat2)
opt2p = optimize(obj2p, Optim.minimizer(opt_Rp1), NelderMead(), Optim.Options(g_tol = 1e-3, time_limit=15))
obj2mc = x->obj_expT(gen_sol_constNS(x, pomc2)[1], fitdat2)
opt2mc = optimize(obj2mc, Optim.minimizer(opt_MC1), NelderMead(), Optim.Options(g_tol = 1e-3, time_limit=10))

sol_Rp2, Rp2_params = gen_sol_Rp(Optim.minimizer(opt2p), pop2);
sol_MC2, MC2_params = gen_sol_constNS(Optim.minimizer(opt2mc), pomc2);

# ----------- Plot separately
begin
plot(u"hr", u"°C", xlim=(0, 70), ylim=(-36, -31), legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3, marker=:square)
modconvtplot!(sol_MC1, sol_Rp1, labels=["MC" "Pikal"],marker=:square)
end
savefig(plotsdir("suc1_T.svg"))

begin
plot(u"hr", u"°C", xlim=(0, 125), ylim=(-36, -32), legend=:bottomright)
@df dat2 exptfplot!(:t, :T1, :T2, :T3, marker=:circle)
modconvtplot!(sol_MC2, sol_Rp2, labels=["MC" "Pikal"],marker=:circle)
end
savefig(plotsdir("suc2_T.svg"))

# --------------- Fit to Cycle 1 and 2 simultaneously

function gen2sol_Rp(fitprm, po1, po2)
    sol1, prm1 = gen_sol_Rp(fitprm, po1)
    sol2, prm2 = gen_sol_Rp(fitprm, po2)
    return sol1, sol2, prm1, prm2
end
function obj12p(fitrpm)
    sol1, sol2, prm1, prm2 = gen2sol_Rp(fitrpm, po_pikal, pop2)
    obj1 = obj_expT(sol1, fitdat1)
    obj2 = obj_expT(sol2, fitdat2)
    return obj1 + obj2
end
function gen2sol_constNS(fitprm, po1, po2)
    sol1, prm1 = gen_sol_constNS(fitprm, po1)
    sol2, prm2 = gen_sol_constNS(fitprm, po2)
    return sol1, sol2, prm1, prm2
end
function obj12mc(fitprm; verbose=false)
    sol1, sol2, prm1, prm2 = gen2sol_constNS(fitprm, po_mc, pomc2)
    obj1 = obj_expT(sol1, fitdat1)
    obj2 = obj_expT(sol2, fitdat2)
    verbose && @info "fit" exp.(fitprm) obj1 obj2
    return obj1 + obj2
end

# g12mc = sqrt.(abs.((Optim.minimizer(opt2mc) .* Optim.minimizer(opt_MC1))))
# g12mc = (Optim.minimizer(opt2mc) .+ Optim.minimizer(opt_MC1))/2
g12mc = Optim.minimizer(opt2mc)
opt12mc = optimize(obj12mc, g12mc, NelderMead(), Optim.Options(g_tol = 1e-3))

g12p = sqrt.(abs.((Optim.minimizer(opt2p) .* Optim.minimizer(opt_Rp1))))
opt12p = optimize(obj12p, g12p, Optim.Options(g_tol = 1e-2, time_limit=10))


sol2_MC1, sol2_MC2, MC1_params2, MC2_params2 = gen2sol_constNS(Optim.minimizer(opt12mc), po_mc, pomc2);
sol2_Rp1, sol2_Rp2, Rp1_params2, Rp2_params2 = gen2sol_Rp(Optim.minimizer(opt12p), po_pikal, pop2);

begin
pl1 = plot(u"hr", u"°C", xlim=(0, 120), ylim=(-36, -31), legend=:bottomright)
# @df dat1 exptfplot!(:t, :T1, :T2, :T3, markershape=:circle)
# @df dat2 exptfplot!(:t, :T1, :T2, :T3, markershape=:square)
plot!(fitdat1)
plot!(fitdat2)
modconvtplot!(sol2_MC1, sol2_Rp1, labels=["MC1" "Pikal"], markershape=:circle)
modconvtplot!(sol2_MC2, sol2_Rp2, labels=["MC2" "Pikal"], markershape=:square)
end
savefig(plotsdir("suc12_T.svg"))

l = range(0u"cm", hf0, length=50)
Rp1 = Rp1_params.Rp.(l)
Rp2 = Rp2_params.Rp.(l)
Rp12 = Rp2_params2.Rp.(l)

begin
colors = [:blue :red :purple :brown]
rpattr = (marker=:utriangle, seriestype=:samplemarkers, linestyle=:dash)
pl2 = poreevolplot(sol_MC1, label="Fit 1, Cycle 1", c=colors[1], marker=:circle)
poreevolplot!(sol2mc, label="Fit 1, Cycle 2", c=colors[1], marker=:square)
poreevolplot!(sol_MC2, label="Fit 2, Cycle 2", c=colors[2], marker=:square)
poreevolplot!(sol2_MC1, label="Fit both, Cycle 1", c=colors[3], marker=:circle)
poreevolplot!(sol2_MC2, label="Fit both, Cycle 2", c=colors[3], marker=:square)
plot!(pl2, l, Rp1, subplot=4, label="Fit 1, 3-parameter Pikal"; c=colors[1], rpattr...)
plot!(pl2, l, Rp2, subplot=4, label="Fit 2, 3-parameter Pikal"; c=colors[2], rpattr...)
plot!(pl2, l, Rp12, subplot=4, label="Fit both, 3-parameter Pikal"; c=colors[3], rpattr...)
plot!(size=(800, 400))
end
savefig(plotsdir("suc12_Rp.svg"))
plot(pl1, pl2, size=(1000, 400))

# -----------------
# Inspect fit results in Rp space rather than temperature space
minteg = cumsum(dat2.mfr[1:end-1].*(dat2.t[2:end]-dat2.t[1:end-1]))
insert!(minteg, 1, 0u"g")
errfac = minteg[end] / (fillvol*(ρ_solution-c_solid))
ld2 = @. errfac*minteg/(A_p*(ρ_solution-c_solid)) |> u"cm"

l2 = hf0 .- sol_MC2[1,:]*u"cm"
Rp_p2 = Rp2_params.Rp.(l)
Rp_mc2 = [FulbrightPores.lyo_mc_dae!(fill(0.0, 4), u, MC2_params, t) for (u,t) in zip(sol_MC2.u, sol_MC2.t)]
Rp_p2 = Rp2_params.Rp.(l)
begin
plot(u"cm", u"cm^2*Torr*hr/g", ylims=(0, 20), xlabel=L"h_d", ylabel=L"R_p", unitformat=(x,y)->"$x [$(latexify(y))]")
plot!(l, Rp_mc2, label="PBE model")
plot!(pl2, l, Rp_p2, subplot=4, label="3-parameter Pikal")
plot!(ld, dat2.Rp, label="TDLAS")
end
savefig(plotsdir("suc2_Rp.svg"))


fits = Optim.minimizer.([opt_MC1, opt2mc, opt12mc]) 
function namefitparam(fitlog, label)
    fit = exp.(fitlog)
    nt = (label = label,
        N = fit[1]*fit_un[1],
        S = fit[2]*fit_un[2],
        k0 = fit[3]*fit_un[3],
        E = fit[4]*fit_un[4],
        C = fit[5]*fit_un[5],
        R0 = fit[6]*fit_un[6])
end
report = Table(map(namefitparam, fits, ["Cycle 1", "Cycle 2", "Both"]))

safesave(datadir("exp_pro", "fit12.html"), pretty_table(HTML, report, header = ["Fit", "N", "S", "k0", "E", "C", "R0"]))
pretty_table(report, header = ["Fit", "N", "S", "k0", "E", "C", "R0"])