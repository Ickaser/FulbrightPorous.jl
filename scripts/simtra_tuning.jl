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

fitdat1 = PrimaryDryFit(collect(dat1.t), 
    (collect(filter(!isnan, dat1.T1)),
    collect(filter(!isnan, dat1.T2[dat1.t .< 60u"hr"])),
    collect(filter(!isnan, dat1.T3)),))

# plot(fitdat1)
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
u0_pikal = ustrip.([u"cm", u"K"], [hf0, Tsh1(0u"s")])
tspan = (0.0, 200.0)
prob = ODEProblem(lyo_1d_dae_f, u0_pikal, tspan, po_pikal; callback=end_drying_callback)
tsol = solve(prob, Rodas3(autodiff=false))

plot(u"hr", u"°C", xlim=(0, 60), ylim=(-40, -30))
modconvtplot!(tsol)
@df dat1 exptfplot!(:t, :T1, :T2, :T3)

function gen_sol_Rp(a123, po, u0, tspan)
    Rp_un = [u"cm^2*Torr*hr/g", u"cm*Torr*hr/g", u"1/cm"]
    new_params = @set po.Rp = RpFormFit((a123.*Rp_un)...)
    prob = ODEProblem(lyo_1d_dae_f, u0, tspan, new_params; callback=end_drying_callback)
    sol = solve(prob, Rodas3(autodiff=false))
    return sol, new_params
end
obj_Rp1 = x->obj_expT(gen_sol_Rp(x, po_pikal, u0_pikal, tspan)[1], fitdat1)

opt_Rp1 = optimize(obj_Rp1, [0.8, 14, 1])

sol_Rp1, Rp1_params = gen_sol_Rp(Optim.minimizer(opt_Rp1), po_pikal, u0_pikal, tspan);
plot(u"hr", u"°C", xlim=(0, 60), ylim=(-40, -30))
modconvtplot!(sol_Rp1)
@df dat1 exptfplot!(:t, :T1, :T2, :T3)

# -------------- Provide guess values for new parameters

# C = 2e9u"m"
# typical diameter: 50 μm
d_typical = 50u"μm"
V_typical = 4/3*π*d_typical^3
N_typical = fillvol/V_typical
# typical total volume: 
N = RampedVariable(N_typical / hf0)
# S = RampedVariable(V_typical^2*N_typical/hf0)
# C = 2e9u"m"
S = RampedVariable(1e0u"cm^5")
C = 2e2u"m"
R0 = 0.01u"cm^2*Torr*hr/g"
k0 = 1e-7u"1/s"
E = 10u"kJ/mol"

μ0_0 = 0.0 
μ2_0 = (0.0u"μm")^6
u0_mc = ustrip.([u"cm", u"K", NoUnits, u"cm^6"], [hf0, Tsh1(0u"s"), μ0_0, μ2_0])

pb_mc = (
    (Rp_g, hf0, c_solid, ρ_solution),
    (K_shf, A_v, A_p),
    (pch1, Tsh1),
    (k0, E, C, N, S, R0)
)
po_mc = ParamObjMC(pb_mc)

fit_un = [u"cm^-1", u"cm^5", u"s^-1", u"J/mol", u"m", u"m/s"]
fitguess = ustrip.(fit_un, [N(0), S(0), k0, E, C, R0])
# fitguess = log.(ustrip.(fit_un, [N(0), S(0), k0, E, C, R0]))


function gen_sol_constNS(fitprm, po, u0, tspan)
    new_params = deepcopy(po)
    # fitprm = exp.(fitprm) # Constrain all to be positive by taking exp
    @reset new_params.N = RampedVariable(fitprm[1]*fit_un[1])
    @reset new_params.S = RampedVariable(fitprm[2]*fit_un[2])
    @reset new_params.k0 = fitprm[3]*fit_un[3] # Constrain k0 to be positive
    @reset new_params.E  = fitprm[4]*fit_un[4]
    @reset new_params.C  = fitprm[5]*fit_un[5]
    @reset new_params.R0 = fitprm[6]*fit_un[6]
    prob = ODEProblem(lyo_mc_dae_f, u0, tspan, new_params; callback=end_drying_callback)
    sol = solve(prob, Rodas3(autodiff=false))
    return sol, new_params
end

tsol, tparams = gen_sol_constNS(fitguess, po_mc, u0_mc, tspan); 
1
modconvtplot(tsol)
plot!(fitdat1)
plot!(Tsh1, tmax=tsol.t[end]*u"hr")
# tparams.N(0) ≈ po_mc.N(0)

obj_MC1 = x->(obj_expT(gen_sol_constNS(x, po_mc, u0_mc, tspan)[1], fitdat1))
@time obj_MC1(fitguess)

opt_MC1 = optimize(obj_MC1, fitguess, Optim.Options(time_limit=5, show_trace=true))

Optim.minimizer(opt_MC1)

sol_MC1, MC1_params = gen_sol_constNS(Optim.minimizer(opt_MC1), po_mc, u0_mc, tspan);
plot(u"hr", u"°C", xlim=(0, 60), ylim=(-40, -30), legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3)
modconvtplot!(sol_MC1, sol_Rp1, labels=["MC" "Pikal"])


# -------------
# Everything from here on out needs to be reviewd
# Especially make sure simulations are all getting right parameters

# --------------------------
# Introduce Cycle 2
# Labeled as case 6 in table 8


dat2raw = CSV.read(datadir("exp_raw", "Simtra_Rp", "Suc_Cycle2.csv"), Table)
dat2 = map(simtra_rename, dat2raw)
filter!(x->x.t .> 0u"hr", dat2)

@df dat1 exptfplot(:t, :T1, :T2, :T3)
@df dat2 exptfplot!(:t, :T1, :T2, :T3, c=:red)
plot!(ylim=(-40, -30))

fitdat2 = PrimaryDryFit(collect(dat2.t), 
    (collect(filter(!isnan, dat2.T1)),
    collect(filter(!isnan, dat2.T2[dat2.t .< 85u"hr"])),
    collect(filter(!isnan, dat2.T3)),))

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

prob2p = ODEProblem(lyo_1d_dae_f, u0_pikal, tspan, pop2; callback=end_drying_callback)
sol2p = solve(prob2p, Rodas3(autodiff=false))

prob2mc = ODEProblem(lyo_mc_dae_f, u0_mc, tspan, pomc2; callback=end_drying_callback)
sol2mc = solve(prob2mc, Rodas3(autodiff=false))


plot(u"hr", u"°C", xlim=(0, 150), ylim=(-36, -31), legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3, marker=:square)
@df dat2 exptfplot!(:t, :T1, :T2, :T3, marker=:circle)
modconvtplot!(sol_MC1, sol_Rp1, labels=["MC" "Pikal"],marker=:square)
modconvtplot!(sol2mc, sol2p, labels=["MC" "Pikal"],marker=:circle)


# ------ Fit to Cycle 2 data separately

# GitHub Copilot guessed the next few lines with very few corrections!
obj2p = x->obj_expT(gen_sol_Rp(x, pop2, u0_pikal, tspan)[1], fitdat2)
opt2p = optimize(obj2p, Optim.minimizer(opt_Rp1), Optim.Options(time_limit=5))
obj2mc = x->obj_expT(gen_sol_constNS(x, pomc2, u0_mc, tspan)[1], fitdat2)
opt2mc = optimize(obj2mc, Optim.minimizer(opt_MC1), Optim.Options(g_tol = 1e-3, time_limit=10))

sol_Rp2, Rp2_params = gen_sol_Rp(Optim.minimizer(opt2p), pop2, u0_pikal, tspan);
sol_MC2, MC2_params = gen_sol_constNS(Optim.minimizer(opt2mc), pomc2, u0_mc, tspan);

# ----------- Plot separately
begin
plot(u"hr", u"°C", xlim=(0, 70), ylim=(-36, -31), legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3, marker=:square)
modconvtplot!(sol_MC1, sol_Rp1, labels=["MC" "Pikal"],marker=:square)
savefig(plotsdir("suc1_T.svg"))
end

begin
plot(u"hr", u"°C", xlim=(0, 125), ylim=(-36, -32), legend=:bottomright)
@df dat2 exptfplot!(:t, :T1, :T2, :T3, marker=:circle)
modconvtplot!(sol_MC2, sol_Rp2, labels=["MC" "Pikal"],marker=:circle)
savefig(plotsdir("suc2_T.svg"))
end

# --------------- Fit to Cycle 1 and 2 simultaneously

function gen2sol_constNS(fitprm, po1, po2, u0, tspan)
    sol1, prm1 = gen_sol_constNS(fitprm, po1, u0, tspan)
    sol2, prm2 = gen_sol_constNS(fitprm, po2, u0, tspan)
    return sol1, sol2, prm1, prm2
end
function obj12mc(fitrpm)
    sol1, sol2, prm1, prm2 = gen2sol_constNS(fitrpm, po_mc, pomc2, u0_mc, tspan)
    obj1 = obj_expT(sol1, fitdat1)
    obj2 = obj_expT(sol2, fitdat2)
    return obj1 + obj2
end
g12mc = sqrt.(abs.((Optim.minimizer(opt2mc) .* Optim.minimizer(opt_MC1))))
opt12mc = optimize(obj12mc, g12mc, Optim.Options(g_tol = 1e-2, time_limit=10))


sol2_MC1, sol2_MC2, MC1_params2, MC2_params2 = gen2sol_constNS(Optim.minimizer(opt12mc), po_mc, pomc2, u0_mc, tspan);

begin
plot(u"hr", u"°C", xlim=(0, 120), ylim=(-36, -31), legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3, markershape=:circle)
@df dat2 exptfplot!(:t, :T1, :T2, :T3, markershape=:square)
modconvtplot!(sol2_MC1, labels="MC", markershape=:circle)
modconvtplot!(sol2_MC2, labels="MC", markershape=:square)
end
savefig(plotsdir("suc12_T.svg"))

# -----------------
# Inspect fit results in Rp space rather than temperature space
minteg = cumsum(dat2.mfr[1:end-1].*(dat2.t[2:end]-dat2.t[1:end-1]))
errfac = minteg[end] / (fillvol*(ρ_solution-c_solid))
ld = @. errfac*minteg/(A_p*(ρ_solution-c_solid)) |> u"cm"

l = hf0 .- sol_MC2[1,:]*u"cm"
Rp_mc2 = [FulbrightPores.lyo_mc_dae!(fill(0.0, 4), u, MC2_params, t) for (u,t) in zip(sol_MC2.u, sol_MC2.t)]
Rp_p2 = Rp2_params.Rp.(l)
begin
plot(u"cm", u"cm^2*Torr*hr/g", ylims=(0, 20), xlabel=L"h_d", ylabel=L"R_p", unitformat=(x,y)->"$x [$(latexify(y))]")
plot!(l, Rp_mc2, label="PBE model")
plot!(l, Rp_p2, label="3-parameter Pikal")
plot!(ld, dat2.Rp, label="TDLAS")
savefig(plotsdir("suc2_Rp.svg"))
end


fits = Optim.minimizer.([opt_MC1, opt2mc, opt12mc]) 
function namefitparam(fit, label)
    nt = (label = label,
        N = fit[1]*fit_un[1],
        S = fit[2]*fit_un[2],
        k0 = fit[3]*fit_un[3],
        E = fit[4]*fit_un[4],
        C = fit[5]*fit_un[5],
        R0 = fit[6]*fit_un[6])
end
report = Table(map(namefitparam, fits, ["Cycle 1", "Cycle 2", "Both"]))

safesave(datadir("exp_pro", "fit12.html"), pretty_table(HTML, report))