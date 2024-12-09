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
@df dat2 plot(:t, :pch_cm, ylim=(0, 200))
plot!(pch2)
Tsh2 = Tsh1 - 5u"K"
@df dat2 plot(:t, :Tsh, xlim=(0, 2))
plot!(Tsh2)


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

function gen_sol_Rp(a123, po)
    Rp_un = [u"cm^2*Torr*hr/g", u"cm*Torr*hr/g", u"1/cm"]
    new_params = @set po.Rp = RpFormFit((a123.*Rp_un)...)
    prob = ODEProblem(new_params)
    sol = solve(prob, Rodas3(autodiff=false))
    return sol, new_params
end
obj_Rp1 = x->obj_expT(gen_sol_Rp(x, pop1)[1], fitdat1)
@time obj_Rp1([0.8, 14, 1])

opt_Rp1 = optimize(obj_Rp1, [0.8, 14, 1], Optim.Options(g_tol=1e-3))

sol1_Rp1, Rp1_params = gen_sol_Rp(Optim.minimizer(opt_Rp1), pop1);
plot(u"hr", u"°C", xlim=(0, 60), ylim=(-40, -30))
modconvtplot!(sol1_Rp1)
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

# tprob = ODEProblem(pomc1)
# tsol = solve(tprob, Rodas3(autodiff=true))
# poreevolplot(tsol)
end

fit_un_arr = [u"cm^-1", u"cm^5", u"s^-1", u"J/mol", u"m", u"m/s"]
fitguess = ustrip.(fit_un_arr, [N(0), S(0), k0, E, C, R0])
fitguesslog = log.(fitguess)
function gen_sol_constNS_arr(fitprmlog, po; verbose=false)
    new_params = deepcopy(po)
    fitprm = exp.(fitprmlog) # Constrain all to be positive by taking exp
    verbose && @info "fitprm" fitprm
    @reset new_params.N = RampedVariable(fitprm[1]*fit_un_arr[1])
    @reset new_params.S = RampedVariable(fitprm[2]*fit_un_arr[2])
    # Arrhenius fit
    k0 = fitprm[3]*fit_un_arr[3]
    E = fitprm[4]*fit_un_arr[4]
    @reset new_params.k_f = T->k0*exp(-E/(R*T))
    @reset new_params.C  = fitprm[5]*fit_un_arr[5]
    @reset new_params.R0 = fitprm[6]*fit_un_arr[6]
    prob = ODEProblem(new_params)
    sol = solve(prob, Rodas3(autodiff=true))
    return sol, new_params
end

fit_un_arb = [u"cm^-1", u"cm^5", fill(u"s^-1", 4)..., u"m", u"m/s"]
fitguess = ustrip.(fit_un_arb, [N(0), S(0), k0/10, k0/5, k0/2, k0, C, R0])
fitguesslog = log.(fitguess)
Tstops = [-45.0, -37, -33, -30]u"°C" .|> u"K"
function gen_sol_constNS_arb(fitprmlog, po; verbose=false)
    new_params = deepcopy(po)
    fitprm = exp.(fitprmlog) # Constrain all to be positive by taking exp
    verbose && @info "fitprm" fitprm
    @reset new_params.N = RampedVariable(fitprm[1]*fit_un_arb[1])
    @reset new_params.S = RampedVariable(fitprm[2]*fit_un_arb[2])
    # Linear interpolant
    # kstops = fitprm[3:6]*fit_un_arb[3]
    # spl = LinearInterpolation(kstops, Tstops, extrapolate=true)
    # Lin interp in log space
    klogs = fitprmlog[3:6]
    spll = LinearInterpolation(klogs, Tstops, extrapolate=true)
    spl = T->exp(spll(T))*fit_un_arb[3]
    # PCHIP interp
    # spl = PCHIPInterpolation(fitprm[3:6]*fit_un_arb[3], Tstops, extrapolate = true)
    @reset new_params.k_f = spl
    @reset new_params.C  = fitprm[7]*fit_un_arb[7]
    @reset new_params.R0 = fitprm[8]*fit_un_arb[8]
    prob = ODEProblem(new_params)
    sol = solve(prob, Rodas3(autodiff=true))
    return sol, new_params
end

# tsol, tparams = gen_sol_constNS_arb(fitguesslog, pomc1); 
# modconvtplot(tsol)
# plot!(fitdat1)

obj_MC1 = x->(issorted(x[3:6]) ? (obj_expT(gen_sol_constNS_arb(x, pomc1)[1], fitdat1)) : NaN)
@time obj_MC1(fitguesslog)

opt_MC1 = optimize(obj_MC1, fitguesslog, NelderMead(), Optim.Options(g_tol=1e-2))
# opt_MC1 = optimize(obj_MC1, fitguesslog, LBFGS(), Optim.Options(g_tol = 1e-2))

exp.(Optim.minimizer(opt_MC1)) .*fit_un_arb

sol1_MC1, MC1_params = gen_sol_constNS_arb(Optim.minimizer(opt_MC1), pomc1);

# ------ Compare existing fit to Cycle 2 data

# ------ Fit to Cycle 2 data separately

obj2p = x->obj_expT(gen_sol_Rp(x, pop2)[1], fitdat2)
opt2p = optimize(obj2p, Optim.minimizer(opt_Rp1), NelderMead(), Optim.Options(g_tol = 1e-3, time_limit=15))
# Enforce monontonic k in objective
obj2mc = x->(issorted(x[3:6]) ? obj_expT(gen_sol_constNS_arb(x, pomc2)[1], fitdat2) : NaN)
opt_MC2 = optimize(obj2mc, Optim.minimizer(opt_MC1), NelderMead(), Optim.Options(g_tol = 1e-3, time_limit=10))

sol2_Rp2, Rp2_params = gen_sol_Rp(Optim.minimizer(opt2p), pop2);
sol2_MC2, MC2_params = gen_sol_constNS_arb(Optim.minimizer(opt_MC2), pomc2);

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

# Cross-check fits
sol1_MC2 = gen_sol_constNS_arb(Optim.minimizer(opt_MC2), pomc1)[1]
sol2_MC1 = gen_sol_constNS_arb(Optim.minimizer(opt_MC1), pomc2)[1]

# --------------- Fit to Cycle 1 and 2 simultaneously

function gen2sol_Rp(fitprm, po1, po2)
    sol1, prm1 = gen_sol_Rp(fitprm, po1)
    sol2, prm2 = gen_sol_Rp(fitprm, po2)
    return sol1, sol2, prm1, prm2
end
function obj12p(fitrpm)
    sol1, sol2, prm1, prm2 = gen2sol_Rp(fitrpm, pop1, pop2)
    obj1 = obj_expT(sol1, fitdat1)
    obj2 = obj_expT(sol2, fitdat2)
    return obj1 + obj2
end
function gen2sol_constNS(fitprm, po1, po2)
    sol1, prm1 = gen_sol_constNS_arb(fitprm, po1)
    sol2, prm2 = gen_sol_constNS_arb(fitprm, po2)
    return sol1, sol2, prm1, prm2
end
function obj12mc(fitprm; verbose=false)
    if !issorted(fitprm[3:6])
        return NaN
    end
    sol1, sol2, prm1, prm2 = gen2sol_constNS(fitprm, pomc1, pomc2)
    obj1 = obj_expT(sol1, fitdat1)
    obj2 = obj_expT(sol2, fitdat2)
    verbose && @info "fit" exp.(fitprm) obj1 obj2
    return obj1 + obj2
end

g12mc = (Optim.minimizer(opt_MC2) .+ Optim.minimizer(opt_MC1))/2
# g12mc = Optim.minimizer(opt_MC2)
opt12mc = optimize(obj12mc, g12mc, NelderMead(), Optim.Options(g_tol = 1e-3))

g12p = sqrt.(abs.((Optim.minimizer(opt2p) .* Optim.minimizer(opt_Rp1))))
opt12p = optimize(obj12p, g12p, Optim.Options(g_tol = 1e-2, time_limit=10))


sol1_MC12, sol2_MC12, MC1_params2, MC2_params2 = gen2sol_constNS(Optim.minimizer(opt12mc), pomc1, pomc2);
sol1_Rp12, sol2_Rp12, Rp1_params2, Rp2_params2 = gen2sol_Rp(Optim.minimizer(opt12p), pop1, pop2);

begin
pl1 = plot(u"hr", u"°C", xlim=(0, 120), ylim=(-36, -31), legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3, markershape=:circle)
@df dat2 exptfplot!(:t, :T1, :T2, :T3, markershape=:square)
# plot!(fitdat1)
# plot!(fitdat2)
modconvtplot!(sol1_MC12, sol1_Rp12, labels=["MC" "Pikal"], markershape=:circle)
modconvtplot!(sol2_MC12, sol2_Rp12, labels=["MC" "Pikal"], markershape=:square)
end
savefig(plotsdir("suc12_T.svg"))

klogs = Optim.minimizer(opt12mc)[3:6]
kstops = exp.(klogs)*fit_un_arb[3]

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
    sum(abs2.(logks .- klogs))
end
p0_wlf = [1e-16, 273-85, 9, 10]
obj_wlf(p0_wlf)
opt_wlf = optimize(obj_wlf, p0_wlf, NelderMead(), Optim.Options(show_trace=true))
plot(Tcheck,  k_wlf.(Tcheck, [Optim.minimizer(opt_wlf).*wlf_un]))
plot!(Tstops,  kstops, yscale=:log10)

# MYEGA temparature dependence of viscosity
function k_myega(T, p)
    k0, K, C = p
    # @info "check" K/T C/T K/T*exp(C/T)
    return k0*exp((K/T)*exp(C/T))
end
myega_un = [u"s^-1", u"K", u"K"]
function obj_myega(p)
    if p[1] < 0 #|| p[2] < 100 || p[3] < 100 
        return NaN
    end
    ks = k_myega.(Tstops, [p.*myega_un])
    # @info "checkobj" ks p
    sum(abs2.(log.(ustrip.(u"s^-1", ks)) - klogs))
end
p0_myega = [1e-10, -243, 2]
obj_myega(p0_myega)
opt_myega = optimize(obj_myega, p0_myega, NelderMead(), Optim.Options(show_trace=true))
plot(Tcheck,  k_myega.(Tcheck, [Optim.minimizer(opt_myega).*myega_un]))
plot!(Tstops,  kstops, yscale=:log10)

# playground temparature dependence of viscosity
function k_play(T, p)
    k0, K, C = p
    # x = ustrip(u"K", T-K)
    # x = ustrip(u"K", T-K)
    x = (T-K)/C
    sig = (tanh(x)+1)/2
    return k0*exp(sig*(T-K)/C)
end
play_un = [u"s^-1", u"K", u"K"]
p0_play = [1e-16, 237, 0.2]
scatter(Tcheck,  k_play.(Tcheck, [p0_play .* play_un]), yscale=:log10)
plot!(Tstops,  kstops, yscale=:log10)
function obj_play(p)
    if p[1] < 0 #|| p[2] < 100 || p[3] < 100 
        return NaN
    end
    ks = k_play.(Tstops, [p.*play_un])
    sum(abs2.(log.(ustrip.(u"s^-1", ks)) - klogs))
end
obj_play(p0_play)
opt_play = optimize(obj_play, p0_play, NelderMead(), Optim.Options(show_trace=true))
plot(Tcheck,  k_play.(Tcheck, [Optim.minimizer(opt_play).*play_un]))
plot!(Tstops,  kstops, yscale=:log10)
play_params1 = @set MC1_params2.k_f = T->k_play(T, Optim.minimizer(opt_play).*play_un)
play_params2 = @set MC2_params2.k_f = T->k_play(T, Optim.minimizer(opt_play).*play_un)
sol1_play = solve(ODEProblem(MC1_params2), Rodas3(autodiff=true))
sol2_play = solve(ODEProblem(MC2_params2), Rodas3(autodiff=true))
poreevolplot(sol1_play)
poreevolplot!(sol2_play)

fit_un_play = [u"cm^-1", u"cm^5", u"s^-1", u"K", u"K", u"m", u"m/s"]
function gen_sol_constNS_play(fitprmlog, po; verbose=false)
    new_params = deepcopy(po)
    fitprm = exp.(fitprmlog) # Constrain all to be positive by taking exp
    verbose && @info "fitprm" fitprm
    @reset new_params.N = RampedVariable(fitprm[1]*fit_un_play[1])
    @reset new_params.S = RampedVariable(fitprm[2]*fit_un_play[2])
    @reset new_params.k_f = T->k_play(T, fitprm[3:5].*fit_un_play[3:5])
    @reset new_params.C  = fitprm[6]*fit_un_play[6]
    @reset new_params.R0 = fitprm[7]*fit_un_play[7]
    prob = ODEProblem(new_params)
    sol = solve(prob, Rodas3(autodiff=true))
    return sol, new_params
end
function gen2sol_constNS_play(fitprm, po1, po2)
    sol1, prm1 = gen_sol_constNS_play(fitprm, po1)
    sol2, prm2 = gen_sol_constNS_play(fitprm, po2)
    return sol1, sol2, prm1, prm2
end
function obj12play(fitprm; verbose=true)
    sol1, sol2, prm1, prm2 = gen2sol_constNS_play(fitprm, MC1_params2, MC2_params2)
    obj1 = obj_expT(sol1, fitdat1)
    obj2 = obj_expT(sol2, fitdat2)
    verbose && @info "fit" exp.(fitprm) obj1 obj2
    return obj1 + obj2
end
p0_playfit = [6.4e7, 1.1e-6, 3e-16, 237, 0.3, 0.02, 2.5e4]
p0_playfitlog = log.(p0_playfit)

opt12play = optimize(obj12play, p0_playfitlog, NelderMead(), Optim.Options(g_tol = 1e-2, time_limit=10))
sol1_play, sol2_play, play1_params, play2_params = gen2sol_constNS_play(Optim.minimizer(opt12play), MC1_params2, MC2_params2);
poreevolplot(sol1_play)
poreevolplot!(sol2_play)

begin
pl1 = plot(u"hr", u"°C", xlim=(0, 120), ylim=(-36, -31), legend=:bottomright)
@df dat1 exptfplot!(:t, :T1, :T2, :T3, markershape=:circle)
@df dat2 exptfplot!(:t, :T1, :T2, :T3, markershape=:square)
# plot!(fitdat1)
# plot!(fitdat2)
modconvtplot!(sol1_play, sol1_Rp12, labels=["MC" "Pikal"], markershape=:circle)
modconvtplot!(sol2_play, sol2_Rp12, labels=["MC" "Pikal"], markershape=:square)
modconvtplot!(sol1_MC12)
modconvtplot!(sol1_MC12)
end
savefig(plotsdir("suc12_T.svg"))

modconvtplot(sol1_play)
modconvtplot!(sol2_play)
modconvtplot!(sol1_MC12)
modconvtplot!(sol1_MC12)

Tcheck = range(extrema(Tstops)..., length=50)
begin
plot(u"°C", u"s^-1", legend=:bottomright, unitformat=:square)
kineticplot!(MC1_params, T=Tcheck, c=:red, label="Fit to case 1")
kineticplot!(MC2_params, T=Tcheck, c=:blue, label="Fit to case 2")
kineticplot!(MC1_params2, T=Tcheck, c=:purple, label="Fit to both")
kineticplot!(play1_params, T=Tcheck, label="2-param sigmoid+exp")
end
savefig(plotsdir("kinetic.svg"))



l = range(0u"cm", hf0, length=50)
Rp1 = Rp1_params.Rp.(l)
Rp2 = Rp2_params.Rp.(l)
Rp12 = Rp2_params2.Rp.(l)

begin
colors = [:red :blue :purple :brown]
rpattr = (marker=:utriangle, seriestype=:samplemarkers, linestyle=:dash)
pl2 = poreevolplot(sol1_MC1, label="Fit 1, Cycle 1", c=colors[1], marker=:circle)
# poreevolplot!(sol2_MC1, label="Fit 1, Cycle 2", c=colors[1], marker=:square)
# poreevolplot!(sol1_MC2, label="Fit 2, Cycle 1", c=colors[2], marker=:circle)
poreevolplot!(sol2_MC2, label="Fit 2, Cycle 2", c=colors[2], marker=:square)
# poreevolplot!(sol1_MC12, label="Fit both, Cycle 1", c=colors[3], marker=:circle)
# poreevolplot!(sol2_MC12, label="Fit both, Cycle 2", c=colors[3], marker=:square)
poreevolplot!(sol1_play, label="Fit both (se), Cycle 1", c=colors[3], marker=:circle)
poreevolplot!(sol2_play, label="Fit both (se), Cycle 2", c=colors[3], marker=:square)
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
ld = @. errfac*minteg/(A_p*(ρ_solution-c_solid)) |> u"cm"

l2 = hf0 .- sol_MC2[1,:]*u"cm"
Rp_p2 = Rp2_params.Rp.(l)
Rp_mc2 = [FulbrightPores.lyo_mc_dae!(fill(0.0, 4), u, MC2_params, t) for (u,t) in zip(sol_MC2.u, sol_MC2.t)]
Rp_p2 = Rp2_params.Rp.(l)
begin
plot(u"cm", u"cm^2*Torr*hr/g", ylims=(0, 20), xlabel=L"h_d", ylabel=L"R_p", unitformat=(x,y)->"$x [$(latexify(y))]")
plot!(l2, Rp_mc2, label="PBE model")
plot!(l2, Rp_p2, label="3-parameter Pikal")
plot!(ld, dat2.Rp, label="TDLAS")
end
savefig(plotsdir("suc2_Rp.svg"))


fits = Optim.minimizer.([opt_MC1, opt_MC2, opt12mc]) 
function namefitparam_arr(fitlog, label)
    fit = exp.(fitlog)
    nt = (label = label,
        N = fit[1]*fit_un_arr[1],
        S = fit[2]*fit_un_arr[2],
        k0 = fit[3]*fit_un_arr[3],
        E = fit[4]*fit_un_arr[4],
        C = fit[5]*fit_un_arr[5],
        R0 = fit[6]*fit_un_arr[6])
end
report = Table(map(namefitparam_arr, fits, ["Cycle 1", "Cycle 2", "Both"]))

safesave(datadir("exp_pro", "fit12.html"), pretty_table(HTML, report, header = ["Fit", "N", "S", "k0", "E", "C", "R0"]))
pretty_table(report, header = ["Fit", "N", "S", "k0", "E", "C", "R0"])