using DrWatson
@quickactivate :FulbrightPores

const LP = LyoPronto

pythonplot()
default(:fontfamily, "Arial")
default(:framestyle, :box)
default(:lw, 3)

# ----------------------------------------
# Set up ODES

vialsize = "6R"
rad_i, rad_o = get_vial_radii(vialsize)
A_p = π*rad_i^2  # cross-sectional area inside the vial
A_v = π*rad_o^2 # vial bottom area


c_solid = 0.05u"g/mL" # g solute / mL solution
ρ_solution = 1u"g/mL" # g/mL total solution density
R0 = 1.4u"cm^2*hr*Torr/g"
A1 = 16.0u"cm*hr*Torr/g"
A2 = 0.0u"1/cm"
Rp = RpFormFit(R0, A1, A2)
Vfill = 5u"mL"
pch = RampedVariable(100u"mTorr")

KC = 2.75e-4u"cal/s/K/cm^2"
KP = 8.93e-4u"cal/s/K/cm^2/Torr"
KD = 0.46u"1/Torr"
K_shf_f = RpFormFit(KC, KP, KD)

T_f0 = 235.15u"K"

h_f0 = Vfill/A_p
ϕ = (ρ_solution - c_solid)/ρ_solution



# controls = Dict{Symbol, RampedVariable}()
pch_rv = RampedVariable(100u"mTorr")
Tsh_rv = RampedVariable([233.15u"K", 283.15u"K"], [0.5u"K/minute"], [])

params_bunch_pikal = Vector{Any}([
    (Rp, h_f0, c_solid, ρ_solution),
    (K_shf_f, A_v, A_p),
    (pch_rv, Tsh_rv), # P_per_vial in W
])
# -----------------

u0_pikal = ustrip.([u"cm", u"K"], [h_f0, T_f0])
# u0 = [ustrip(u"g*K", mpTp0), ustrip(u"K", Tw0), ustrip(u"cm", L0)]

tspan = (0.0, ustrip(u"hr", 80u"hr"))
p_pikal = tuple(params_bunch_pikal...)

pprob = ODEProblem(lyo_1d_dae_f, u0_pikal, tspan, p_pikal; callback=end_drying_callback)
@time psol = solve(pprob, Rodas3(autodiff=false))

# --------------

k0 = 0u"1/s"
E = 0u"J/mol"
C = 2e2u"m"
# typical diameter: 50 μm
d_typical = 50u"μm"
V_typical = 4/3*π*d_typical^3
N_typical = Vfill/V_typical
# typical total volume: 
N = RampedVariable(N_typical / h_f0)
S = RampedVariable(1e-1u"cm^5")
# R0 = 0.1u"cm^2*Torr*hr/g"

μ2_0 = (100u"μm")^6
u0_mc = ustrip.([u"cm", u"K", NoUnits, u"cm^6"], [h_f0, T_f0, 1.0e2, μ2_0])

params_bunch_mc = Vector{Any}([
    (Rp, h_f0, c_solid, ρ_solution),
    (K_shf_f, A_v, A_p),
    (pch_rv, Tsh_rv), 
    (k0, E, C, N, S, R0)
])

p_mc = tuple(params_bunch_mc...)
po_mc = ParamObjMC(params_bunch_mc)

prob = ODEProblem(lyo_mc_dae_f, u0_mc, tspan, p_mc; callback=end_drying_callback)
@time sol = solve(prob, Rodas3(autodiff=false))

plot(u"hr", u"°C")
modconvtplot!(sol, psol, labels=[L"$T_f$, PBE." L"$T_f$, trad."])
plot!(ylabel="Temp [°C]")
savefig(plotsdir("quicktest_T.svg"))

function calc_Rp(hf, Tf, μ0, μ2)
    hd = (h_f0-hf*u"cm")
    Tf *= u"K"
    μ2 *= u"cm^6"
    ϕ = (ρ_solution - c_solid)/ρ_solution
    μ1 = hd*ϕ*A_p
    # μ32 = μ0^(-1//8)*μ1^(3//4)*μ2^(3//8)
    # Rp = R0 + C * sqrt(R*Tf/Mw) * μ1*hd^(3/2)/μ32/ϕ
    Rp = R0 + C * η*R*Tf/Mw/((calc_psub(Tf) + pch(0))/2) * μ1*hd/μ2/ϕ
    return hd, Rp
end


ls = range(0u"m", h_f0, length=50)
plot(sol, idxs=(calc_Rp, 1, 2, 3, 4), label="Rp from PBE")
plot!(ls, Rp.(ls), label="Rp as usual", title="Mannitol R_p")
savefig(plotsdir("quicktest_Rp.svg"))

fitdat = PrimaryDryFit(psol.t*u"hr", (psol[2,:]*u"K",))

plot(fitdat)
