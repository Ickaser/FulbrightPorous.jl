let
# ------------------
# Simtra data
# Kv fitting to FIg 7
NLSP = NonlinearLeastSquaresProblem

pKvdat = [ 50.0 4.18; 75.0 5.25;
100. 6.21; 125. 7.08; 150. 7.81;
175. 8.53; 200. 9.09; 250. 10.0; ] # WebPlotDigitized from Fig 7 of IntJPharm
pd = pKvdat[:,1]*u"mTorr"
Kvd = pKvdat[:,2]*u"J/hr/cm^2/K"
function fitKv(KKK) 
    Kc, Kp, Kd = KKK .* (u"W/m^2/K", u"W/m^2/K/Torr", u"1/Torr")
    Kv = RpFormFit(Kc, Kp, Kd)
    err = @. ustrip(u"W/m^2/K", Kv(pd) - Kvd)
    # return sum(abs2.(err)), Kv
    return err, Kv
end
optKv = solve(NLSP((x,p)->fitKv(x)[1], [10.0, 100, 1]), LevenbergMarquardt())
K_shf = fitKv(optKv.u)[2]

safesave(datadir("exp_pro", "simtra_kv.jld2"), @strdict K_shf)

end
# --------------------
# Gravimetric analysis for amber 10mL vials in MicroFD
Kv_all_1, Δm_1, int_dT_theory_1, int_dT_data_1 = let
using QuadGK

Ap, Av = π .* get_vial_radii("10R").^2 
hf0 = 2u"mL"/Ap
# Read the CSV file into a table
grav_raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-01-27-09_MFD_ISW_1.csv"), Table)
grav_dat = map(x->(i=x.Name, mi=x.var"Mass before (g)"*u"g",
    mf=x.var"Mass after (g)"*u"g", mf2=x.var"Mass Day After"*u"g",
    tc=!ismissing(x.Thermocouple)), 
    grav_raw)
filter!(x->x.i>0, grav_dat)
grav_notc = filter(x->!x.tc, grav_dat)

# @df grav_dat scatter(:i, :mf-:mi, c=:tc.+1)
# @df grav_dat scatter!(:i, :mf2-:mi, c=:tc.+1)

grav_Tdat = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-01-27-09_MFD_ISW_2.csv"), header=8, Table)
grav_Tdat = map(x->microfd_fulbright_columnrename(x, grav_Tdat.CycleTime[1]), grav_Tdat)


grav_pd = filter(x->(0u"mTorr" < x.pch_cm < 110u"mTorr"), grav_Tdat)
grav_pd.t .-= grav_pd.t[1]


Tsh = RampedVariable([228.15u"K", 253.15u"K"], 0.5u"K/minute")
pch = RampedVariable(100u"mTorr")
Tsub = LyoPronto.calc_Tsub(pch(0))
integrand = t-> max(0u"K", Tsh(t) - Tsub)
int_dT_theory = quadgk(integrand, 0.0u"hr", 4.0u"hr")[1]

Tf_i = @df grav_pd LinearInterpolation(:T1, :t)
Tsh_i = @df grav_pd LinearInterpolation(:Tsh_i, :t)
int_dT_data = quadgk(t->(Tsh_i(t) - Tf_i(t)), grav_pd.t[begin], grav_pd.t[end])[1]

# When using the Tf from thermocouple, we are already imitating very closely the theory for Kv.
# The temperature drop from bottom of ice to top is not zero, but if we explicitly account for it after the fact,
# the resulting Tsub is very very close to the frost point of the chamber pressure

Δm = [x.mi-x.mf2 for x in grav_notc]
mean_Δm = sum(Δm)/length(Δm)


# begin
# pl1 = @df grav_pd plot(:t, :Tsh_i, c=1)
# @df grav_pd plot!(:t, :T1, c=1)
# @df grav_pd plot!(:t, :T2, c=1)
# plot!(Tsh, tmax=4u"hr", c=2)
# plot!([0u"hr", 4u"hr"], [Tsub,Tsub], c=2)
# end

Kv = ConstPhysProp(LyoPronto.ΔH*mean_Δm/int_dT_data/Av |> u"W/m^2/K")

# This is slightly wrong--not accounting for decreasing hf--but want to see it anyway
Tsub_i = @. (grav_pd.T1 - grav_pd.Tsh_i)*Kv.val*Av/LyoPronto.k_ice/Ap*hf0 + grav_pd.T1

begin
plT = @df grav_pd plot(:t, :Tsh_i)
@df grav_pd plot!(:t, :T1)
@df grav_pd plot!(:t, :T2)
plT2 = @df grav_pd plot(:t, :T1, c=2)
@df grav_pd plot!(:t, Tsub_i, c=4)
hline!([Tsub], c=5)
plp = @df grav_pd plot(:t, :pch_cm, ylim=(0, 200))
@df grav_pd plot!(:t, :pch_pir)
plot(plT, plT2, plp, link=:x, layout=(3,1))
end

Kv_all = @. LyoPronto.ΔH*Δm/int_dT_data/Av|> u"W/m^2/K"


safesave(datadir("exp_pro", "gravimetric_kv.jld2"), @strdict Kv)

Kv_all, Δm, int_dT_theory, int_dT_data
end


target_Δm = 1u"g"
Tsh_hi = RampedVariable([233.15u"K", 283.15u"K"], 1u"K/minute")
integ_hi = t-> max(0u"K", Tsh_hi(t) - Tsub)
int_dT_time = t-> quadgk(integ_hi, 0.0u"hr", t)[1]
nlprob = IntervalNonlinearProblem{false}((0.0, 4.0)) do t, _
    dT = int_dT_time(t*u"hr")
    Δm = Kv(0)*dT*Av/LyoPronto.ΔH
    return ustrip(u"g", Δm - target_Δm)
end
t_target = solve(nlprob, ITP()).u[1]


Kv_all_2, Δm_2, int_dT_theory_2, int_dT_data_2 = let
using QuadGK

Av = π * get_vial_radii("10R")[2]^2 
# Read the CSV file into a table
grav_raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-06-11-10_MFD_ISW_1.csv"), Table)
grav_dat = map(x->(i=x.var"Vial num", mi=x.var"Mass before"*u"g",
    mf=x.var"Mass immediately after"*u"g", mf2=x.var"Mass later"*u"g",
    tc=x.TC), 
    grav_raw)
grav_notc = filter(x->!x.tc, grav_dat)

# @df grav_dat scatter(:i, :mf-:mi, c=:tc.+1)
# @df grav_dat scatter!(:i, :mf2-:mi, c=:tc.+3)
# @df grav_dat histogram(:mf-:mi, bins=10, alpha=0.5)
# @df grav_dat histogram!(:mf2-:mi, bins=10, alpha=0.5)

grav_Tdat = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-06-11-10_MFD_ISW_2.csv"), header=8, Table)
grav_Tdat = map(x->microfd_fulbright_columnrename(x, grav_Tdat.CycleTime[1]), grav_Tdat)


grav_pd = filter(x->(0u"mTorr" < x.pch_cm < 110u"mTorr"), grav_Tdat)
grav_pd.t .-= grav_pd.t[1]

begin
plT = @df grav_pd plot(:t, :Tsh_i)
@df grav_pd plot!(:t, :T1)
@df grav_pd plot!(:t, :T2)
plp = @df grav_pd plot(:t, :pch_cm, ylim=(0, 200))
@df grav_pd plot!(:t, :pch_pir)
plot(plT, plp, link=:x, layout=(2,1))
end

Tf_i = @df grav_pd LinearInterpolation(:T2, :t)
Tsh_i = @df grav_pd LinearInterpolation(:Tsh_i, :t)
int_dT_data = quadgk(t->(Tsh_i(t) - Tf_i(t)), grav_pd.t[begin], grav_pd.t[end])[1]

Δm = [x.mi-x.mf for x in grav_notc]
mean_Δm = sum(Δm)/length(Δm)

# Replace this with an integration across the data once I have it...
Tsh = RampedVariable([233.15u"K", 283.15u"K"], 1u"K/minute")
pch = RampedVariable(100u"mTorr")
Tsub = LyoPronto.calc_Tsub(pch(0))
integrand = t-> max(0u"K", Tsh(t) - Tsub)
int_dT_theory = quadgk(integrand, 0.0u"hr", grav_pd.t[end])[1]

# begin
# pl1 = @df grav_pd plot(:t, :Tsh_i, c=1)
# @df grav_pd plot!(:t, :T1, c=1)
# @df grav_pd plot!(:t, :T2, c=1)
# plot!(Tsh, tmax=grav_pd.t[end], c=2)
# plot!([0u"hr", grav_pd.t[end]], [Tsub,Tsub], c=2)
# end

Kv = ConstPhysProp(LyoPronto.ΔH*mean_Δm/int_dT_data/Av |> u"W/m^2/K")
Kv_all = LyoPronto.ΔH*Δm/int_dT_data/Av .|> u"W/m^2/K"

safesave(datadir("exp_pro", "gravimetric_kvhi.jld2"), @strdict Kv)

Kv_all, Δm, int_dT_theory, int_dT_data
end

begin
Kv_fit = 18.997379673262678u"W*K^-1*m^-2"
boxplot([Kv_all_1, Kv_all_2], linewidth=1, label=["Gravimetric 1" "Gravimetric 2"])
scatter!(fill(1, size(Kv_all_1)), Kv_all_1, c=:black, ma=0.5, label="")
scatter!(fill(2, size(Kv_all_2)), Kv_all_2, c=:black, ma=0.5, label="")
hline!([Kv_fit], label="fit to data", lw=3, c=3)
plot!(size=(200,400), xwiden=1.5, xticks=[], ylabel="K_v", unitformat=latexify)
plot!(legend=:outerbottom)
end
savefig(plotsdir("gravimetric_vs_fit_fulbright.svg"))
