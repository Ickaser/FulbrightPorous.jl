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
let
using QuadGK

Av = π * get_vial_radii("10R")[2]^2 
# Read the CSV file into a table
grav_raw = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-01-27-09_MFD_ISW_1.csv"), Table)
grav_dat = map(x->(i=x.Name, mi=x.var"Mass before (g)"*u"g",
    mf=x.var"Mass after (g)"*u"g", mf2=x.var"Mass Day After"*u"g",
    tc=!ismissing(x.Thermocouple)), 
    grav_raw)
filter!(x->x.i>0, grav_dat)
grav_notc = filter(x->!x.tc, grav_dat)

@df grav_dat scatter(:i, :mf-:mi, c=:tc.+1)
@df grav_dat scatter!(:i, :mf2-:mi, c=:tc.+1)

grav_Tdat = CSV.read(datadir("exp_raw", "Fulbright2025", "2025-01-27-09_MFD_ISW_2.csv"), header=8, Table)
grav_Tdat = map(x->microfd_fulbright_columnrename(x, grav_Tdat.CycleTime[1]), grav_Tdat)


istart = findfirst(0u"mTorr" .< grav_Tdat.pch_cm .< 110u"mTorr")
iend = findlast(0u"mTorr" .< grav_Tdat.pch_cm .< 110u"mTorr")
grav_pd = grav_Tdat[istart:iend]

begin
plT = @df grav_pd plot(:t, :Tsh_i)
@df grav_pd plot!(:t, :T1)
@df grav_pd plot!(:t, :T2)
plp = @df grav_pd plot(:t, :pch_cm, ylim=(0, 200))
@df grav_pd plot!(:t, :pch_pir)
plot(plT, plp, link=:x, layout=(2,1))
end

Tf_i = @df grav_pd LinearInterpolation(:T1, :t)
Tsh_i = @df grav_pd LinearInterpolation(:Tsh_i, :t)
int_dT_data = quadgk(t->(Tsh_i(t) - Tf_i(t)), grav_pd.t[begin], grav_pd.t[end])[1]

Δm = [x.mi-x.mf2 for x in grav_notc]
mean_Δm = sum(Δm)/length(Δm)

# Replace this with an integration across the data once I have it...
Tsh = RampedVariable([228.15u"K", 253.15u"K"], 1u"K/minute")
pch = RampedVariable(100u"mTorr")
nlprob = NonlinearProblem((T,_)->ustrip(u"mTorr", calc_psub(T*u"K")-pch(0)), 240.0)
Tsub = solve(nlprob).u[1]*u"K"
integrand = t-> max(0u"K", Tsh(t) - Tsub)
int_dT_theory = quadgk(integrand, 0.0u"hr", 4.0u"hr")[1]

Kv = ConstPhysProp(LyoPronto.ΔH*mean_Δm/int_dT_theory/Av |> u"W/m^2/K")

safesave(datadir("exp_pro", "gravimetric_kv.jld2"), @strdict Kv)

end
