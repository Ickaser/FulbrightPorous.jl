@quickactivate :FulbrightPores
pythonplot()
const NLSP = NonlinearLeastSquaresProblem
const LM = LevenbergMarquardt
const LMa = LevenbergMarquardt(autodiff=AutoFiniteDiff())
const LMf = LevenbergMarquardt(autodiff=AutoForwardDiff())

# Information from Fissore and paper
Kv = RpFormFit(7.78u"W/m^2/K", 1.44u"W/m^2/K/Pa", 0.02u"Pa^-1")

const psub_Tt = 273.16u"K" # Triple point T, kelvin
const psub_pt = 611.654_771u"Pa" # Triple point p, Pa
const psub_ei = [
 20.996_966_510_789_7,
  3.724_374_782_713_62,
-13.920_548_321_552_4,
 29.698_876_501_356_6,
-40.197_239_263_594_4,
 29.788_048_105_021_5,
 -9.130_509_635_477_21,
]
"""
    calc_psub_highorder(T)

Compute pressure (in Pascals) of sublimation at temperature `T` in Kelvin

From Feistel and Wagner, 2006
"""
function calc_psub_highorder(T)
    if clamp(T, 20.0u"K", 273.0u"K") != T
        @warn "Invalid temperature for sublimation pressure correlation. Clamped to 20 and 273 K" T
    end
    η = sum(psub_ei .* (T/psub_Tt) .^ (0:6))
    lnπ = 1.5log(T/psub_Tt) + (1 - psub_Tt/T) * η
    return exp(lnπ) * psub_pt
end

function calc_pwi(T)
    
    t1 = -9.09718*(psub_Tt/T - 1)
    t2 = -3.56654*log10(psub_Tt/T)
    t3 =  0.876793*(1 - T/psub_Tt)
    t4 = log10(6.1071)
    return 100u"Pa" * exp(t1 + t2 + t3 + t4)
end
Tran = range(223.15u"K", 273.15u"K", length=51)
begin
plot(u"°C", u"Pa", xlabel="T_{sub}", ylabel="p_{sub}", unitformat=latexsquareunitlabel, legend=:topleft)
plot!(Tran, calc_pwi, label="data nella tesi")
plot!(Tran, calc_psub, label="used in LyoPronto")
plot!(Tran, calc_psub_highorder, label="high order from literature")
plot!(size=(400,300), legend=:topleft)
end
savefig(plotsdir("checking_psub.png"))

pch1 = RampedVariable(5u"Pa")
pch2 = RampedVariable(5u"Pa")
pch3 = RampedVariable(20u"Pa")
Ap = Av = π*(7.125u"mm")^2
hf0 = 1.5u"mL"/Ap |> u"cm"
c_solid = 0.05u"g/mL"
ρ_solution = 1u"g/mL"
# Rpg = RpFormFit(0.5u"cm^2/Torr/hr/g", 14u"cm/Torr/hr/g", 0.1u"1/cm")
Rp_t1 = RpFormFit(2.11e4u"m/s", 1.62e8u"s^-1", 1.55e3u"m^-1")
Rp_t2 = RpFormFit(1.10e4u"m/s", 6.19e7u"s^-1", 5.00e3u"m^-1")
Rp_t3 = RpFormFit(5.05e3u"m/s", 2.77e7u"s^-1", 1.49e3u"m^-1")

# Unit add functions
fissore_units1(row) = (t=row.t*u"hr", Tf=row.Tf*u"°C")
fissore_units2(row) = (t=row.t*u"hr", Tsh=row.Tsh*u"°C")

# -------------------

# Read data for Cycle 1
dat1a = map(fissore_units1, CSV.read(datadir("exp_raw", "Fissore2015", "1_Tf.csv"), Table))
dat1b = map(fissore_units2, CSV.read(datadir("exp_raw", "Fissore2015", "1_Tsh.csv"), Table))
# @df dat1a plot(:t, :Tf)

fitdat1 = @df filter(x->(0u"s"<x.t<23u"hr"), dat1a) PrimaryDryFit(:t, (:Tf,), 25u"hr")
plot(fitdat1)

Tsh1 = RampedVariable([-40, -20, -25, -30]u"°C".|>u"K", [1, 1, 1]u"K/minute", [1, 3.5]u"hr")
# @df dat1b plot(:t, :Tsh, xlim=(-2, 10), ylim=(-50, -18))
# plot!(Tsh1)
 
# WebPlotDigitized for Cycle 1 pressure rise test
lRp_mat1 = [
0.05294117647058827 0.13064516129032233
0.6004524886877828 0.6048387096774193
0.7452488687782804 1.3499999999999996
0.8696832579185519 2.061290322580645
0.9285067873303167 2.75
0.9941176470588236 3.337096774193548
1.0054298642533936 3.879032258064517
1.0280542986425338 4.466129032258065
1.0303167420814479 5.019354838709677
1.0303167420814479 5.583870967741936
1.04841628959276 6.103225806451613
1.050678733031674 6.622580645161291
1.0665158371040724 7.108064516129032
1.0755656108597285 7.6161290322580655
1.1027149321266967 8.011290322580646
]
# Digitized with x and y switched, accidentally, so fix that in the table
lRp_1 = Table(map(x->(l=x[2]*u"cm", Rp=x[1]*1e5u"m/s"), eachrow(lRp_mat1)))
# @df lRp_1 scatter(:l, :Rp)
po_1p = ParamObjPikal([
    (Rp_t1, hf0, c_solid, ρ_solution),
    (Kv, Av, Ap),
    (pch1, Tsh1)])


# --------------------
# Read data for Cycle 2
dat2a = map(fissore_units1, CSV.read(datadir("exp_raw", "Fissore2015", "2_Tf.csv"), Table))
dat2b = map(fissore_units2, CSV.read(datadir("exp_raw", "Fissore2015", "2_Tsh.csv"), Table))
# @df dat2a plot(:t, :Tf)
fitdat2 = @df filter(x->(0u"s"<x.t<15u"hr"), dat2a) PrimaryDryFit(:t, (:Tf,), 17u"hr")
plot(fitdat2)

Tsh2 = RampedVariable([-40, -20]u"°C".|>u"K", 1u"K/minute")
@df dat2b plot(:t, :Tsh, xlim=(-2, 10), ylim=(-50, -18))
plot!(Tsh2)

po_2p = ParamObjPikal([
    (Rp_t2, hf0, c_solid, ρ_solution),
    (Kv, Av, Ap),
    (pch2, Tsh2)])

lRp_mat2 = [
0.4455958549222799 0.4276816608996541
0.694300518134715 0.5003460207612458
0.9430051813471501 0.6560553633217995
1.243523316062176 0.7079584775086506
1.5025906735751295 0.7910034602076126
1.782383419689119 0.853287197231834
2.072538860103627 0.7847750865051903
2.3212435233160624 0.7972318339100347
2.5803108808290154 0.8512110726643599
3.1139896373056986 0.8740484429065745
3.3937823834196887 0.8771626297577855
3.937823834196892 0.8989619377162631
4.487046632124354 0.9238754325259517
5.000000000000001 0.9176470588235295
5.782383419689119 0.8792387543252596
6.0310880829015545 0.9134948096885814
6.772020725388601 0.8750865051903115
7.222797927461141 0.9000000000000001
7.647668393782382 0.8740484429065745
8.036269430051815 0.8948096885813149
]
lRp_2 = Table(map(x->(l=x[1]*u"cm", Rp=x[2]*1e5u"m/s"), eachrow(lRp_mat2)))

# --------------
# Read data for Cycle 3
dat3a = map(fissore_units1, CSV.read(datadir("exp_raw", "Fissore2015", "3_Tf.csv"), Table))
dat3b = map(fissore_units2, CSV.read(datadir("exp_raw", "Fissore2015", "3_Tsh.csv"), Table))
@df dat3a plot(:t, :Tf)

fitdat3 = @df filter(x->(0u"s"<x.t<9u"hr"), dat3a) PrimaryDryFit(:t, (:Tf,), 12u"hr")
plot(fitdat3)

# Match a ramped variable to the reported Tsh
Tsh3 = RampedVariable([-43, -20, -30, -25]u"°C".|>u"K", [1, -0.5, 1]u"K/minute", [0, 5.7]u"hr")
# @df dat3b plot(:t, :Tsh, xlim=(-2, 10), ylim=(-50, -15))
# plot!(Tsh3)

po_3p = ParamObjPikal([
    (Rp_t3, hf0, c_solid, ρ_solution),
    (Kv, Av, Ap),
    (pch3, Tsh3)])



# WebPlotDigitized for Cycle 3 pressure rise test
lRp_mat3 = [
0.5695623775310256 0.11204188481675437
0.992815153494448 0.1445026178010469
1.588504245591117 0.15916230366492146
2.544741998693664 0.206282722513089
3.1979098628347487 0.2282722513089004
3.8667537557152185 0.2261780104712039
4.519921619856302 0.23350785340314117
5.141737426518616 0.19999999999999996
5.612018288700196 0.23769633507853405
6.082299150881776 0.19476439790575895
6.599608099281515 0.13193717277486883
6.991508817766166 0.2178010471204186
7.378184193337687 0.17277486910994755
7.770084911822338 0.12460732984293199
8.047028086218159 0.14345549738219932
]
lRp_3 = Table(map(x->(l=x[1]*u"cm", Rp=x[2]*1e5u"m/s"), eachrow(lRp_mat3)))

# ------------
# Compare all three Rps, to check digitization

begin
plot(u"cm", u"m/s", xlabel="l", ylabel="R_p", unitformat=latexsquareunitlabel, legend=:topleft)
# blankplot_lRp()
# plot!(xlim=(0,8.1), ylim=(0, 1.2e5))
@df lRp_1 scatter!(:l, :Rp, msize=8, marker=:square, label="Cycle 1")
@df lRp_2 scatter!(:l, :Rp, msize=8, marker=:circle, label="Cycle 2")
@df lRp_3 scatter!(:l, :Rp, msize=8, marker=:utriangle, label="Cycle 3")
plot!(l_sm, Rp_t1.(l_sm), c=1, label="3PP, Tesi")
plot!(l_sm, Rp_t2.(l_sm), c=2, label="3PP, Tesi")
plot!(l_sm, Rp_t3.(l_sm), c=3, label="3PP, Tesi")
end
savefig(plotsdir("checking_Rp_fissore.png"))

# ----------------------
# Use values of Rp reported in thesis to check simulation
sol_1t = solve(ODEProblem(po_1p), Rodas3())
sol_2t = solve(ODEProblem(po_2p), Rodas3())
sol_3t = solve(ODEProblem(po_3p), Rodas3())

begin
# first
pl1 = blankplot_hrC(xlim=(-0, 50), ylim=(-50, -20))
@df dat1a plot!(:t, :Tf, label="Tf, exp")
modconvtplot!(sol_1t)
@df dat1b plot!(:t, :Tsh, label="Tsh data")
plot!(Tsh1, label="Tsh interpolator")
# second
pl2 = blankplot_hrC(xlim=(-0, 20), ylim=(-50, -15))
@df dat2a plot!(:t, :Tf, label="Tf, exp")
modconvtplot!(sol_2t)
@df dat2b plot!(:t, :Tsh, label="Tsh data")
plot!(Tsh2, label="Tsh interpolator")
# third
pl3 = blankplot_hrC(xlim=(-0, 20), ylim=(-45, -15))
@df dat3a plot!(:t, :Tf, label="Tf, exp")
modconvtplot!(sol_3t)
@df dat3b plot!(:t, :Tsh, label="Tsh data")
plot!(Tsh3, label="Tsh interpolator")
# combine
plot(pl1, pl2, pl3, layout=(3,1), size=(600, 600))
end
savefig(plotsdir("checking_T_fissore.png"))

# -----------------------

# Traditional approach: fit Rp to the pressure rise data
Rp_un = [u"cm^2*Torr*hr/g", u"cm*Torr*hr/g", u"cm^-1"]
function Rpfit(A012, dat)
    Rp = RpFormFit((A012.*Rp_un)...)
    error = ustrip.(u"cm^2*Torr*hr/g"^2, sum(abs2, Rp.(dat.l) .- dat.Rp))
    # @show Rp
    return [error]
end
p0 = [0.1, 0.5, 0.1]
Rpg = RpFormFit((p0.*Rp_un)...)
optRp1 = solve(NLSP(Rpfit, p0, p=lRp_1), LMf, show_trace=Val(true))
optRp2 = solve(NLSP(Rpfit, p0, p=lRp_2), LMf, show_trace=Val(true))
optRp3 = solve(NLSP(Rpfit, p0, p=lRp_3), LMf)
Rp1 = RpFormFit((optRp1.u.*Rp_un)...)
Rp2 = RpFormFit((optRp2.u.*Rp_un)...)
Rp3 = RpFormFit((optRp3.u.*Rp_un)...)

# -----------------------
# Try again with the interpolations I fit:
po_1Rp = @set po_1p.Rp = Rp1
po_2Rp = @set po_2p.Rp = Rp2
po_3Rp = @set po_3p.Rp = Rp3
sol_1Rp = solve(ODEProblem(po_1Rp), Rodas3())
sol_2Rp = solve(ODEProblem(po_2Rp), Rodas3())
sol_3Rp = solve(ODEProblem(po_3Rp), Rodas3())

begin
# first
pl1 = blankplot_hrC(xlim=(-0, 50), ylim=(-50, -20))
@df dat1a plot!(:t, :Tf, label="Tf, exp")
modconvtplot!(sol_1Rp)
@df dat1b plot!(:t, :Tsh, label="Tsh data")
plot!(Tsh1, label="Tsh interpolator")
# second
pl2 = blankplot_hrC(xlim=(-0, 20), ylim=(-50, -15))
@df dat2a plot!(:t, :Tf, label="Tf, exp")
modconvtplot!(sol_2Rp)
@df dat2b plot!(:t, :Tsh, label="Tsh data")
plot!(Tsh2, label="Tsh interpolator")
# third
pl3 = blankplot_hrC(xlim=(-0, 20), ylim=(-45, -15))
@df dat3a plot!(:t, :Tf, label="Tf, exp")
modconvtplot!(sol_3Rp)
@df dat3b plot!(:t, :Tsh, label="Tsh data")
plot!(Tsh3, label="Tsh interpolator")
# combine
plot(pl1, pl2, pl3, layout=(3,1), size=(600, 600))
end
savefig(plotsdir("checking_T2_fissore.png"))


# ------------------------------

# Fitting to temperatures
function gen_sol_Rp(a012, po)
    new_params = @set po.Rp = RpFormFit((a012.*Rp_un)...)
    prob = ODEProblem(new_params)
    sol = solve(prob, Rodas3())
    return sol
end
# obj_Rp(x, p) = [any(x .< 0) ? obj_expT(gen_sol_Rp(x, p[1]), p[2]) : NaN]
obj_Rp(x, p) = [obj_expT(gen_sol_Rp(x, p[1]), p[2])]
function gen_sol_KRp(Ka012, po)
    new_params = @set po.Rp = RpFormFit((Ka012[2:4].*Rp_un)...)
    @reset new_params.Kshf = ConstPhysProp(Ka012[1]*u"W/m^2/K")
    prob = ODEProblem(new_params)
    sol = solve(prob, Rodas3())
    return sol
end
# obj_KRp(x, p) = [any(x[2:4] .< 0) ? obj_expT(gen_sol_KRp(x, p[1]), p[2]) : NaN]
obj_KRp(x, p) = [obj_expT(gen_sol_KRp(x, p[1]), p[2])]



l_sm = range(0u"cm", 8u"cm", length=100)

begin
blankplot_lRp()
plot!(xlim=(0,8.1), ylim=(-0.1, 3), legend=:topleft)
@df lRp_1 scatter!(:l, :Rp, c=1, marker=:square, label="Cycle 1, PR")
plot!(l_sm, Rp1.(l_sm), c=1, label="Cycle 1 3PP fit to PR")
plot!(l_sm, Rp_t1.(l_sm), c=1, ls=:dash, label="thesis fit")
@df lRp_2 scatter!(:l, :Rp, c=2, marker=:circle, label="Cycle 2, PR")
plot!(l_sm, Rp2.(l_sm), c=2, label="Cycle 2 3PP fit to PR")
plot!(l_sm, Rp_t2.(l_sm), c=2, ls=:dash, label="thesis fit")
@df lRp_3 scatter!(:l, :Rp, c=3, marker=:utriangle, label="Cycle 3, PR")
plot!(l_sm, Rp3.(l_sm), c=3, label="Cycle 3 3PP fit to PR")
plot!(l_sm, Rp_t3.(l_sm), c=3, ls=:dash, label="thesis fit")
end
savefig(plotsdir("FissorePisano2015_3PP_Rp_PRfits_tesi.png"))


Rpg_nondim = ustrip.(Rp_un, [Rp1.R0, Rp1.A1, Rp1.A2])
KRpg_nondim = [14, 2, 14, 0.1]
solt = gen_sol_KRp(KRpg_nondim, po_1p)
modconvtplot(solt)
plot!(fitdat1)
opt_KRp1 = solve(NLSP(obj_KRp, KRpg_nondim, p=(po_1p, fitdat1)), LM(), show_trace=Val(true))
opt_KRp2 = solve(NLSP(obj_KRp, KRpg_nondim, p=(po_2p, fitdat2)), LM(), show_trace=Val(true))
opt_KRp3 = solve(NLSP(obj_KRp, KRpg_nondim, p=(po_3p, fitdat3)), LM(), show_trace=Val(true))
sol_1pK = gen_sol_KRp(opt_KRp1.u, po_1p)
sol_2pK = gen_sol_KRp(opt_KRp2.u, po_2p)
sol_3pK = gen_sol_KRp(opt_KRp3.u, po_3p)
opt_Rp1 = solve(NLSP(obj_Rp, Rpg_nondim, p=(po_1p, fitdat1)), LM(), show_trace=Val(true))
opt_Rp2 = solve(NLSP(obj_Rp, Rpg_nondim, p=(po_2p, fitdat2)), LM(), show_trace=Val(true))
opt_Rp3 = solve(NLSP(obj_Rp, Rpg_nondim, p=(po_3p, fitdat3)), LM(), show_trace=Val(true))
sol_1p = gen_sol_Rp(opt_Rp1.u, po_1p)
sol_2p = gen_sol_Rp(opt_Rp2.u, po_2p)
sol_3p = gen_sol_Rp(opt_Rp3.u, po_3p)
sol_1p.prob.p.Kshf
sol_2p.prob.p.Kshf
sol_3p.prob.p.Kshf

sols = [sol_1p, sol_2p, sol_3p, sol_1pK, sol_2pK, sol_3pK]
nms = ["1 Rp", "2 Rp", "3 Rp", "1 KvRp", "2 KvRp", "3 KvRp"]
prmstab = Table(map(sols, nms) do sol, nm
    Kv = sol.prob.p.Kshf(sol.prob.p.pch(0))
    Rp = sol.prob.p.Rp
    Rpmax = Rp(hf0)
    (;R0, A1, A2) = Rp
    return (case=nm, Kv=Kv, R0=R0, A1=A1, A2=A2, Rpmax=Rpmax)
end )

begin
labs = permutedims(["3PP Rp", "3PP Kv-Rp"])
pl1 = blankplot_hrC(legend=:bottomright)
modconvtplot!(sol_1p, sol_1pK, labels=labs)
plot!(fitdat1)
pl2 = blankplot_hrC(legend=:bottomright)
modconvtplot!(sol_2p, sol_2pK, labels=labs)
plot!(fitdat2)
pl3 = blankplot_hrC(legend=:bottomright)
modconvtplot!(sol_3p, sol_3pK, labels=labs)
plot!(fitdat3)
plot(pl1, pl2, pl3, layout=(3,1), size=(600, 600))
end
savefig(plotsdir("FissorePisano2015_3PP_KRp_Tfits.png"))

begin
blankplot_lRp()
plot!(xlim=(0,8.1), ylim=(-0.1, 10), legend=:topleft)
@df lRp_1 scatter!(:l, :Rp, c=1, marker=:square, label="Cycle 1, PR")
plot!(l_sm, Rp1.(l_sm), c=1, label="Cycle 1 3PP fit to PR")
plot!(l_sm, sol_1pK.prob.p.Rp.(l_sm), c=1, ls=:dash, label="Cycle 1 3PP fit to T")
@df lRp_2 scatter!(:l, :Rp, c=2, marker=:circle, label="Cycle 2, PR")
plot!(l_sm, Rp2.(l_sm), c=2, label="Cycle 2 3PP fit to PR")
plot!(l_sm, sol_2pK.prob.p.Rp.(l_sm), c=2, ls=:dash, label="Cycle 2 3PP fit to T")
@df lRp_3 scatter!(:l, :Rp, c=3, marker=:utriangle, label="Cycle 3, PR")
plot!(l_sm, Rp3.(l_sm), c=3, label="Cycle 3 3PP fit to PR")
plot!(l_sm, sol_3pK.prob.p.Rp.(l_sm), c=3, ls=:dash, label="Cycle 3 3PP fit to T")
end
savefig(plotsdir("FissorePisano2015_3PP_Rp_PRfits.png"))