fp_plot_defaults()
Latexify.set_default(labelformat=:square)

ful_suc_fits = load(datadir("exp_pro", "fulbright_sucrose_fitresults.jld2"), "fits")

k_ful_sigm_5 = ful_suc_fits[1].kf
k_ful_sigm_10 = ful_suc_fits[2].kf
k_ful_sigm = ful_suc_fits[3].kf

sim_suc_cases = load(datadir("exp_pro", "simtra_sucrose_fitresults.jld2"), "mcfit_cases")
sim_suc_fit = load(datadir("exp_pro", "simtra_sucrose_fitresults.jld2"), "mcfits")

k_sim_sigm = transform(sim_suc_cases[11].trans, sim_suc_fit[11]).k_f
k_sim_intp = transform(sim_suc_cases[3].trans, sim_suc_fit[3]).k_f
k_sim_myeg = transform(sim_suc_cases[7].trans, sim_suc_fit[7]).k_f

begin
pl1 = plot(u"°C", u"s^-1", xlabel="T\\;", ylabel="k \\, ", unitformat=latexify)
kineticplot!(k_sim_sigm, label="BSA+Sucrose, sigmoid", c=4)
kineticplot!(k_sim_intp, label="BSA+Sucrose, interpolation", c=4, marker=:square)
# kineticplot!(k_sim_myeg, label="BSA+Sucrose, MYEGA", c=4, marker=:diamond)
kineticplot!(k_ful_sigm, label="Sucrose (5% & 10%), sigmoid", c=5, )
kineticplot!(k_ful_sigm_5, label="Sucrose (5%), sigmoid", c=5, linestyle=:dash)
kineticplot!(k_ful_sigm_10, label="Sucrose (10%), sigmoid", c=5, linestyle=:dot)
plot!(pl1, legend=(0.66, 0.5), size=(650,300), right_margin=100Plots.px)
plot!(pl1, bottom_margin=20Plots.px, left_margin=12Plots.px)
# pl2 = deepcopy(pl1)
# plot!(pl2, xlim=(-38, -32), ylim=(1e-13, 1e-10), legend=:none, ylabel="")
# plot(pl1, pl2, layout=@layout([a b{0.3w}]))
end
savefig(plotsdir("rate_constants.svg"))
savefig(plotsdir("rate_constants.pdf"))
# Fix the legend before exporting to PDF

k_arrhenius = T->1e-7u"s^-1" * exp(-2000u"K"/T) # Near-maximum slope for Arrhenius
k_arrhenius(-25u"°C"|> u"K")/ k_arrhenius(-45u"°C" |> u"K")  
k_sim_intp(-25u"°C"|> u"K")/ k_sim_intp(-45u"°C" |> u"K")