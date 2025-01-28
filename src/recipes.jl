export blankplot_hrC
blankplot_hrC(;kwargs...) = plot(u"hr", u"°C", xlabel="Time", ylabel="Temperature", unitformat=:square;kwargs...)

# ------------------------
# Visualize the pore size evolution
export blankplot_lRp
blankplot_lRp(;kwargs...) = plot(u"cm", u"cm^2*Torr*hr/g", xlabel=L"h_d", ylabel=L"R_p", unitformat=(x,y)->"$x [$(latexify(y))]";kwargs...)

@userplot PoreEvolPlot
@recipe function f(pep::PoreEvolPlot)
    sol = pep.args[1]
    t = sol.t*u"hr"
    hf = sol[1,:]*u"cm"
    μ0 = sol[3,:]
    μ2 = sol[4,:]*u"mm^6" .|> u"mm^6"
    po = sol.prob.p
    hd = (po.hf0 .- hf) .|> u"cm"
    ϕ = (po.ρ_solution - po.c_solid)/po.ρ_solution
    μ1 = hd*ϕ*po.Ap .|> u"mm^3"
    mean = μ1 ./ μ0
    var = μ2 ./ μ0 .- mean.^2
    Rp = map(sol.t, sol.u) do t,u
        lyo_mc_dae!(fill(0.0, 4), u, sol.prob.p, t) |> u"cm^2*Torr*hr/g"
    end
    seriestype --> :samplemarkers
    layout --> @layout [[a ; b; c;] d]
    @series begin
        subplot --> 1
        ylabel --> "number\nof pores" #*", number of pores"
        xlabel --> ""
        legend --> :none
        xformatter -->  _->""
        return hd, μ0
    end
    @series begin
        subplot --> 2
        ylabel --> "typical \npore size\n" #*", total pore volume"
        xlabel --> ""
        legend --> :none
        xformatter -->  _->""
        return hd, mean
    end
    @series begin
        subplot --> 3
        ylabel --> "pore size \nvariance\n" #*", pore size second moment"
        xlabel --> "Dry layer thickness"
        legend --> :none
        return hd, var
    end
    @series begin
        subplot --> 4
        ylabel --> "Rp"
        xlabel --> "Dry layer thickness"
        legend --> :bottomright
        # return MCRpPlot((sol,))
        return hd, Rp
    end
end

@userplot MCRpPlot
@recipe function f(mrp::MCRpPlot)
    sol = mrp.args[1]
    t = sol.t*u"hr"

    Rp = map(sol.t, sol.u) do ti,ui
        lyo_mc_dae!(fill(0.0, 4), ui, sol.prob.p, ti) |> u"m/s"
    end
    @series begin
        ylabel --> "Rp"
        label --> "\$R_p\$"
        t, Rp
    end
end

function calc_l(dat, hf0)
    minteg = @df dat cumsum(:mfr[1:end-1].*(:t[2:end]-:t[1:end-1]))
    insert!(minteg, 1, 0u"g")
    frac = minteg ./ minteg[end]
    ld = frac .*hf0 
end
function calc_lRp(sol)
    l = (sol[1,1] .- sol[1,:])*u"cm"
    Rp_mc = [FulbrightPores.lyo_mc_dae!(fill(0.0, 4), u, sol.prob.p, t) for (u,t) in zip(sol.u, sol.t)]
    return l, Rp_mc
end

@userplot PBMTDLAS
@recipe function f(pt::PBMTDLAS)
    sol, dat, Rpp = pt.args
    l, Rp_mc = calc_lRp(sol)
    ld = calc_l(dat, sol.prob.p.hf0)
    Rpd = dat.Rp
    Rp_pikal = Rpp.(ld)
    # unitformat --> :square
    # xlabel --> "Dry layer thickness"
    # ylabel --> "Rp"
    @series begin
        c --> 1
        seriestype := :line
        label --> "PBM"
        return l, Rp_mc
    end
    @series begin
        c --> 2
        seriestype := :line
        label --> "3PP"
        return ld, Rp_pikal
    end
    @series begin
        c --> 3
        seriestype := :line
        label --> "TDLAS"
        return ld, Rpd
    end
end

# ------------------------

@userplot KineticPlot
@recipe function f(kp::KineticPlot; T = range(230, 250, length=101)*u"K")
    prm = kp.args[1]
    if prm isa ParamObjMC
        kf = prm.k_f
    elseif prm isa ODESolution
        kf = prm.prob.p.k_f
    end
    if kf isa DataInterpolations.AbstractInterpolation
        @series begin
            seriestype := :scatter
            primary := false
            return kf.t, kf.u
        end
    # elseif kf isa InterpK
    #     @series begin
    #         seriestype := :scatter
    #         primary := false
    #         return kf.spl.t, kf.spl.u
    #     end
    end
    @series begin
        seriestype := :samplemarkers
        yscale --> :log10
        label --> L"k_f"
        return T, kf.(T)
    end
end

# ----------------------------
using Distributions
@userplot BetaPlot
@recipe function f(bp::BetaPlot)
    sol = bp.args[1]
    ts = range(1, sol.t[end], length=5)
    for ti in ts
        hf = sol(ti, idxs=1)*u"cm"
        μ0 = sol(ti, idxs=3)
        μ2 = sol(ti, idxs=4)*u"mm^6" .|> u"mm^6"
        po = sol.prob.p
        hd = (po.hf0 .- hf) .|> u"cm"
        ϕ = (po.ρ_solution - po.c_solid)/po.ρ_solution
        μ1 = hd*ϕ*po.Ap .|> u"mm^3"
        
        # mean = μ1 ./ μ0
        # var = μ2 ./ μ0 .- mean.^2
        maxV = ϕ*po.hf0*po.Ap |> u"mm^3"
        mean = μ1/μ0/maxV
        var = μ2/μ0/maxV^2 - mean^2
        v = mean*(1-mean)/var - 1
        α = mean*v
        β = (1-mean)*v
        xs = logrange(1e-36, 1, length=100)
        vs = xs .* maxV
        βs = pdf(Beta(α, β), xs)
        @info "check" α β mean var βs
        @series begin
            label --> "t = $(round(ti, digits=2))"
            vs, βs
        end
    end
end