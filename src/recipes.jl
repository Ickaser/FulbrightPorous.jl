export blankplot_hrC
blankplot_hrC(;kwargs...) = plot(u"hr", u"°C", xlabel="Time", ylabel="Temperature", unitformat=:square;kwargs...)

# ------------------------
# Visualize the pore size evolution
export blankplot_lRp
blankplot_lRp(;kwargs...) = plot(u"cm", u"cm^2*Torr*hr/g", xlabel=L"h_d", ylabel=L"R_p", unitformat=(x,y)->"$x [$(latexify(y))]";kwargs...)

public calc_moms
public calc_moms_NS

@userplot PoreEvolPlot
@recipe function f(pep::PoreEvolPlot; kind=:moms)
    sol = pep.args[1]
    if kind != :moms && kind != :stats
        @warn "kind should be either :moms or :stats, defaulting to :moms" kind
        kind = :moms
    end
    t = sol.t*u"hr"
    hf = sol[1,:]*u"cm"
    μ0 = sol[3,:]
    μ2 = sol[4,:]*u"mm^6" .|> u"mm^6"
    po = sol.prob.p
    hd = (po.hf0 .- hf) .|> u"cm"
    ϕ = (po.ρsolution - po.csolid)/po.ρsolution
    μ1 = hd*ϕ*po.Ap .|> u"mm^3"
    mean = μ1 ./ μ0
    var = μ2 ./ μ0 .- mean.^2
    Rp = map(sol.t, sol.u) do t,u
        lyo_mc_dae!(fill(0.0, 4), u, sol.prob.p, t, Val(true)) |> u"cm^2*Torr*hr/g"
    end
    seriestype --> :samplemarkers
    layout --> @layout [[a ; b; c;] d]
    @series begin
        subplot --> 1
        xlabel --> ""
        legend --> :none
        xformatter -->  _->""
        if kind == :stats
            ylabel --> "number\nof pores" #*", number of pores"
        elseif kind == :moms
            ylabel --> "μ0"
        end
        return hd, μ0
    end
    @series begin
        subplot --> 2
        xlabel --> ""
        legend --> :none
        xformatter -->  _->""
        if kind == :stats
            ylabel --> "typical \npore size\n" #*", total pore volume"
            return hd, mean
        elseif kind == :moms
            ylabel --> "μ1" 
            return hd, μ1
        end
    end
    @series begin
        subplot --> 3
        xlabel --> "Dry layer thickness"
        legend --> :none
        if kind == :stats
            ylabel --> "pore size \nvariance\n" #*", pore size second moment"
            return hd, var
        elseif kind == :moms
            ylabel --> "μ2"
            return hd, μ2
        end
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
    po = sol.prob.p
    hf = sol[1,:]*u"cm"
    hd = (po.hf0 .- hf) .|> u"cm"
    Rp = map(sol.t, sol.u) do ti,ui
        lyo_mc_dae!(fill(0.0, 4), ui, sol.prob.p, ti, Val(true)) |> u"m/s"
    end
    @series begin
        ylabel --> "Rp"
        label --> "\$R_p\$"
        hd, Rp
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
    Rp_mc = [FulbrightPores.lyo_mc_dae!(fill(0.0, 4), u, sol.prob.p, t, Val(true)) for (u,t) in zip(sol.u, sol.t)]
    return l, Rp_mc
end

@userplot PBMTDLAS
@recipe function f(pt::PBMTDLAS)
    pbmsol, dat, Rpp = pt.args
    l, Rp_mc = calc_lRp(pbmsol)
    ld = calc_l(dat, pbmsol.prob.p.hf0)
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
        seriestype := :samplemarkers
        marker --> :o
        step --> length(Rpd) ÷ 20
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
    elseif prm isa InterpK || prm isa Sigmoid_k || prm isa MYEGA_k
        kf = prm
    else
        error("Kinetic plot expects either a ParamObjMC or ODESolution")
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
function calc_moms(u, po)
    hf = u[1]*u"cm"
    μ0 = u[3]
    μ2 = u[4]*u"mm^6" .|> u"mm^6"
    hd = (po.hf0 .- hf) .|> u"cm"
    ϕ = 1 - po.csolid/po.ρsolution
    μ1 = hd*ϕ*po.Ap .|> u"mm^3"
    mean = μ1/μ0
    var = μ2/μ0  - mean^2
    return mean, var, po
end
function calc_moms_NS(po)
    ϕ = 1 - po.csolid/po.ρsolution
    mean = ϕ/po.N(0) 
    var = (po.S(0)/po.N(0) - ϕ^2/po.N(0)^2)
    return mean, var, po
end
function make_dist(mean, var, po, disttype)
    ϕ = 1 - po.csolid/po.ρsolution
    maxV = ϕ*po.hf0*po.Ap |> u"mm^3"
    mean = mean/maxV |> NoUnits
    var = var/maxV^2 |> NoUnits
    if disttype == :gamma
        θ = var/mean
        α = mean/θ
        return Gamma(α, θ), maxV
    elseif disttype == :normal
        return Truncated(Normal(mean, sqrt(var)), 0, Inf), maxV
    elseif disttype == :beta
        v = mean*(1-mean)/var - 1
        α = mean*v
        β = (1-mean)*v
        return Beta(α, β), maxV
    else
        error("disttype should be :gamma, :normal, or :beta")
    end
end
@userplot DistPlot
@recipe function f(dp::DistPlot; disttype=:normal)
    mean, var, po = dp.args
    dist, maxV = make_dist(mean, var, po, disttype)
    xs = logrange(max(1e-16, quantile(dist, 1e-4)), quantile(dist, 1-1e-4), length=100)
    vs = xs .* maxV
    ps = pdf(dist, xs)
    @series begin
        vs, ps
    end
end
function DistPlot(tupsol::Tuple{ODESolution})
    sol = tupsol[1]
    po = sol.prob.p
    mean, var, po = calc_moms(sol[end], po)
    return DistPlot((mean, var, po))
end

@userplot TimeDistPlot
@recipe function f(tdp::TimeDistPlot; disttype=:normal)
    sol = tdp.args[1]
    po = sol.prob.p

    # COmpute new-pore distribution
    mean, var, po = calc_moms_NS(po)
    @series begin
        disttype := disttype
        DistPlot((mean, var, po))
    end

    # Compute pore distribution at given time points
    ts = range(1, sol.t[end], length=5)
    for ti in ts
        mean, var, po = calc_moms(sol(ti), po)
        @series begin
            label --> "t = $(round(ti, digits=2))"
            disttype := disttype
            DistPlot((mean, var, po))
        end
    end
end