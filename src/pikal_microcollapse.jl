export lyo_mc_dae!, lyo_mc_dae_f
export ParamObjMC


const R = Unitful.R
import LyoPronto: Mw, k_ice, ΔHsub
const η = LyoPronto.μ_vap # Water vapor viscosity
export R, ΔHsub, Mw, η

_getval(x) = hasproperty(x, :value) ? x.value : x

"""
    $(SIGNATURES)
    Calculate the mass transfer, heat transfer, and other quantities for the microcollapse-augmented Pikal model.
"""
@inline function calc_md_Q_mc(u, params, t)

    (; hf0, csolid, ρsolution,
    Kshf, Av, Ap,
    pch, Tsh, 
    k_f, C, Rs, ) = params
    
    td = t*u"hr" # Dimensional time
    hf = u[1]*u"cm"
    Tf = u[2]*u"K"
    μ0 = u[3]
    μ2 = u[4]*u"mm^6"

    hd = hf0 - hf
    ϕ = (ρsolution - csolid)/ρsolution
    μ1 = hd*ϕ*Ap
    Q_shf = Av*Kshf(pch(td))*(Tsh(td) - Tf)
    Tsub = Tf - Q_shf/k_ice/Ap*hf
    k = k_f(Tsub)

    porefac = μ1^2*μ0/μ2^3 # Knudsen
    # porefac = μ1/μ2 # Viscous but double check this

    # Compute Rp, checking for negative pore structure first
    # # Direct comparison with 0 interacts weirdly with ForwardDiff: cares about derivatives, not only value
    # # if _getval(μ1) <= 0u"cm^3" || _getval(μ0) <= 0.0 || _getval(μ2) <= 0u"mm^6" 
    # if isnan(porefac) || (hasproperty(porefac, :value) ? porefac.value : porefac) <= 0.0u"mm^-12" # Manually check for ForwardDiff
    #     Rp = 0.0u"m/s" # No pores to resist flow
    # else
    #     Rp = C*sqrt(R*Tsub/Mw)/ϕ*hd * sqrt(sqrt(sqrt(porefac))) |> u"m/s" 
    #     # Rp = C*η*R*Tf/Mw/((calc_psub(Tsub)+pch(td))/2) * hd/ϕ*porefac # Viscous but double check
    # end

    # Compute Rp with NaNs for aphysical inputs, then check for NaN 
    # This happens when μ0, μ1, or μ2 are negative or zero
    Rp = C*NaNMath.sqrt(R*Tsub/Mw)/ϕ*hd * sqrt(sqrt(NaNMath.sqrt(porefac))) |> u"m/s" 
    if !isfinite(Rp)
        Rp = 0.0u"m/s" # No pores to resist flow
    end

    Rsp = Rs + Rp
    md = - Ap*(LyoPronto.calc_psub(Tsub)-pch(td))/Rsp
    return (; md, Q_shf, k, Tsub, Rp, Rsp, μ1)
end

"""
    $(SIGNATURES)
"""
function lyo_mc_dae!(du, u, params, t)
    hf = u[1]*u"cm"
    Tf = u[2]*u"K"
    μ0 = u[3]
    if Tf < 0u"K"# || Tsub < 0u"K"
        du .= NaN
    end
    (; md, Q_shf, k, μ1) = calc_md_Q_mc(u, params, t)
    (; csolid, ρsolution, Ap, hf0, N, S) = params

    hd = hf0 - hf
    
    Q_sub = uconvert(u"W", md*ΔHsub)

    dhf_dt = min(0.0u"cm/hr", md/(ρsolution-csolid)/Ap |>u"cm/hr" ) # Cap dhf_dt at 0: no desublimation
    dμ0_dt = -k/2*μ0^2 - N(hd)*Ap*dhf_dt
    dμ2_dt = k*μ1^2 - S(hd)*Ap*dhf_dt # dhf_dt < 0 for drying proceeding

    du[1] = ustrip(u"cm/hr", dhf_dt)
    du[2] = ustrip(u"W", Q_sub + Q_shf)
    du[3] = ustrip(u"hr^-1", dμ0_dt)
    du[4] = ustrip(u"mm^6/hr", dμ2_dt)

    return 
end

const lyo_mc_mm = Diagonal([1.0, 0, 1, 1])

@doc raw"""
    lyo_mc_dae_f = ODEFunction(lyo_mc_dae!, mass_matrix=lyo_mc_mm)

Compute the right hand side function for the microcollapse-augmented Pikal model.

The DAE system which is the Pikal model plus pore moment evolution (3 ODEs, one nonlinear algebraic equation for pseudosteady conditions)
is here treated as a constant-mass-matrix implicit ODE system.
The implementation is in [`lyo_mc_dae!`](@ref)

The initial conditions `u0 = [h_f, Tf, μ0, μ2]` should be unitless, but are internally assigned to be in `[cm, K, NoUnits, mm^6]`.
The unitless time is taken to be in hours, so derivatives are given in unitless `[cm/hr, K/hr, 1/hr, mm^6/hr]`.

`params` is a ParamObjMC:
```
params = ParamObjMC((
    (hf0, csolid, ρsolution),
    (Kshf, Av, Ap),
    (pch, Tsh) ,
    (k_f, C, N, S, Rs)  
))
```
where those listed following are callables returning `Quantity`s, and the rest are `Quantity`s.
See [`RpFormFit`](@ref LyoPronto.RpFormFit) and [`RampedVariable`](@ref LyoPronto.RampedVariable) for convenience types that can help with the callables.
- `N(hd)` with `hd` a length returns source term for number of pores per volume (as a Unitful quantity)
- `S(hd)` with `hd` a length returns source term for pore size second moment per volume (as a Unitful quantity)
- `Kshf(p)` with `p` a pressure returns heat transfer coefficient (as a Unitful quantity).
- `Tsh(t)`, `pch(t)` return shelf temperature and chamber pressure respectively at time `t`.
- `k_f(T)` returns the rate constant for pore coalescence at temperature `T`.
"""
const lyo_mc_dae_f = ODEFunction{true, SciMLBase.FullSpecialize}(lyo_mc_dae!, mass_matrix=lyo_mc_mm)

@concrete terse struct ParamObjMC <: LyoPronto.ParamObj
    hf0
    csolid
    ρsolution
    Kshf
    Av
    Ap
    pch
    Tsh
    k_f
    C
    N
    S
    Rs
end

function ParamObjMC(tuptup)
    return ParamObjMC(tuptup[1]..., tuptup[2]...,
        tuptup[3]..., tuptup[4]...)
end
# Base.size(::ParamObjMC) = (4,)
Base.length(::ParamObjMC) = 14
# function Base.getindex(po::ParamObjMC, i)
#     if i == 1
#         return (po.hf0, po.csolid, po.ρsolution)
#     elseif i == 2
#         return (po.Kshf, po.Av, po.Ap)
#     elseif i == 3
#         return (po.pch, po.Tsh)
#     elseif i == 4
#         return (po.k_f, po.C, po.N, po.S, po.Rs)
#     else
#         error(BoundsError, "Attempt to access LyoPronto.ParamsObjRF at index $i. Only indices 1 to 4 allowed")
#     end
# end
function Base.show(io::IO, po::ParamObjMC) 
    compact = get(io, :compact, false)::Bool
    if compact
        return print(io, "ParamObjMC(...)")
    end
    str = "ParamObjMC( 
    $(po.hf0), 
    $(po.csolid), 
    $(po.ρsolution), 
    $(po.Kshf), 
    $(po.Av), 
    $(po.Ap), 
    $(po.pch), 
    $(po.Tsh),
    $(po.k_f), 
    $(po.C), 
    $(po.N), 
    $(po.S), 
    $(po.Rs) 
    )"
    return print(io, str)
end

function LyoPronto.calc_u0(po::ParamObjMC)
    # ustrip.([u"cm", u"K", u"mm^6", NoUnits], [po.hf0, float(po.Tsh(0u"s")), 0.0u"mm^6", 0.0])
    [ustrip(u"cm", po.hf0), ustrip(u"K", float(po.Tsh(0u"s"))), 0.0, 0.0]
end
LyoPronto.get_tstops(po::ParamObjMC) = get_tstops((po.Tsh, po.pch))
LyoPronto.get_t0(po::ParamObjMC) = LyoPronto.get_t0(po.Tsh, po.pch)*1.001

function ODEProblem(po::ParamObjMC; u0=calc_u0(po), tspan=(0.0, 200.0))
    tstops = get_tstops(po)
    t0 = LyoPronto.get_t0(po)
    tspan = @reset tspan[1] = t0
    return ODEProblem{true}(lyo_mc_dae_f, u0, tspan, po; 
        tstops=tstops, dt=0.1,
        callback=end_drying_callback, initializealg=BrownFullBasicInit())
end
