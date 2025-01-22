export lyo_mc_dae!, lyo_mc_dae_f
export ParamObjMC


const R = Unitful.R
import LyoPronto: Mw, k_ice, ΔHsub
const η = LyoPronto.μ_vap # Water vapor viscosity
export R, ΔHsub, Mw, η

@doc raw"""
"""
function lyo_mc_dae!(du, u, params, t)

    hf0, c_solid, ρ_solution = params[1]
    Kshf, Av, Ap, = params[2]
    pch, Tsh = params[3] 
    k_f, C, N, S, Rs = params[4] 
    
    td = t*u"hr" # Dimensional time
    hf = u[1]*u"cm"
    Tf = u[2]*u"K"
    μ0 = u[3]
    μ2 = u[4]*u"cm^6"

    # exit hatch for bad parameter guesses
    if Tf < 0u"K"
        return fill(NaN, 4)
    end

    k = k_f(Tf)
    hd = hf0 - hf
    ϕ = (ρ_solution - c_solid)/ρ_solution
    μ1 = hd*ϕ*Ap
    Qshf = Av*Kshf(pch(td))*(Tsh(td) - Tf)
    Tsub = Tf - Qshf/k_ice/Ap*hf

    if μ1 < 0u"cm^3" || μ0 < 0.0 || μ2 <= 0u"cm^6"
        Rp = 0u"m/s" # No pores to resist flow
    else
        # μ32 = μ0^(-1//8)*μ1^(3//4)*μ2^(3//8)
        # Rp = R0 + C * sqrt(R*Tf/Mw) * μ1*hd/μ32/ϕ # Knudsen
        Rp = C*sqrt(R*Tf/Mw)/ϕ*hd * sqrt(sqrt(sqrt(μ1^2*μ0/μ2^3))) # Knudsen, rewritten to avoid runtime dispatch
        # Rp = R0 + C*η*R*Tf/Mw/((calc_psub(Tsub)+pch(td))/2) * μ1*hd/μ2/ϕ # Viscous
    end
    Rsp = Rs + Rp
    dmdt = - Ap*(calc_psub(Tsub)-pch(td))/Rsp
    Qsub = uconvert(u"W", dmdt*ΔHsub)

    dhf_dt = min(0u"cm/hr", dmdt/(ρ_solution-c_solid)/Ap) # Cap dhf_dt at 0: no desublimation
    dμ0_dt = -k/2*μ0^2 - N(hd)*Ap*dhf_dt
    dμ2_dt = k*μ1^2 - S(hd)*Ap*dhf_dt # dhf_dt < 0 for drying proceeding

    du[1] = ustrip(u"cm/hr", dhf_dt)
    du[2] = ustrip(u"W", Qsub + Qshf)
    du[3] = ustrip(u"hr^-1", dμ0_dt)
    du[4] = ustrip(u"cm^6/hr", dμ2_dt)

    return Rsp
end

const lyo_mc_mm = [1.0 0.0 0.0 0.0; 
                   0.0 0.0 0.0 0.0; 
                   0.0 0.0 1.0 0.0;
                   0.0 0.0 0.0 1.0]
# lyo_mc_mm2 = [1.0 0.0 0.0 0.0; 
#                    0.0 0.0 0.0 0.0; 
#                    0.0 0.0 1.0 0.0;
#                    0.0 0.0 0.0 1.0]

@doc raw"""
    lyo_mc_dae_f = ODEFunction(lyo_mc_dae!, mass_matrix=lyo_mc_mm)

Compute the right hand side function for the microcollapse-augmented Pikal model.

The DAE system which is the Pikal model plus pore moment evolution (3 ODEs, one nonlinear algebraic equation for pseudosteady conditions)
is here treated as a constant-mass-matrix implicit ODE system.
The implementation is in [`lyo_mc_dae!`](@ref)

The initial conditions `u0 = [h_f, Tf, μ0, μ2]` should be unitless, but are internally assigned to be in `[cm, K, NoUnits, cm^6]`.
The unitless time is taken to be in hours, so derivatives are given in unitless `[cm/hr, K/hr, 1/hr, cm^6/hr]`.

`params` has the form:
```
params = (
    (hf0, c_solid, ρ_solution),
    (Kshf, Av, Ap),
    (pch, Tsh) ,
    (k_f, C, N, S, Rx)  
)
```
where those listed following are callables returning `Quantity`s, and the rest are `Quantity`s.
See [`RpFormFit`](@ref LyoPronto.RpFormFit) and [`RampedVariable`](@ref LyoPronto.RampedVariable) for convenience types that can help with the callables.
- `N(hd)` with `hd` a length returns source term for number of pores per height (as a Unitful quantity)
- `S(hd)` with `hd` a length returns source term for pore size second moment per height (as a Unitful quantity)
- `Kshf(p)` with `p` a pressure returns heat transfer coefficient (as a Unitful quantity).
- `Tsh(t)`, `pch(t)` return shelf temperature and chamber pressure respectively at time `t`.
- `k_f(T)` returns the rate constant for pore coalescence at temperature `T`.
"""
lyo_mc_dae_f = ODEFunction(lyo_mc_dae!, mass_matrix=lyo_mc_mm)

struct ParamObjMC{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13}
    hf0::T1
    c_solid::T2
    ρ_solution::T3
    Kshf::T4
    Av::T5
    Ap::T6
    pch::T7
    Tsh::T8
    k_f::T9
    C::T10
    N::T11
    S::T12
    Rs::T13
end
function ParamObjMC(tuptup)
    return ParamObjMC(tuptup[1]..., tuptup[2]...,
        tuptup[3]..., tuptup[4]...)
end
Base.size(po::ParamObjMC) = (4,)
Base.length(po::ParamObjMC) = 14
function Base.getindex(po::ParamObjMC, i)
    if i == 1
        return (po.hf0, po.c_solid, po.ρ_solution)
    elseif i == 2
        return (po.Kshf, po.Av, po.Ap)
    elseif i == 3
        return (po.pch, po.Tsh)
    elseif i == 4
        return (po.k_f, po.C, po.N, po.S, po.Rs)
    else
        error(BoundsError, "Attempt to access LyoPronto.ParamsObjRF at index $i. Only indices 1 to 4 allowed")
    end
end
function Base.show(io::IO, po::ParamObjMC) 
    str = "ParamObjMC( 
    $(po.hf0), 
    $(po.c_solid), 
    $(po.ρ_solution), 
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

function ODEProblem(po::ParamObjMC; u0=nothing, tspan=nothing)
    if isnothing(tspan)
        tspan = (0.0, 200.0)
    end
    if isnothing(u0)
        u0 = ustrip.([u"cm", u"K", u"cm^6", NoUnits], [po.hf0, po.Tsh(0u"s"), 0.0u"cm^6", 0.0])
    end
    tstops = [0.0]
    for control in [po.Tsh, po.pch]
        if control isa RampedVariable && !isnothing(control.timestops)
            tstops = vcat(tstops, ustrip.(u"hr", control.timestops))
        elseif control isa LinearInterpolation
            tstops = vcat(tstops, ustrip.(u"hr", control.t))
        end
    end
    return ODEProblem(lyo_mc_dae_f, u0, tspan, po; 
        tstops=tstops, callback=end_drying_callback)
end
