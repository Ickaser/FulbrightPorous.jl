export lyo_mc_dae!, lyo_mc_dae_f
export ParamObjMC


const R = Unitful.R
import LyoPronto: Mw, k_ice, ΔHsub
const η = LyoPronto.μ_vap # Water vapor viscosity
export R, ΔHsub, Mw, η

@doc raw"""
"""
function lyo_mc_dae!(du, u, params, t)

    _, hf0, c_solid, ρ_solution = params[1]
    Kshf, Av, Ap, = params[2]
    pch, Tsh = params[3] 
    k_f, C, N, S, R0 = params[4] 
    
    td = t*u"hr" # Dimensional time
    hf = u[1]*u"cm"
    Tf = u[2]*u"K"
    μ0 = u[3]
    μ2 = u[4]*u"cm^6"
    if Tf < 0u"K"
        return [NaN, NaN, NaN]
    end

    k = k_f(Tf)
    hd = hf0 - hf
    ϕ = (ρ_solution - c_solid)/ρ_solution
    μ1 = hd*ϕ*Ap
    pchl = pch(td)
    Qshf = Av*Kshf(pchl)*(Tsh(td) - Tf)
    Tsub = Tf - Qshf/k_ice/Ap*hf

    if μ1 <= 0u"cm^3" || μ0 <= 0.0 || μ2 <= 0u"cm^6"
        Rp = R0
    else
        μ32 = μ0^(-1//8)*μ1^(3//4)*μ2^(3//8)
        Rp = R0 + C * sqrt(R*Tf/Mw) * μ1*hd/μ32/ϕ # Knudsen
        # Rp = R0 + C*η*R*Tf/Mw/((calc_psub(Tsub)+pch(td))/2) * μ1*hd/μ2/ϕ # Viscous
    end
    dmdt = - Ap*(calc_psub(Tsub)-pch(td))/Rp
    Qsub = uconvert(u"W", dmdt*ΔHsub)

    dhf_dt = min(0u"cm/hr", dmdt/(ρ_solution-c_solid)/Ap) # Cap dhf_dt at 0: no desublimation

    du[1] = ustrip(u"cm/hr", dhf_dt)
    du[2] = ustrip(u"W", Qsub + Qshf)
    du[3] = ustrip(u"hr^-1", -k/2*μ0^2 - N(hd)*dhf_dt)
    du[4] = ustrip(u"cm^6/hr", k*μ1^2 - S(hd)*dhf_dt)

    return Rp
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

The DAE system which is the Pikal model plus pore moment evolution (2 ODEs, one nonlinear algebraic equation for pseudosteady conditions)
is here treated as a constant-mass-matrix implicit ODE system.
The implementation is in [`lyo_mc_dae!`](@ref)

The initial conditions `u0 = [h_f, Tf, μ2]` should be unitless, but are internally assigned to be in `[cm, K, cm^6]`.
The unitless time is taken to be in hours, so derivatives are given in unitless `[cm/hr, K/hr, cm^6/hr]`.

`params` has the form:
```
params = (
    (_, hf0, c_solid, ρ_solution),
    (Kshf, Av, Ap),
    (pch, Tsh) ,
    (k0, E, C, N, S, R0),
)
```
where those listed following are callables returning `Quantity`s, and the rest are `Quantity`s.
See [`RpFormFit`](@ref LyoPronto.RpFormFit) and [`RampedVariable`](@ref LyoPronto.RampedVariable) for convenience types that can help with the callables.
- `S(hd)` with `hd` a length returns source term for pore size second moment (as a Unitful quantity)
- `Kshf(p)` with `p` a pressure returns heat transfer coefficient (as a Unitful quantity).
- `Tsh(t)`, `pch(t)` return shelf temperature and chamber pressure respectively at time `t`.
"""
lyo_mc_dae_f = ODEFunction(lyo_mc_dae!, mass_matrix=lyo_mc_mm)


struct ParamObjMC
    _
    hf0
    c_solid
    ρ_solution
    Kshf
    Av
    Ap
    pch
    Tsh
    k_f
    C
    N
    S
    R0
end
function ParamObjMC(tuptup)
    return ParamObjMC(tuptup[1]..., tuptup[2]...,
        tuptup[3]..., tuptup[4]...)
end
Base.size(po::ParamObjMC) = (4,)
function Base.getindex(po::ParamObjMC, i)
    if i == 1
        return (po._, po.hf0, po.c_solid, po.ρ_solution)
    elseif i == 2
        return (po.Kshf, po.Av, po.Ap)
    elseif i == 3
        return (po.pch, po.Tsh)
    elseif i == 4
        return (po.k_f, po.C, po.N, po.S, po.R0)
    else
        error(BoundsError, "Attempt to access LyoPronto.ParamsObjRF at index $i. Only indices 1 to 4 allowed")
    end
end

function ODEProblem(po::ParamObjMC)
    return ODEProblem(lyo_mc_dae_f, 
    ustrip.([u"cm", u"K", u"cm^6", NoUnits], 
    [po.hf0, po.Tsh(0u"s"), 0.0u"cm^6", 0.0]),
    (0.0, 300.0), po; callback=end_drying_callback)
end
