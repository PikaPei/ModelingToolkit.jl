using ModelingToolkit, OrdinaryDiffEq, Test

@parameters t σ ρ β
@variables x(t) y(t) z(t) a(t)
@derivatives D'~t

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ a*y - β*z,
       0 ~ x - a]

lorenz1 = ODESystem(eqs,t,[x,y,z,a],[σ,ρ,β],name=:lorenz1)

import ModelingToolkit: Constant, get_variables, lhss, rhss

# This decides what will be the final internal state of the system.
function get_nonaliases(sys)
    get_variables.(filter(!iszero, lhss(sys.eqs))) |> Iterators.flatten |> collect
end

function eliminate_aliases(sys, keep=nothing)
    if keep === nothing
        keep = get_nonaliases(sys)
    end
    eqs = sys.eqs
    idxs = findall(x->x.lhs isa Constant && iszero(x.lhs), eqs)
    isempty(idxs) && return sys

    to_eliminate = []
    for i in idxs
        append!(to_eliminate, filter(x->!any(isequal(x), keep), get_variables(eqs[i])))
    end
    subs = solve_for(eqs[idxs], to_eliminate)
    subdict = Dict(lhss(subs) .=> rhss(subs))
    eqs′ = eqs[setdiff(1:length(eqs), idxs)]
    eqs′ = lhss(eqs′) .~ substitute.(rhss(eqs′), (subdict,))
    return eqs′, subs
end

eqs, outputs = eliminate_aliases(lorenz1, [x,y,z])
eqs_staying = lorenz1.eqs[1:end-1]

@test isequal(eqs, lhss(eqs_staying) .~ substitute.(rhss(eqs_staying), (Dict(a=>x),)))
@test isequal(outputs, [a ~ x])

eqs_auto, outputs_auto = eliminate_aliases(lorenz1)
@test isequal(eqs_auto, eqs)
@test isequal(outputs_auto, outputs_auto)

using Test


hasanyvars(expr, vars) = any(occursin.(vars, (expr,)))

function solve_for(eqs, vars)
    @assert length(eqs) >= length(vars)
    @assert all(iszero(eq.lhs) for eq in eqs)
    neweqs = []
    for (i, eq) in enumerate(eqs)
        rhs = eq.rhs
        if rhs.op == (-)
            if any(isequal(rhs.args[1]), vars) && any(isequal(rhs.args[2]), vars)
                push!(neweqs, rhs.args[1] ~ rhs.args[2]) # pick one?
                warn("todo")
            elseif any(isequal(rhs.args[1]), vars)
                push!(neweqs, rhs.args[1] ~ rhs.args[2])
            elseif any(isequal(rhs.args[2]), vars)
                push!(neweqs, rhs.args[2] ~ rhs.args[1])
            else
                error("todo 2")
            end
        elseif rhs.op == (+)
            eqs[i] = 0 ~ rhs.args[1] - (-rhs.args[2])
        else
            error("todo")
        end
    end
    if length(neweqs) >= length(vars)
        return neweqs
    else
        # substitute
        eqs′ = Equation.(0, substitute.(rhss(eqs), (Dict(lhss(neweqs) .=> rhss(neweqs),))))
        solve_for(eqs′, vars)
    end
end

