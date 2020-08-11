using ModelingToolkit, OrdinaryDiffEq, Test

@parameters t σ ρ β
@variables x(t) y(t) z(t) a(t)
@derivatives D'~t

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ a*y - β*z,
       0 ~ x - a]

lorenz1 = ODESystem(eqs,t,[x,y,z,a],[σ,ρ,β],name=:lorenz1)

eqs_auto, outputs_auto = ModelingToolkit.eliminate_aliases(lorenz1, [x,y,z])
@test isequal(eqs, lhss(eqs_staying) .~ substitute.(rhss(eqs_staying), (Dict(a=>x),)))
@test isequal(outputs, [a ~ x])

eqs_auto, outputs_auto = ModelingToolkit.eliminate_aliases(lorenz1)
@test isequal(eqs_auto, eqs)
@test isequal(outputs_auto, outputs_auto)
