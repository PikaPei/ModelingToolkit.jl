using ModelingToolkit, OrdinaryDiffEq, Test

test_equal(a, b) = @test isequal(simplify(a), simplify(b))

@parameters t σ ρ β F(t)
@variables x(t) y(t) z(t) a(t) u(t)
@derivatives D'~t

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ a*y - β*z,
       0 ~ x - a]

lorenz1 = ODESystem(eqs,t,[x,y,z,a],[σ,ρ,β],name=:lorenz1)

eqs_auto, outputs_auto = ModelingToolkit.eliminate_aliases(lorenz1, [x,y,z])

eqs_staying = lorenz1.eqs[1:end-1]
test_equal.(eqs_auto, lhss(eqs_staying) .~ substitute.(rhss(eqs_staying), (Dict(a=>x),)))
@test isequal(outputs, [a ~ x])

@test length(equations(lorenz1)) == 3
# @test length(states(lorenz1)) == 3

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

test_equal.(equations(lorenz1), eqs)

#=
test_lorenz1_aliased = ODESystem(eqs,t,[x,y,z],[σ,ρ,β],observed=[a ~ x],name=:lorenz1)

@test isequal(ODESystem(equations(lorenz1), observed=observed(lorenz1), name=:lorenz1),
              ODESystem(eqs,t,[x,y,z],[σ,ρ,β],observed=[a ~ x],name=:lorenz1))
=#
# Multi-System Reduction

eqs1 = [D(x) ~ σ*(y-x) + F,
       D(y) ~ x*(ρ-z)-u,
       D(z) ~ x*y - β*z]
aliases = [u ~ x + y - z]
lorenz1 = ODESystem(eqs1,pins=[F],observed=aliases,name=:lorenz1)

eqs2 = [D(x) ~ F,
       D(y) ~ x*(ρ-z)-x,
       D(z) ~ x*y - β*z]
aliases2 = [u ~ x - y - z]
lorenz2 = ODESystem(eqs2,pins=[F],observed=aliases2,name=:lorenz2)

connections = [lorenz1.F ~ lorenz2.u,
               lorenz2.F ~ lorenz1.u]
connected = ODESystem([0 ~ a + lorenz1.x - lorenz2.y],t,[a],[],observed=connections,systems=[lorenz1,lorenz2])

# Reduced Unflattened System

connections2 = [lorenz1.F ~ lorenz2.u,
                lorenz2.F ~ lorenz1.u,
                a ~ -lorenz1.x + lorenz2.y]
connected = ODESystem(Equation[],t,[],[],observed=connections2,systems=[lorenz1,lorenz2])

# Reduced Flattened System

flatten(sys::ModelingToolkit.AbstractSystem) = ODESystem(ModelingToolkit.allequations(sys),
                                                         independent_variable(sys),
                                                         states(sys),
                                                         parameters(sys))
flattened_system = flatten(connected)

#states(flattened_system) == [
#        lorenz1.x
#        lorenz1.y
#        lorenz1.z
#        lorenz2.x
#        lorenz2.y
#        lorenz2.z
#]

# parameters(reduced_flattened_system) == [
#         lorenz1.σ
#         lorenz1.ρ
#         lorenz1.β
#         lorenz2.σ
#         lorenz2.ρ
#         lorenz2.β
# ]

equations(flattened_system) == [
       D(lorenz1.x) ~ lorenz1.σ*(lorenz1.y-lorenz1.x) + lorenz2.x - lorenz2.y - lorenz2.z,
       D(lorenz1.y) ~ lorenz1.x*(ρ-z)-lorenz1.x - lorenz1.y + lorenz1.z,
       D(lorenz1.z) ~ lorenz1.x*lorenz1.y - lorenz1.β*lorenz1.z,
       D(lorenz2.x) ~ lorenz1.x + lorenz1.y - lorenz1.z,
       D(lorenz2.y) ~ lorenz2.x*(lorenz2.ρ-lorenz2.z)-lorenz2.x,
       D(lorenz2.z) ~ lorenz2.x*lorenz2.y - lorenz2.β*lorenz2.z]

observed(reduced_flattened_system) == [
        lorenz1.F ~ lorenz2.u
        lorenz2.F ~ lorenz1.u
        a ~ -lorenz1.x + lorenz2.y
]
