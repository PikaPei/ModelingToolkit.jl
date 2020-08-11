"""
```julia
calculate_tgrad(sys::AbstractSystem)
```

Calculate the time gradient of a system.

Returns a vector of [`Expression`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_tgrad end

"""
```julia
calculate_gradient(sys::AbstractSystem)
```

Calculate the gradient of a scalar system.

Returns a vector of [`Expression`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_gradient end

"""
```julia
calculate_jacobian(sys::AbstractSystem)
```

Calculate the jacobian matrix of a system.

Returns a matrix of [`Expression`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_jacobian end

"""
```julia
calculate_factorized_W(sys::AbstractSystem)
```

Calculate the factorized W-matrix of a system.

Returns a matrix of [`Expression`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_factorized_W end

"""
```julia
calculate_hessian(sys::AbstractSystem)
```

Calculate the hessian matrix of a scalar system.

Returns a matrix of [`Expression`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_hessian end

"""
```julia
generate_tgrad(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; kwargs...)
```

Generates a function for the time gradient of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_tgrad end

"""
```julia
generate_gradient(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; kwargs...)
```

Generates a function for the gradient of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_gradient end

"""
```julia
generate_jacobian(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; sparse = false, kwargs...)
```

Generates a function for the jacobian matrix matrix of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_jacobian end

"""
```julia
generate_factorized_W(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; sparse = false, kwargs...)
```

Generates a function for the factorized W-matrix matrix of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_factorized_W end

"""
```julia
generate_hessian(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; sparse = false, kwargs...)
```

Generates a function for the hessian matrix matrix of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_hessian end

"""
```julia
generate_function(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; kwargs...)
```

Generate a function to evaluate the system's equations.
"""
function generate_function end

function Base.getproperty(sys::AbstractSystem, name::Symbol)

    if name ∈ fieldnames(typeof(sys))
        return getfield(sys,name)
    elseif !isempty(sys.systems)
        i = findfirst(x->x.name==name,sys.systems)
        if i !== nothing
            return rename(sys.systems[i],renamespace(sys.name,name))
        end
    end

    i = findfirst(x->x.name==name,sys.states)
    if i !== nothing
        x = rename(sys.states[i],renamespace(sys.name,name))
        if :iv ∈ fieldnames(typeof(sys))
            return x(getfield(sys,:iv)())
        else
            return x()
        end
    end

    if :ps ∈ fieldnames(typeof(sys))
        i = findfirst(x->x.name==name,sys.ps)
        if i !== nothing
            return rename(sys.ps[i],renamespace(sys.name,name))()
        end
    end

    if :observed ∈ fieldnames(typeof(sys))
        i = findfirst(x->convert(Variable,x.lhs).name==name,sys.observed)
        if i !== nothing
            return rename(convert(Variable,sys.observed[i].lhs),renamespace(sys.name,name))(getfield(sys,:iv)())
        end
    end

    throw(error("Variable $name does not exist"))
end

renamespace(namespace,name) = Symbol(namespace,:₊,name)

function namespace_variables(sys::AbstractSystem)
    [rename(x,renamespace(sys.name,x.name)) for x in states(sys)]
end

function namespace_parameters(sys::AbstractSystem)
    [rename(x,renamespace(sys.name,x.name)) for x in parameters(sys)]
end

function namespace_pins(sys::AbstractSystem)
    [rename(x,renamespace(sys.name,x.name)) for x in pins(sys)]
end

namespace_equations(sys::AbstractSystem) = namespace_equation.(equations(sys),sys.name,sys.iv.name)


function namespace_equation(eq::Equation,name,ivname)
    _lhs = namespace_operation(eq.lhs,name,ivname)
    _rhs = namespace_operation(eq.rhs,name,ivname)
    _lhs ~ _rhs
end

function namespace_operation(O::Operation,name,ivname)
    if O.op isa Variable && O.op.name != ivname
        Operation(rename(O.op,renamespace(name,O.op.name)),namespace_operation.(O.args,name,ivname))
    else
        Operation(O.op,namespace_operation.(O.args,name,ivname))
    end
end
namespace_operation(O::Constant,name,ivname) = O

independent_variable(sys::AbstractSystem) = sys.iv
states(sys::AbstractSystem) = unique(isempty(sys.systems) ? setdiff(sys.states, convert.(Variable,sys.pins)) : [sys.states;reduce(vcat,namespace_variables.(sys.systems))])
parameters(sys::AbstractSystem) = isempty(sys.systems) ? sys.ps : [sys.ps;reduce(vcat,namespace_parameters.(sys.systems))]
pins(sys::AbstractSystem) = isempty(sys.systems) ? sys.pins : [sys.pins;reduce(vcat,namespace_pins.(sys.systems))]
function observed(sys::AbstractSystem; keep=nothing)
    allequations(sys; keep=keep, diffeqs=false, alias_zero_lhs=false)
end

function states(sys::AbstractSystem,name::Symbol)
    x = sys.states[findfirst(x->x.name==name,sys.states)]
    rename(x,renamespace(sys.name,x.name))(sys.iv())
end

function parameters(sys::AbstractSystem,name::Symbol)
    x = sys.ps[findfirst(x->x.name==name,sys.ps)]
    rename(x,renamespace(sys.name,x.name))()
end

function pins(sys::AbstractSystem,name::Symbol)
    x = sys.pins[findfirst(x->x.name==name,sys.ps)]
    rename(x,renamespace(sys.name,x.name))(sys.iv())
end

lhss(xs) = map(x->x.lhs, xs)
rhss(xs) = map(x->x.rhs, xs)

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
    isempty(idxs) && return sys.eqs, sys.observed

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

function make_lhs_0(eq)
    eq.lhs isa Constant && iszero(eq.lhs) && return eq
    0 ~ eq.lhs - eq.rhs
end

function allequations(sys::ModelingToolkit.AbstractSystem; keep=nothing, diffeqs=true, aliases=true, alias_zero_lhs=true)
    myeqs, outputs = eliminate_aliases(sys, keep)
    eqs = [myeqs;
           reduce(vcat,
                  namespace_equations.(sys.systems);
                  init=Equation[])]

    vcat(diffeqs ? eqs : [],
         aliases ? (alias_zero_lhs ? make_lhs_0.(outputs) : outputs) : [])
end

function equations(sys::ModelingToolkit.AbstractSystem; keep=nothing)
    allequations(sys; keep=keep, aliases=false)
end

function equations(sys::ModelingToolkit.AbstractSystem; keep=nothing)
    allequations(sys; keep=keep, aliases=false)
end

function states(sys::AbstractSystem,args...)
    name = last(args)
    extra_names = reduce(Symbol,[Symbol(:₊,x.name) for x in args[1:end-1]])
    newname = renamespace(extra_names,name)
    rename(x,renamespace(sys.name,newname))(sys.iv())
end

function parameters(sys::AbstractSystem,args...)
    name = last(args)
    extra_names = reduce(Symbol,[Symbol(:₊,x.name) for x in args[1:end-1]])
    newname = renamespace(extra_names,name)
    rename(x,renamespace(sys.name,newname))()
end

function islinear(sys::AbstractSystem)
    rhs = [eq.rhs for eq ∈ equations(sys)]

    iv = sys.iv()
    dvs = [dv(iv) for dv ∈ states(sys)]

    all(islinear(r, dvs) for r in rhs)
end

function pins(sys::AbstractSystem,args...)
    name = last(args)
    extra_names = reduce(Symbol,[Symbol(:₊,x.name) for x in args[1:end-1]])
    newname = renamespace(extra_names,name)
    rename(x,renamespace(sys.name,newname))(sys.iv())
end

struct AbstractSysToExpr
    sys::AbstractSystem
    states::Vector{Variable}
end
AbstractSysToExpr(sys) = AbstractSysToExpr(sys,states(sys))
function (f::AbstractSysToExpr)(O::Operation)
    any(isequal(O), f.states) && return O.op.name  # variables
    if isa(O.op, Variable)
        isempty(O.args) && return O.op.name  # 0-ary parameters
        return build_expr(:call, Any[O.op.name; f.(O.args)])
    end
    return build_expr(:call, Any[O.op; f.(O.args)])
end
(f::AbstractSysToExpr)(x) = convert(Expr, x)
