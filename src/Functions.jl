module NDFunctions

using LinearAlgebra
using SparseArrays

import Base.convert
import Base.promote_rule

#=
This module needs to be documented. Together with Constructors
it forms the backbone of the core API, users should provide us
with Arrays of VertexFunction and EdgeFunction as well as a graph and that's it.
=#

export StaticVertex
export StaticEdge
export ODEVertex
export ODEEdge
export VertexFunction
export EdgeFunction
export DDEVertex
export DDEEdge

#=All the structs contain the field sym, this has not been used yet in the main
constructor network_dynamics. Shouldn't be too hard though. More of a design question. StaticVertex
not implemented yet. =#

@Base.kwdef struct StaticVertex{F}
    f!::F # (v, e_s, e_t, p, t) -> nothing
    dim::Int # number of dimensions of x
    sym=[:v for i in 1:dim]::Array{Symbol, 1} # Symbols for the dimensions
end


@Base.kwdef struct StaticEdge{F}
    f!::F # (e, v_s, v_t, p, t) -> nothing
    dim::Int # number of dimensions of x
    sym::Array{Symbol, 1} = [:e for i in 1:dim] # Symbols for the dimensions
end


@Base.kwdef struct ODEVertex{F, T}
    f!::F # The function with signature (dx, x, e_s, e_t, p, t) -> nothing
    dim::Int # number of dimensions of x
    mass_matrix::T = I # Mass matrix for the equation
    sym::Array{Symbol, 1} = [:v for i in 1:dim] # Symbols for the dimensions
end


@Base.kwdef struct ODEEdge{F, T}
    f!::F # The function with signature (dx, x, e_s, e_t, p, t) -> nothing
    dim::Int # number of dimensions of x
    mass_matrix::T = I # Mass matrix for the equation
    sym=[:e for i in 1:dim] # Symbols for the dimensions
end


@Base.kwdef struct DDEVertex{F, T}
    f!::F # The function with signature (dv, v, h_v, e_s, e_t, h_e_s, h_e_d, p, t) -> nothing where h is the history function
    dim::Int # number of dimensions
    mass_matrix::T = I # Mass matrix for the equation
    sym=[:v for i in 1:dim] # Symbols for the dimensions
    tau_s=[] # Array of Delays for the incoming currents of different variables
    tau_d=[] # Array of Delays for the outgoing currents of different variables
end


@Base.kwdef struct DDEEdge{F, T}
    f!::F # The function with signature (de, e, h_e, v_s, v_d, h_s, h_d, p, t)
    dim::Int # number of variables
    mass_matrix::T = I # Mass matrix for the equation
    sym::Array{Symbol, 1} = [:v for i in 1:dim] # Symbols for the dimensions
end

# =========================================================
# =========================================================

# Now for some type magic and automatic promotion stuff...

struct ODE_from_Static{F}
    f!::F
end
function (ofs::ODE_from_Static)(dx,x,args...)
    # If mass matrix = 0 the differential equation sets dx = 0.
    # To set x to the value calculated by f! we first write the value calculated
    # by f! into dx, then subtract x. This leads to the  constraint
    # 0 = - x + f(...)
    # where f(...) denotes the value that f!(a, ...) writes into a.
    ofs.f!(dx,args...)
    dx .-= x
    nothing
end

function ODEVertex(sv::StaticVertex)
    ODEVertex(ODE_from_Static(sv.f!), sv.dim, 0., sv.sym)
end

function ODEEdge(se::StaticEdge)
    ODEEdge(ODE_from_Static(se.f!), se.dim, 0., se.sym)
end


const VertexFunction = Union{ODEVertex, StaticVertex, DDEVertex}
const EdgeFunction = Union{ODEEdge, StaticEdge, DDEEdge}

convert(::Type{ODEVertex{ODE_from_Static{F}, T}}, x::StaticVertex{T}) where F where T = ODEVertex(x)
promote_rule(::Type{ODEVertex{ODE_from_Static{F}, T}}, ::Type{StaticVertex{F}}) where F where T = ODEVertex{ODE_from_Static{F} ,T}
promote_rule(::Type{ODEVertex}, ::Type{StaticVertex}) = ODEVertex

# investigate whether convert(::Type(ODEEdge), se::StaticEdge) = ODEEdge(se)
# is useful, check out PowerDynamics for what it does...



end
