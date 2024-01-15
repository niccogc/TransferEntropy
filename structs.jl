abstract type Kernels end
struct Frequency <: Kernels end
Freq = Frequency()
abstract type Dynamical_models end
struct Dyn_mod <: Dynamical_models end
Dyn = Dyn_mod()
abstract type TreeNode end
#f and g are the functions that evolve system x and y, a and b are the parameters of the system, x0 and y0 are the initial conditions
struct Dynamical_model <: Dynamical_models
    x::Vector{Vector{Float64}}
    y::Vector{Vector{Float64}}
    f::Function
    g::Function
    a::Array{Any,1}
    b::Array{Any,1}
    x0::Function
    y0::Function
end

# x and y are the collection of values that each variable can take in time, if the variables are representation of the same system then they will be only a vector
# n is the amount of time we observed the system
#X and Y are n realization of the dynamics
struct Time_Series
    n::Int64
    x::Vector{Vector{Float64}}
    y::Vector{Vector{Float64}}
    X::Vector{Vector{Float64}}
    Y::Vector{Vector{Float64}}
    Dynamic::Dynamical_models
end

struct TE
    k::Int64
    l::Int64
    n::Int64
    value::Float64
    t::Int64
    Dynamic::Dynamical_models
    Kernel::Kernels
    Prob::Matrix{Float64}
end

abstract type operators end

struct density_prob
    prob::Float64
    mat::Matrix{ComplexF64}
    result::Int64
end

struct pauli_operator <: operators
    mat::Matrix{ComplexF64}
    base::String
    evec::Vector{Vector{ComplexF64}}
    eval::Vector{ComplexF64}
end

struct channel <: operators
    mat::Matrix{ComplexF64}
    base::String
end

struct measure <: operators
    proj::Vector{Matrix{ComplexF64}}
    base::String
    value::Vector{Int64}
end

struct operator <: operators
    mat::Matrix{ComplexF64}
    base::String
    evec::Vector{Vector{ComplexF64}}
    eval::Vector{ComplexF64}
end
struct state_vector
    vec::Vector{ComplexF64}
    base::String
end

struct QT_Dynamical_model <: Dynamical_models
    x::Vector{Vector{Float64}}
    y::Vector{Vector{Float64}}
    state::state_vector
    operators::Vector{operators}
end

struct density_matrix <: operators
    mat::Matrix{ComplexF64}
    base::String
end

mutable struct Binary_tree <: TreeNode
    idx::Int64
    parent::Union{TreeNode, Nothing}
    value::density_prob
    left::Union{TreeNode, Nothing}
    right::Union{TreeNode, Nothing}
end

struct Tree
    nodes::Vector{TreeNode}
end