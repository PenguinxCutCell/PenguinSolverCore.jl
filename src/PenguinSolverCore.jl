module PenguinSolverCore

using LinearAlgebra
using LinearSolve
using SparseArrays

export LinearSystem, assemble!, solve!, theta_step!, assemble_matrix!, assemble_rhs!
export CoupledBlock, CoupledProblem, CouplingMap
export OneWayCoupling, TwoWayCoupling
export solve_coupled!, step_coupled!
export get_coupling_field, set_coupling_field!
export advance_steady!, advance_unsteady!, coupling_residual

mutable struct LinearSystem{T}
    A::SparseMatrixCSC{T,Int}
    b::Vector{T}
    x::Vector{T}
    cache
    last_dt::Union{Nothing,T}
end

function LinearSystem(A::SparseMatrixCSC{T,Int}, b::Vector{T}; x::Union{Nothing,Vector{T}}=nothing) where {T}
    n = size(A, 2)
    length(b) == size(A, 1) || throw(DimensionMismatch("rhs length does not match matrix rows"))
    xvec = x === nothing ? zeros(T, n) : x
    length(xvec) == n || throw(DimensionMismatch("solution length does not match matrix cols"))
    return LinearSystem{T}(A, b, xvec, nothing, nothing)
end

function LinearSystem{T}(n::Int) where {T}
    A = spzeros(T, n, n)
    b = zeros(T, n)
    x = zeros(T, n)
    return LinearSystem{T}(A, b, x, nothing, nothing)
end

function assemble!(sys::LinearSystem, model, t, dt)
    throw(MethodError(assemble!, (sys, model, t, dt)))
end

function assemble_matrix!(sys::LinearSystem, model, t, dt, θ)
    throw(MethodError(assemble_matrix!, (sys, model, t, dt, θ)))
end

function assemble_rhs!(sys::LinearSystem, model, uⁿ, t, dt, θ)
    throw(MethodError(assemble_rhs!, (sys, model, uⁿ, t, dt, θ)))
end

function _solve_direct!(sys::LinearSystem{T}; reuse_factorization::Bool=true) where {T}
    if !reuse_factorization || !(sys.cache isa NamedTuple && get(sys.cache, :kind, nothing) === :direct)
        F = lu(sys.A)
        sys.cache = (kind=:direct, factor=F)
    end
    ldiv!(sys.x, sys.cache.factor, sys.b)
    return sys.x
end

function _solve_linearsolve!(sys::LinearSystem{T}; reuse_factorization::Bool=true, kwargs...) where {T}
    if !reuse_factorization || !(sys.cache isa NamedTuple && get(sys.cache, :kind, nothing) === :linearsolve)
        prob = LinearProblem(sys.A, sys.b)
        lcache = init(prob; kwargs...)
        sol = LinearSolve.solve!(lcache)
        sys.cache = (kind=:linearsolve, cache=lcache)
        copyto!(sys.x, sol.u)
        return sys.x
    end
    sol = LinearSolve.solve!(sys.cache.cache; b=sys.b)
    copyto!(sys.x, sol.u)
    return sys.x
end

function solve!(sys::LinearSystem; method::Symbol=:direct, reuse_factorization::Bool=true, kwargs...)
    if method === :direct
        return _solve_direct!(sys; reuse_factorization=reuse_factorization)
    elseif method === :linearsolve
        return _solve_linearsolve!(sys; reuse_factorization=reuse_factorization, kwargs...)
    end
    throw(ArgumentError("unknown solve method `$method`; expected :direct or :linearsolve"))
end

function theta_step!(sys::LinearSystem, model, uⁿ, t, dt, θ;
                     method::Symbol=:direct,
                     reuse_factorization::Bool=true,
                     kwargs...)
    dt_changed = sys.last_dt === nothing || sys.last_dt != dt

    if dt_changed && hasmethod(assemble_matrix!, Tuple{typeof(sys), typeof(model), typeof(t), typeof(dt), typeof(θ)})
        assemble_matrix!(sys, model, t, dt, θ)
    elseif dt_changed
        assemble!(sys, model, t, dt)
    end

    if hasmethod(assemble_rhs!, Tuple{typeof(sys), typeof(model), typeof(uⁿ), typeof(t), typeof(dt), typeof(θ)})
        assemble_rhs!(sys, model, uⁿ, t, dt, θ)
    elseif !dt_changed
        assemble!(sys, model, t, dt)
    end

    sys.last_dt = dt
    return solve!(sys; method=method, reuse_factorization=reuse_factorization, kwargs...)
end

include("coupling/types.jl")
include("coupling/interfaces.jl")
include("coupling/transfer.jl")
include("coupling/relaxation.jl")
include("coupling/history.jl")
include("coupling/one_way.jl")
include("coupling/two_way.jl")

end
