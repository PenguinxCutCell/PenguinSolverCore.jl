"""
    AbstractCouplingMode

Abstract supertype for coupling orchestration modes.
"""
abstract type AbstractCouplingMode end

function _check_unique_symbols(values::AbstractVector{Symbol}, context::AbstractString)
    length(values) == length(unique(values)) ||
        throw(ArgumentError("$context expects unique symbols, got $(collect(values))."))
    return values
end

"""
    OneWayCoupling(order)

Sequential, one-pass block orchestration. Blocks are advanced once in `order`
without outer fixed-point iterations.
"""
struct OneWayCoupling <: AbstractCouplingMode
    order::Vector{Symbol}

    function OneWayCoupling(order::Vector{Symbol})
        isempty(order) && throw(ArgumentError("OneWayCoupling order cannot be empty."))
        _check_unique_symbols(order, "OneWayCoupling")
        return new(order)
    end
end

OneWayCoupling(order::AbstractVector{Symbol}) = OneWayCoupling(collect(order))
OneWayCoupling(order::Symbol...) = OneWayCoupling(collect(order))

"""
    TwoWayCoupling(order; maxiter, atol, rtol, relaxation, norm_fields, sweep, verbose)

Outer Picard fixed-point coupling with optional per-block under-relaxation.
"""
struct TwoWayCoupling{T<:Real} <: AbstractCouplingMode
    order::Vector{Symbol}
    maxiter::Int
    atol::T
    rtol::T
    relaxation::Dict{Symbol,T}
    norm_fields::Vector{Symbol}
    sweep::Symbol
    verbose::Bool
end

function TwoWayCoupling(order::AbstractVector{Symbol};
                        maxiter::Integer=50,
                        atol::Real=1e-8,
                        rtol::Real=1e-6,
                        relaxation::AbstractDict{Symbol,<:Real}=Dict{Symbol,Float64}(),
                        norm_fields::AbstractVector{Symbol}=Symbol[],
                        sweep::Symbol=:GaussSeidel,
                        verbose::Bool=false)
    isempty(order) && throw(ArgumentError("TwoWayCoupling order cannot be empty."))
    _check_unique_symbols(order, "TwoWayCoupling")
    maxiter > 0 || throw(ArgumentError("TwoWayCoupling maxiter must be > 0."))
    sweep in (:GaussSeidel, :Jacobi) ||
        throw(ArgumentError("TwoWayCoupling sweep must be :GaussSeidel or :Jacobi, got `$sweep`."))

    T = Float64
    T = promote_type(T, typeof(float(atol)), typeof(float(rtol)))
    for ω in values(relaxation)
        T = promote_type(T, typeof(float(ω)))
    end

    relax = Dict{Symbol,T}()
    for (name, ω_raw) in pairs(relaxation)
        ω = convert(T, ω_raw)
        (zero(T) < ω <= one(T)) ||
            throw(ArgumentError("relaxation[$name] must satisfy 0 < ω <= 1, got $ω."))
        relax[name] = ω
    end

    norm_fields_vec = unique(collect(norm_fields))
    return TwoWayCoupling{T}(collect(order), Int(maxiter), convert(T, atol), convert(T, rtol),
                             relax, norm_fields_vec, sweep, verbose)
end

TwoWayCoupling(order::Symbol...; kwargs...) = TwoWayCoupling(collect(order); kwargs...)

"""
    CoupledBlock(name, model, state; cache=nothing, metadata=Dict())
    CoupledBlock(name, model; init=nothing, cache=nothing, metadata=Dict())

Container for one independently solvable block in a coupled orchestration.
"""
mutable struct CoupledBlock{M,S,C}
    name::Symbol
    model::M
    state::S
    cache::C
    metadata
end

function CoupledBlock(name::Symbol, model, state;
                      cache=nothing,
                      metadata=Dict{Symbol,Any}())
    return CoupledBlock{typeof(model),typeof(state),typeof(cache)}(name, model, state, cache, metadata)
end

function CoupledBlock(name::Symbol, model;
                      init=nothing,
                      cache=nothing,
                      metadata=Dict{Symbol,Any}())
    state = initialize_state(model, init)
    return CoupledBlock(name, model, state; cache=cache, metadata=metadata)
end

_identity_coupling_apply(data, from_block, to_block, problem) = data

"""
    CouplingMap(from, to, field; apply!=identity)

Directed coupling transfer from `from` block to `to` block for one named `field`.
`apply!` can transform data before injection.
"""
struct CouplingMap{F}
    from::Symbol
    to::Symbol
    field::Symbol
    apply!::F
end

function CouplingMap(from::Symbol, to::Symbol, field::Symbol; apply! = _identity_coupling_apply)
    return CouplingMap(from, to, field, apply!)
end

"""
    CoupledProblem(blocks, coupling_mode; maps=CouplingMap[])

Bundle of coupled blocks, orchestration mode, and directed coupling maps.
"""
struct CoupledProblem{B,M}
    blocks::Vector{B}
    coupling_mode::M
    maps::Vector{CouplingMap}
end

function CoupledProblem(blocks::AbstractVector{<:CoupledBlock},
                        coupling_mode::AbstractCouplingMode;
                        maps::AbstractVector{<:CouplingMap}=CouplingMap[])
    block_vec = CoupledBlock[blocks...]
    map_vec = CouplingMap[maps...]
    return CoupledProblem{CoupledBlock,typeof(coupling_mode)}(block_vec, coupling_mode, map_vec)
end
