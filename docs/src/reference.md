```@meta
CurrentModule = PenguinSolverCore
```

# API Reference

## Linear-System Core

```julia
LinearSystem{T}
```

Mutable container for linear solves:

- `A::SparseMatrixCSC{T,Int}`
- `b::Vector{T}`
- `x::Vector{T}`
- `cache`
- `last_dt::Union{Nothing,T}`

Constructors:

```julia
LinearSystem(A::SparseMatrixCSC{T,Int}, b::Vector{T}; x=nothing)
LinearSystem{T}(n::Int)
```

Assembly hooks:

```julia
assemble!(sys::LinearSystem, model, t, dt)
assemble_matrix!(sys::LinearSystem, model, t, dt, θ)
assemble_rhs!(sys::LinearSystem, model, uⁿ, t, dt, θ)
```

Solver entrypoints:

```julia
solve!(sys::LinearSystem; method::Symbol = :direct, reuse_factorization::Bool = true, kwargs...)
theta_step!(sys::LinearSystem, model, uⁿ, t, dt, θ; method=:direct, reuse_factorization=true, kwargs...)
```

## Generic Block Coupling

### Public Types

```julia
AbstractCouplingMode
OneWayCoupling(order)
TwoWayCoupling(order; maxiter, atol, rtol, relaxation, norm_fields, sweep, verbose)
CoupledBlock(name, model, state; cache=nothing, metadata=Dict())
CoupledBlock(name, model; init=nothing, cache=nothing, metadata=Dict())
CouplingMap(from, to, field; apply! = identity)
CoupledProblem(blocks, coupling_mode; maps=CouplingMap[])
```

### Public Coupled Solve Entrypoints

```julia
solve_coupled!(problem::CoupledProblem; return_history=false, store_history=false, kwargs...)
step_coupled!(problem::CoupledProblem, t, dt; return_history=false, store_history=false, kwargs...)
```

Behavior:

- `OneWayCoupling`: single ordered pass, with map transfers after each block solve.
- `TwoWayCoupling`: outer Picard loop with residual convergence checks and optional per-block under-relaxation.

If `return_history=true`, both functions return `(problem, history)` where `history` is a `CouplingHistory` record of outer-iteration residuals.

### Extension Hooks for Physics Packages

Implement these methods for your block/model types:

```julia
initialize_state(model, init)
advance_steady!(block::CoupledBlock; kwargs...)
advance_unsteady!(block::CoupledBlock, t, dt; kwargs...)
get_coupling_field(block::CoupledBlock, ::Val{field}) where field
set_coupling_field!(block::CoupledBlock, ::Val{field}, data) where field
copy_state(state)
coupling_residual(old_block::CoupledBlock, new_block::CoupledBlock, fields)
block_summary(block::CoupledBlock)
```

Fallback definitions throw `MethodError` so missing extensions are explicit.

### Coupling History

`CouplingHistory` stores:

- `iterations::Vector{Int}`
- `residuals::Vector{Float64}`
- `block_residuals::Vector{Dict{Symbol,Float64}}`

`TwoWayCoupling` convergence uses absolute and relative criteria:

```julia
global_residual <= atol + rtol * baseline
```

where `baseline` is the first outer residual.

## Canonical Docstrings

```@docs
AbstractCouplingMode
OneWayCoupling
TwoWayCoupling
CoupledBlock
CouplingMap
CoupledProblem
CouplingHistory
initialize_state
advance_steady!
advance_unsteady!
get_coupling_field
set_coupling_field!
copy_state
coupling_residual
block_summary
solve_coupled!
step_coupled!
apply_coupling_map!
apply_coupling_maps!
state_distance
underrelax!
compute_block_residual
compute_global_residual
```
