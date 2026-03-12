```@meta
CurrentModule = PenguinSolverCore
```

# PenguinSolverCore.jl

## Overview

`PenguinSolverCore.jl` provides reusable solver-core infrastructure with two layers:

- Linear-system plumbing: `LinearSystem`, assembly hooks, `solve!`, `theta_step!`
- Generic block coupling orchestration:
  - `OneWayCoupling` for passive sequential field transfer
  - `TwoWayCoupling` for outer Picard fixed-point iterations

The package is intentionally physics-agnostic. It does not encode Darcy/Stokes/transport
model logic. Physics packages remain responsible for their own assembly and block solves,
and extend SolverCore via coupling hooks.

## Generic Block Coupling

Coupled orchestration is built around:

- `CoupledBlock`: one independently solvable block
- `CouplingMap`: named-field transfer (`:velocity`, `:temperature`, etc.)
- `CoupledProblem`: blocks + coupling mode + maps

Extension hooks to overload in physics packages:

- `advance_steady!` / `advance_unsteady!`
- `get_coupling_field` / `set_coupling_field!`
- `copy_state`
- optional `coupling_residual`

### Minimal Toy Example

```@example
using PenguinSolverCore

struct ProducerModel end
mutable struct ProducerState
    velocity::Float64
end

struct FollowerModel
    gain::Float64
end
mutable struct FollowerState
    velocity::Float64
    pressure::Float64
end

PenguinSolverCore.initialize_state(::ProducerModel, init) = ProducerState(0.0)
PenguinSolverCore.initialize_state(::FollowerModel, init) = FollowerState(0.0, 0.0)
PenguinSolverCore.copy_state(s::ProducerState) = ProducerState(s.velocity)
PenguinSolverCore.copy_state(s::FollowerState) = FollowerState(s.velocity, s.pressure)

function PenguinSolverCore.advance_steady!(block::CoupledBlock{M,S,C}; kwargs...) where {M<:ProducerModel,S<:ProducerState,C}
    block.state.velocity = 2.0
    return block
end

function PenguinSolverCore.advance_steady!(block::CoupledBlock{M,S,C}; kwargs...) where {M<:FollowerModel,S<:FollowerState,C}
    block.state.pressure = block.model.gain * block.state.velocity
    return block
end

PenguinSolverCore.get_coupling_field(block::CoupledBlock{M,S,C}, ::Val{:velocity}) where {M<:ProducerModel,S<:ProducerState,C} = block.state.velocity
function PenguinSolverCore.set_coupling_field!(block::CoupledBlock{M,S,C}, ::Val{:velocity}, data) where {M<:FollowerModel,S<:FollowerState,C}
    block.state.velocity = data
    return block
end

flow = CoupledBlock(:flow, ProducerModel(); init=nothing)
transport = CoupledBlock(:transport, FollowerModel(3.0); init=nothing)

problem = CoupledProblem(
    [flow, transport],
    OneWayCoupling([:flow, :transport]);
    maps=[CouplingMap(:flow, :transport, :velocity)],
)

solve_coupled!(problem)
transport.state.pressure
```

`TwoWayCoupling` uses the same hooks, but with an outer fixed-point loop, residual checks,
and optional under-relaxation.

## Main Sections

- [API Reference](reference.md)
- [Recent Updates](updates.md)
