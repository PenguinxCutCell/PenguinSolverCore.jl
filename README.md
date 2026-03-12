# PenguinSolverCore.jl

[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://PenguinxCutCell.github.io/PenguinSolverCore.jl/dev)
![CI](https://github.com/PenguinxCutCell/PenguinSolverCore.jl/workflows/CI/badge.svg)
[![codecov](https://codecov.io/gh/PenguinxCutCell/PenguinSolverCore.jl/graph/badge.svg?token=Q50hXtjAKk)](https://codecov.io/gh/PenguinxCutCell/PenguinSolverCore.jl)

`PenguinSolverCore.jl` provides reusable solver orchestration primitives for PenguinxCutCell packages.
It now includes:
- linear-system assembly/solve plumbing
- generic multiphysics block coupling (`OneWayCoupling` and `TwoWayCoupling`)

## Features

- `LinearSystem{T}` container for sparse matrix `A`, RHS `b`, solution `x`, and solve cache
- Extendable assembly hooks: `assemble!`, `assemble_matrix!`, `assemble_rhs!`
- Solve entrypoint with method selection:
  - `:direct` (LU factorization)
  - `:linearsolve` (via `LinearSolve.jl`)
- `theta_step!` helper with matrix reuse behavior when `dt` is unchanged
- Physics-agnostic coupling orchestration:
  - `OneWayCoupling` for passive sequential coupling
  - `TwoWayCoupling` for outer Picard fixed-point coupling
  - Named-field transfer via `CouplingMap` + extension hooks

## Recent Updates (March 2026)

- Refactored package API around `LinearSystem` + assembly/solve hooks.
- Added `solve!` backend selection (`:direct` and `:linearsolve`) and cache reuse controls.
- Added `theta_step!` skeleton to support theta-method stepping with split matrix/RHS assembly.
- Simplified docs layout and moved API details into `docs/src/reference.md`.

## Installation

After cloning the repo,

```julia
using Pkg
Pkg.dev("PenguinSolverCore")
```

## Quick Example

```julia
using SparseArrays
using PenguinSolverCore

struct ToyModel end

function PenguinSolverCore.assemble_matrix!(sys::LinearSystem, ::ToyModel, t, dt, θ)
    n = length(sys.b)
    sys.A = spdiagm(0 => fill(2.0 + dt + θ, n))
    return sys
end

function PenguinSolverCore.assemble_rhs!(sys::LinearSystem, ::ToyModel, uⁿ, t, dt, θ)
    sys.b .= uⁿ .+ t .+ dt .+ θ
    return sys
end

sys = LinearSystem(spdiagm(0 => ones(2)), zeros(2))
model = ToyModel()
u0 = [1.0, 3.0]

x = theta_step!(sys, model, u0, 0.0, 0.1, 0.5; method = :direct)
```

## Generic Block Coupling

`PenguinSolverCore.jl` remains physics-agnostic: it does not assemble Darcy/Stokes/transport equations.
Physics packages provide block behavior by overloading hooks such as:

- `advance_steady!` / `advance_unsteady!`
- `get_coupling_field` / `set_coupling_field!`
- `copy_state` (for outer-loop residual snapshots)

Minimal toy example:

```julia
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
# transport.state.pressure == 6.0
```

For mutually coupled blocks, replace `OneWayCoupling(...)` with:

```julia
TwoWayCoupling([:A, :B]; maxiter=50, atol=1e-8, rtol=1e-6, relaxation=Dict(:A=>0.7))
```

and call `solve_coupled!` (steady) or `step_coupled!` (unsteady).

## Documentation

- Local build:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

- Dev docs site: <https://PenguinxCutCell.github.io/PenguinSolverCore.jl/dev>

## License

This package is part of the PenguinxCutCell project.
