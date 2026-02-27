# PenguinSolverCore.jl

[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://PenguinxCutCell.github.io/PenguinSolverCore.jl/dev)

Core solver plumbing shared by Penguin physics packages.

## Features

- Solver-agnostic update scheduling (`AtTimes`, `Periodic`, `EveryStep`)
- Updater contract with explicit invalidation semantics (`:rhs_only`, `:matrix`, `:geometry`)
- Coalesced rebuild path via `InvalidationCache` and `apply_scheduled_updates!`
- Optional SciML bridge (package extension) for `ODEFunction`, `ODEProblem`, and callbacks
- Support for mass matrices in SciML problems

Recommended updater convention:
- use `:rhs_only` for boundary/interface value changes (for example time-varying Dirichlet or Robin `g`)
- use `:matrix` when coefficients/factorizations change
- use `:geometry` when operator topology/geometry changes

## Installation

After cloning the repo,

```julia
using Pkg
Pkg.dev("PenguinSolverCore")
```

## Quick Example

```julia
using PenguinSolverCore

mutable struct ToySystem <: AbstractSystem
    val::Float64
    cache::InvalidationCache
    updates::UpdateManager
    rebuild_calls::Int
end

ToySystem() = ToySystem(0.0, InvalidationCache(), UpdateManager(), 0)
rhs!(du, sys::ToySystem, u, p, t) = (du[1] = -u[1] + sys.val)
rebuild!(sys::ToySystem, u, p, t) = (sys.rebuild_calls += 1)

mutable struct StepUpdater <: AbstractUpdater
    fired::Vector{Float64}
end
StepUpdater() = StepUpdater(Float64[])

function update!(upd::StepUpdater, sys::ToySystem, u, p, t)
    push!(upd.fired, t)
    sys.val = 1.0
    return :matrix
end

sys = ToySystem()
add_update!(sys, AtTimes([0.5]), StepUpdater())

u = [0.0]
for (k, t) in enumerate(0.0:0.25:1.0)
    apply_scheduled_updates!(sys, u, nothing, t; step = k - 1)
end
```

## Documentation

- Local build:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

- Dev docs site: <https://PenguinxCutCell.github.io/PenguinSolverCore.jl/dev>

## License

This package is part of the PenguinxCutCell project.
