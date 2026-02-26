```@meta
CurrentModule = PenguinSolverCore
```

# PenguinSolverCore.jl

## Installation

```julia
using Pkg
Pkg.add("PenguinSolverCore")
```

## Overview

`PenguinSolverCore.jl` provides shared solver plumbing for Penguin physics packages:

- Common system interface (`rhs!`, `mass_matrix`, `apply_updates!`, `rebuild!`)
- Solver-agnostic update scheduling (`AtTimes`, `Periodic`, `EveryStep`)
- Cache invalidation + rebuild triggering (`:rhs_only`, `:matrix`, `:geometry`)
- Optional SciML bridge for `ODEFunction`, `ODEProblem`, and callback mapping

## Quick Example (Core Loop)

```@example
using PenguinSolverCore
import PenguinSolverCore: rhs!, rebuild!, update!

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
    sys.val = 2.0
    return :matrix
end

sys = ToySystem()
upd = StepUpdater()
add_update!(sys, AtTimes([0.5]), upd)

u = [0.0]
for (k, t) in enumerate(0.0:0.25:1.0)
    apply_scheduled_updates!(sys, u, nothing, t; step = k - 1)
end

(sys.val, sys.rebuild_calls, upd.fired)
```

## SciML Bridge

If `SciMLBase` is installed, the package extension provides:

- `sciml_odefunction(sys; jac, jac_prototype)`
- `sciml_callbackset(sys; include_every_step = false)`
- `sciml_odeproblem(sys, u0, tspan; p, callback, include_every_step)`

Discrete update schedules are translated to `SciMLBase.DiscreteCallback`s and use `add_tstop!` for exact event times.

## Main Sections

- [Callbacks and Scheduling](callbacks.md)
- [API Reference](reference.md)
