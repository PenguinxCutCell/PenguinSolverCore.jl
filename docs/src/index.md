```@meta
CurrentModule = PenguinSolverCore
```

# PenguinSolverCore.jl

## Overview

`PenguinSolverCore.jl` currently focuses on linear-system assembly and solve infrastructure:

- `LinearSystem` as a mutable state container (`A`, `b`, `x`, solve cache)
- Assembly hooks for custom models (`assemble!`, `assemble_matrix!`, `assemble_rhs!`)
- Solve methods backed by direct LU or `LinearSolve.jl`
- `theta_step!` helper for time-stepping skeletons with `dt`-aware matrix assembly

## Recent Updates (March 2026)

- Package internals were streamlined to a lean linear-solver core.
- `solve!` now supports `method = :direct` or `method = :linearsolve`.
- `theta_step!` now tracks `last_dt` and only rebuilds matrix terms when needed.
- Documentation was rebuilt around the current API surface.

## Quick Example

```@example
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
x = theta_step!(sys, ToyModel(), [1.0, 3.0], 0.0, 0.1, 0.5)

x
```

## Main Sections

- [API Reference](reference.md)
- [Recent Updates](updates.md)
