# PenguinSolverCore.jl

[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://PenguinxCutCell.github.io/PenguinSolverCore.jl/dev)

`PenguinSolverCore.jl` provides reusable linear-system assembly and solve plumbing for PenguinxCutCell packages.

## Features

- `LinearSystem{T}` container for sparse matrix `A`, RHS `b`, solution `x`, and solve cache
- Extendable assembly hooks: `assemble!`, `assemble_matrix!`, `assemble_rhs!`
- Solve entrypoint with method selection:
  - `:direct` (LU factorization)
  - `:linearsolve` (via `LinearSolve.jl`)
- `theta_step!` helper with matrix reuse behavior when `dt` is unchanged

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

## Documentation

- Local build:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

- Dev docs site: <https://PenguinxCutCell.github.io/PenguinSolverCore.jl/dev>

## License

This package is part of the PenguinxCutCell project.
