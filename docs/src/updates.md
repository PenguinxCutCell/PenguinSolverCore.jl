# Recent Updates

## March 2026

- Reworked package scope around sparse linear-system workflows.
- Added explicit solver backend selection in `solve!`:
  - `:direct` using LU factorization and optional factorization reuse.
  - `:linearsolve` using `LinearSolve.jl` cache initialization and reuse.
- Added `theta_step!` orchestration with:
  - matrix assembly updates when `dt` changes,
  - RHS-only updates every step,
  - stored `last_dt` for reuse decisions.
- Added a generic, physics-agnostic block coupling layer:
  - `OneWayCoupling` for passive sequential coupling.
  - `TwoWayCoupling` for Picard outer iterations with optional under-relaxation.
  - Named-field transfer via `CouplingMap` and extension hooks.
- Added toy-only tests and docs examples for coupling orchestration.
