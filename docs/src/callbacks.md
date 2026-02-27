```@meta
CurrentModule = PenguinSolverCore
```

# Callbacks and Scheduling

## Update Model

An updater is any type that implements:

```julia
update!(updater, sys, u, p, t) -> Symbol
```

Allowed return symbols:

- `:nothing`  : no change
- `:rhs_only` : RHS-only change, no rebuild
- `:matrix`   : matrix/factorization invalidation
- `:geometry` : geometry invalidation (rebuild required)

Use `UpdateManager` plus `add_update!` to register events.

### Recommended conventions

For consistency across Penguin physics packages:

- Box Dirichlet value updates (including time-varying `Ref` payloads): `:rhs_only`
- Interface RHS updates such as Robin/jump `g(x,t)`: `:rhs_only`
- Interface coefficient updates that change constraint matrices (`a`, `b`, `b1`, `b2`, `α1`, `α2`): `:matrix`
- Geometry/moment/operator topology changes: `:geometry`

## Schedules

- `AtTimes(ts)` : fire at specific times (internally normalized/sorted)
- `Periodic(dt, t0)` : fire from `t0` every `dt`
- `EveryStep()` : fire every step (core loop; optional in SciML mode)

## Core Driver

`apply_scheduled_updates!` evaluates schedules, executes updaters, applies invalidation flags, and calls `rebuild!` exactly once per call if needed.

```julia
apply_scheduled_updates!(sys, mgr, u, p, t; step = 0)
```

## SciML Behavior

With the `SciMLBase` extension loaded:

- `AtTimes` and `Periodic` are mapped to `DiscreteCallback`
- callback initialization adds `tstops` for exact stepping at event times
- due updates at integrator start are also applied during callback initialization
- user callback can be composed with the core callback via `CallbackSet`
