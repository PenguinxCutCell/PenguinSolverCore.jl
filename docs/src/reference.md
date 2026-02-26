```@meta
CurrentModule = PenguinSolverCore
```

# API Reference

## Core System Hooks

```@docs
rhs!
mass_matrix
apply_updates!
rebuild!
```

## Updater Hook

```@docs
update!
```

## Exported Types and Utilities

- `AbstractSystem`
- `AbstractUpdater`
- `AbstractSchedule`
- `EveryStep`, `AtTimes`, `Periodic`
- `UpdateEvent`, `UpdateManager`
- `normalize_times`, `should_fire`
- `add_update!`, `clear_updates!`, `has_due_updates`
- `InvalidationCache`, `invalidate!`, `clear_invalidations!`, `needs_rebuild`
- `apply_scheduled_updates!`
- `sciml_odefunction`, `sciml_callbackset`, `sciml_odeproblem` (requires `SciMLBase` extension)

Use `?name` in the Julia REPL for available inline docs and method signatures.
