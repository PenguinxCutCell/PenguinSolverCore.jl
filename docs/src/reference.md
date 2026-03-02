```@meta
CurrentModule = PenguinSolverCore
```

# API Reference

## Core Type

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

## Assembly Hooks

Implement these methods for your model/system type:

```julia
assemble!(sys::LinearSystem, model, t, dt)
assemble_matrix!(sys::LinearSystem, model, t, dt, θ)
assemble_rhs!(sys::LinearSystem, model, uⁿ, t, dt, θ)
```

Default definitions throw `MethodError` to make missing implementations explicit.

## Solvers

```julia
solve!(sys::LinearSystem; method::Symbol = :direct, reuse_factorization::Bool = true, kwargs...)
```

- `method = :direct` uses `lu(sys.A)` and cached factorization reuse.
- `method = :linearsolve` uses `LinearSolve.jl` cache initialization/reuse.

Returns the updated solution vector `sys.x`.

## Theta Step Helper

```julia
theta_step!(sys::LinearSystem, model, uⁿ, t, dt, θ;
            method::Symbol = :direct,
            reuse_factorization::Bool = true,
            kwargs...)
```

Behavior:

- Reassembles matrix terms when `dt` changes (`assemble_matrix!` preferred, fallback `assemble!`).
- Reassembles RHS each step (`assemble_rhs!` preferred, fallback `assemble!` when matrix unchanged).
- Stores `sys.last_dt = dt` and then calls `solve!`.
