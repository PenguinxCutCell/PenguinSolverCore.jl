abstract type AbstractSystem end

"""
    rhs!(du, sys, u, p, t)

Semidiscrete residual evaluation hook. Physics packages are expected to implement this.
"""
function rhs!(du, sys::AbstractSystem, u, p, t)
    throw(MethodError(rhs!, (du, sys, u, p, t)))
end

"""
    mass_matrix(sys)

Optional mass matrix for semidiscretizations. Return `nothing` for explicit ODE form.
"""
mass_matrix(sys::AbstractSystem) = nothing

"""
    apply_updates!(sys, u, p, t)

Optional hook to apply smooth time dependence and post-update bookkeeping.
"""
apply_updates!(sys::AbstractSystem, u, p, t) = nothing

"""
    rebuild!(sys, u, p, t)

Optional hook to rebuild matrices/factorizations after discrete invalidations.
"""
rebuild!(sys::AbstractSystem, u, p, t) = nothing
