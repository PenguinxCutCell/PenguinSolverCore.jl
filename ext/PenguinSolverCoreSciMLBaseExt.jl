module PenguinSolverCoreSciMLBaseExt

using PenguinSolverCore
import SciMLBase

import PenguinSolverCore: sciml_odefunction, sciml_odeproblem, sciml_callbackset

const SCI_ML_ATOL = 100 * eps(Float64)

function sciml_odefunction(sys::PenguinSolverCore.AbstractSystem; jac = nothing, jac_prototype = nothing)
    f! = (du, u, p, t) -> PenguinSolverCore.rhs!(du, sys, u, p, t)

    kwargs = (mass_matrix = PenguinSolverCore.mass_matrix(sys),)
    jac === nothing || (kwargs = merge(kwargs, (; jac = jac)))
    jac_prototype === nothing || (kwargs = merge(kwargs, (; jac_prototype = jac_prototype)))

    return SciMLBase.ODEFunction(f!; kwargs...)
end

function _in_tspan(t, tspan; atol::Real = SCI_ML_ATOL)
    t0, tf = tspan
    if tf >= t0
        return (t >= t0 - atol) && (t <= tf + atol)
    end
    return (t <= t0 + atol) && (t >= tf - atol)
end

function _add_periodic_tstops!(integrator, schedule::PenguinSolverCore.Periodic)
    t0, tf = integrator.sol.prob.tspan
    (isfinite(t0) && isfinite(tf)) || return nothing
    tf >= t0 || return nothing

    t = schedule.t0
    while t < t0 - SCI_ML_ATOL
        t += schedule.dt
    end

    while t <= tf + SCI_ML_ATOL
        SciMLBase.add_tstop!(integrator, t)
        t += schedule.dt
    end

    return nothing
end

function _initialize_scheduled_tstops!(integrator, mgr::PenguinSolverCore.UpdateManager)
    tspan = integrator.sol.prob.tspan
    for ev in mgr.events
        schedule = ev.schedule
        if schedule isa PenguinSolverCore.AtTimes
            for t in schedule.ts
                _in_tspan(t, tspan) || continue
                SciMLBase.add_tstop!(integrator, t)
            end
        elseif schedule isa PenguinSolverCore.Periodic
            _add_periodic_tstops!(integrator, schedule)
        end
    end
    return nothing
end

function _combine_callbacks(core_cb, user_cb)
    core_cb === nothing && return user_cb
    user_cb === nothing && return core_cb
    return SciMLBase.CallbackSet(core_cb, user_cb)
end

function sciml_callbackset(sys::PenguinSolverCore.AbstractSystem; include_every_step::Bool = false, atol::Real = SCI_ML_ATOL)
    mgr = PenguinSolverCore.maybe_update_manager(sys)
    mgr === nothing && return nothing
    isempty(mgr.events) && return nothing

    callbacks = Any[]

    has_scheduled = any(ev -> (ev.schedule isa PenguinSolverCore.AtTimes || ev.schedule isa PenguinSolverCore.Periodic), mgr.events)
    if has_scheduled
        filter_scheduled = schedule -> (schedule isa PenguinSolverCore.AtTimes || schedule isa PenguinSolverCore.Periodic)

        condition = (u, t, integrator) -> PenguinSolverCore.has_due_updates(
            mgr,
            t,
            integrator.iter,
            integrator;
            schedule_filter = filter_scheduled,
            atol = atol,
        )

        affect! = integrator -> PenguinSolverCore.apply_scheduled_updates!(
            sys,
            mgr,
            integrator.u,
            integrator.p,
            integrator.t;
            step = integrator.iter,
            integrator = integrator,
            schedule_filter = filter_scheduled,
            atol = atol,
        )

        initialize = function (cb, u, t, integrator)
            _initialize_scheduled_tstops!(integrator, mgr)

            if PenguinSolverCore.has_due_updates(
                mgr,
                integrator.t,
                integrator.iter,
                integrator;
                schedule_filter = filter_scheduled,
                atol = atol,
            )
                PenguinSolverCore.apply_scheduled_updates!(
                    sys,
                    mgr,
                    integrator.u,
                    integrator.p,
                    integrator.t;
                    step = integrator.iter,
                    integrator = integrator,
                    schedule_filter = filter_scheduled,
                    atol = atol,
                )
            end
            return nothing
        end

        push!(callbacks, SciMLBase.DiscreteCallback(condition, affect!; initialize = initialize))
    end

    if include_every_step
        has_every_step = any(ev -> ev.schedule isa PenguinSolverCore.EveryStep, mgr.events)
        if has_every_step
            filter_every = schedule -> schedule isa PenguinSolverCore.EveryStep

            condition = (u, t, integrator) -> PenguinSolverCore.has_due_updates(
                mgr,
                t,
                integrator.iter,
                integrator;
                schedule_filter = filter_every,
                atol = atol,
            )

            affect! = integrator -> PenguinSolverCore.apply_scheduled_updates!(
                sys,
                mgr,
                integrator.u,
                integrator.p,
                integrator.t;
                step = integrator.iter,
                integrator = integrator,
                schedule_filter = filter_every,
                atol = atol,
            )

            push!(callbacks, SciMLBase.DiscreteCallback(condition, affect!))
        end
    end

    isempty(callbacks) && return nothing
    return length(callbacks) == 1 ? first(callbacks) : SciMLBase.CallbackSet(callbacks...)
end

function sciml_odeproblem(
    sys::PenguinSolverCore.AbstractSystem,
    u0,
    tspan;
    p = nothing,
    callback = nothing,
    include_every_step::Bool = false,
    kwargs...,
)
    f = sciml_odefunction(sys)
    core_cb = sciml_callbackset(sys; include_every_step = include_every_step)
    merged_cb = _combine_callbacks(core_cb, callback)

    if merged_cb === nothing
        return SciMLBase.ODEProblem(f, u0, tspan, p; kwargs...)
    end
    return SciMLBase.ODEProblem(f, u0, tspan, p; callback = merged_cb, kwargs...)
end

end
