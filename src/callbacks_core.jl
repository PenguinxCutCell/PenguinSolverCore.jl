abstract type AbstractUpdater end

"""
    update!(updater, sys, u, p, t) -> Symbol

Mutate `sys` (or `p`) in-place and return one of:
`:nothing`, `:rhs_only`, `:matrix`, or `:geometry`.
"""
function update!(upd::AbstractUpdater, sys, u, p, t)
    throw(MethodError(update!, (upd, sys, u, p, t)))
end

struct UpdateEvent{S<:AbstractSchedule,U<:AbstractUpdater}
    schedule::S
    updater::U
end

mutable struct UpdateManager
    events::Vector{UpdateEvent}
    last_fire::Vector{Float64}
    function UpdateManager(events::Vector{UpdateEvent} = UpdateEvent[])
        new(events, fill(-Inf, length(events)))
    end
end

function maybe_update_manager(sys)
    hasproperty(sys, :updates) || return nothing
    mgr = getproperty(sys, :updates)
    mgr isa UpdateManager || throw(ArgumentError("system `$(typeof(sys))` has `updates` field but it is not an UpdateManager"))
    return mgr
end

function update_manager(sys)
    mgr = maybe_update_manager(sys)
    mgr === nothing && throw(ArgumentError("system `$(typeof(sys))` has no `updates::UpdateManager` field"))
    return mgr
end

function add_update!(mgr::UpdateManager, schedule::AbstractSchedule, updater::AbstractUpdater)
    push!(mgr.events, UpdateEvent(schedule, updater))
    push!(mgr.last_fire, -Inf)
    return mgr
end

add_update!(sys, schedule::AbstractSchedule, updater::AbstractUpdater) = add_update!(update_manager(sys), schedule, updater)

function clear_updates!(mgr::UpdateManager)
    empty!(mgr.events)
    empty!(mgr.last_fire)
    return mgr
end

clear_updates!(sys) = clear_updates!(update_manager(sys))

function _event_should_fire(schedule::AbstractSchedule, last_fire::Float64, t, step::Int, integrator; atol::Real)
    if schedule isa Periodic
        return should_fire(schedule, t, step, integrator; last_fire = last_fire, atol = atol)
    elseif schedule isa AtTimes
        if isfinite(last_fire) && isapprox(last_fire, t; atol = atol, rtol = 0.0)
            return false
        end
        return should_fire(schedule, t, step, integrator; atol = atol)
    end
    return should_fire(schedule, t, step, integrator; atol = atol)
end

function has_due_updates(
    mgr::UpdateManager,
    t,
    step::Int = 0,
    integrator = nothing;
    schedule_filter::Function = _ -> true,
    atol::Real = DEFAULT_TIME_ATOL,
)
    for k in eachindex(mgr.events)
        ev = mgr.events[k]
        schedule_filter(ev.schedule) || continue
        _event_should_fire(ev.schedule, mgr.last_fire[k], t, step, integrator; atol = atol) && return true
    end
    return false
end

function _process_change!(cache::Union{InvalidationCache,Nothing}, change::Symbol)
    if change === :matrix
        cache === nothing || invalidate!(cache, Val(:matrix))
        return true
    elseif change === :geometry
        cache === nothing || invalidate!(cache, Val(:geometry))
        return true
    elseif change === :rhs_only || change === :nothing
        return false
    end
    throw(ArgumentError("unknown change flag `$change`; expected :rhs_only, :matrix, :geometry, or :nothing"))
end

function apply_scheduled_updates!(
    sys,
    mgr::UpdateManager,
    u,
    p,
    t;
    step::Int = 0,
    integrator = nothing,
    schedule_filter::Function = _ -> true,
    atol::Real = DEFAULT_TIME_ATOL,
)
    cache = maybe_invalidation_cache(sys)
    needs_local_rebuild = false

    for k in eachindex(mgr.events)
        ev = mgr.events[k]
        schedule_filter(ev.schedule) || continue
        _event_should_fire(ev.schedule, mgr.last_fire[k], t, step, integrator; atol = atol) || continue

        change = update!(ev.updater, sys, u, p, t)
        needs_local_rebuild |= _process_change!(cache, change)
        mgr.last_fire[k] = Float64(t)
    end

    apply_updates!(sys, u, p, t)

    needs_rebuild_now = needs_local_rebuild
    if cache !== nothing
        needs_rebuild_now |= needs_rebuild(cache)
    end

    if needs_rebuild_now
        rebuild!(sys, u, p, t)
    end

    cache === nothing || clear_invalidations!(cache)
    return nothing
end

apply_scheduled_updates!(sys, u, p, t; kwargs...) = apply_scheduled_updates!(sys, update_manager(sys), u, p, t; kwargs...)
