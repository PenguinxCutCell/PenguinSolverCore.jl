module PenguinSolverCore

include("interfaces.jl")
include("cache.jl")
include("schedule.jl")
include("callbacks_core.jl")
include("state.jl")
include("dofs.jl")

export AbstractSystem,
    rhs!,
    mass_matrix,
    apply_updates!,
    rebuild!,
    AbstractUpdater,
    update!,
    AbstractSchedule,
    EveryStep,
    AtTimes,
    Periodic,
    normalize_times,
    should_fire,
    UpdateEvent,
    UpdateManager,
    maybe_update_manager,
    update_manager,
    add_update!,
    clear_updates!,
    has_due_updates,
    InvalidationCache,
    invalidate!,
    clear_invalidations!,
    needs_rebuild,
    maybe_invalidation_cache,
    apply_scheduled_updates!,
    StateView,
    pack_state,
    unpack_state,
    DofMap,
    restrict,
    prolong!,
    sciml_odefunction,
    sciml_odeproblem,
    sciml_callbackset

function sciml_odefunction(args...; kwargs...)
    throw(ArgumentError("`sciml_odefunction` requires SciMLBase. Add `SciMLBase` to the active environment to enable the extension."))
end

function sciml_odeproblem(args...; kwargs...)
    throw(ArgumentError("`sciml_odeproblem` requires SciMLBase. Add `SciMLBase` to the active environment to enable the extension."))
end

function sciml_callbackset(args...; kwargs...)
    throw(ArgumentError("`sciml_callbackset` requires SciMLBase. Add `SciMLBase` to the active environment to enable the extension."))
end

end
