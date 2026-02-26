const DEFAULT_TIME_ATOL = 100 * eps(Float64)

abstract type AbstractSchedule end

struct EveryStep <: AbstractSchedule end

struct AtTimes{T<:Real} <: AbstractSchedule
    ts::Vector{T}
    function AtTimes{T}(ts::Vector{T}) where {T<:Real}
        new{T}(normalize_times(ts))
    end
end

AtTimes(ts::AbstractVector{T}) where {T<:Real} = AtTimes{T}(collect(ts))
AtTimes(ts::Real...) = AtTimes(collect(ts))

struct Periodic{T<:Real} <: AbstractSchedule
    dt::T
    t0::T
    function Periodic{T}(dt::T, t0::T) where {T<:Real}
        dt > zero(T) || throw(ArgumentError("`dt` must be > 0"))
        new{T}(dt, t0)
    end
end

Periodic(dt::T, t0::T = zero(T)) where {T<:Real} = Periodic{T}(dt, t0)

function normalize_times(ts::AbstractVector{T}) where {T<:Real}
    out = collect(ts)
    isempty(out) && return out
    sort!(out)
    unique!(out)
    return out
end

should_fire(::EveryStep, t, step, integrator_or_nothing; kwargs...) = true

function should_fire(schedule::AtTimes, t, step, integrator_or_nothing; atol::Real = DEFAULT_TIME_ATOL)
    isempty(schedule.ts) && return false

    idx = searchsortedfirst(schedule.ts, t)
    if idx <= length(schedule.ts) && isapprox(t, schedule.ts[idx]; atol = atol, rtol = 0.0)
        return true
    end
    if idx > 1 && isapprox(t, schedule.ts[idx - 1]; atol = atol, rtol = 0.0)
        return true
    end
    return false
end

function should_fire(schedule::Periodic, t, step, integrator_or_nothing; last_fire::Real = -Inf, atol::Real = DEFAULT_TIME_ATOL)
    if t < schedule.t0 - atol
        return false
    end
    if !isfinite(last_fire)
        return t >= schedule.t0 - atol
    end
    return (t - last_fire) >= (schedule.dt - atol)
end
