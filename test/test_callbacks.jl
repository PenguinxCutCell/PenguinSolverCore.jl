using Test
using PenguinSolverCore
import PenguinSolverCore: rhs!, rebuild!, update!

mutable struct ToySystem <: AbstractSystem
    val::Float64
    cache::InvalidationCache
    updates::UpdateManager
    rebuild_calls::Int
end

ToySystem(val::Float64 = 0.0) = ToySystem(val, InvalidationCache(), UpdateManager(), 0)

function rhs!(du, sys::ToySystem, u, p, t)
    du[1] = -u[1] + sys.val
    return du
end

rebuild!(sys::ToySystem, u, p, t) = (sys.rebuild_calls += 1)

mutable struct ThresholdUpdater <: AbstractUpdater
    trigger::Float64
    fired::Vector{Float64}
end

ThresholdUpdater(trigger::Float64) = ThresholdUpdater(trigger, Float64[])

function update!(upd::ThresholdUpdater, sys::ToySystem, u, p, t)
    push!(upd.fired, t)
    old = sys.val
    sys.val = t
    return (old <= upd.trigger && sys.val > upd.trigger) ? :matrix : :rhs_only
end

mutable struct PulseUpdater <: AbstractUpdater
    fired::Vector{Float64}
end

PulseUpdater() = PulseUpdater(Float64[])

function update!(upd::PulseUpdater, sys::ToySystem, u, p, t)
    push!(upd.fired, t)
    sys.val += 1.0
    return :matrix
end

@testset "Core schedules + invalidation" begin
    sys = ToySystem()
    at_upd = ThresholdUpdater(0.6)
    add_update!(sys, AtTimes([0.25, 0.50, 0.75]), at_upd)

    u = [0.0]
    for (k, t) in enumerate(0.0:0.25:1.0)
        apply_scheduled_updates!(sys, u, nothing, t; step = k - 1)
    end

    @test at_upd.fired == [0.25, 0.50, 0.75]
    @test sys.rebuild_calls == 1
    @test !sys.cache.dirty_matrix
    @test !sys.cache.dirty_geometry
end

@testset "Periodic schedule fires on cadence" begin
    sys = ToySystem()
    per_upd = PulseUpdater()
    add_update!(sys, Periodic(0.50, 0.0), per_upd)

    u = [0.0]
    for (k, t) in enumerate(0.0:0.25:1.0)
        apply_scheduled_updates!(sys, u, nothing, t; step = k - 1)
    end

    @test per_upd.fired == [0.0, 0.5, 1.0]
    @test sys.rebuild_calls == 3
end

mutable struct AlwaysMatrixUpdater <: AbstractUpdater
    fired::Vector{Float64}
end
AlwaysMatrixUpdater() = AlwaysMatrixUpdater(Float64[])
function update!(upd::AlwaysMatrixUpdater, sys::ToySystem, u, p, t)
    push!(upd.fired, t)
    sys.val += 1
    return :matrix
end

@testset "Rebuild contract: coalesce multiple matrix updates at same t" begin
    sys = ToySystem()
    u = [0.0]

    upd1 = AlwaysMatrixUpdater()
    upd2 = AlwaysMatrixUpdater()
    upd3 = AlwaysMatrixUpdater()

    add_update!(sys, AtTimes([0.5]), upd1)
    add_update!(sys, AtTimes([0.5]), upd2)
    add_update!(sys, AtTimes([0.5]), upd3)

    apply_scheduled_updates!(sys, u, nothing, 0.5; step=0)

    @test sys.rebuild_calls == 1
    @test upd1.fired == [0.5]
    @test upd2.fired == [0.5]
    @test upd3.fired == [0.5]
end

mutable struct AlwaysRhsUpdater <: AbstractUpdater
    fired::Vector{Float64}
end
AlwaysRhsUpdater() = AlwaysRhsUpdater(Float64[])
function update!(upd::AlwaysRhsUpdater, sys::ToySystem, u, p, t)
    push!(upd.fired, t)
    sys.val += 0.1
    return :rhs_only
end

@testset "Rebuild contract: rhs_only does not trigger rebuild; mixed triggers once" begin
    sys = ToySystem()
    u = [0.0]

    u_rhs = AlwaysRhsUpdater()
    u_mat = AlwaysMatrixUpdater()
    add_update!(sys, AtTimes([0.5]), u_rhs)
    add_update!(sys, AtTimes([0.5]), u_mat)

    apply_scheduled_updates!(sys, u, nothing, 0.5; step=0)

    @test sys.rebuild_calls == 1
    @test u_rhs.fired == [0.5]
    @test u_mat.fired == [0.5]
end

@testset "Rebuild contract: rhs_only never rebuilds" begin
    sys = ToySystem()
    u = [0.0]

    u_rhs1 = AlwaysRhsUpdater()
    u_rhs2 = AlwaysRhsUpdater()
    add_update!(sys, AtTimes([0.25, 0.5]), u_rhs1)
    add_update!(sys, Periodic(0.5, 0.0), u_rhs2)

    for (k,t) in enumerate(0.0:0.25:1.0)
        apply_scheduled_updates!(sys, u, nothing, t; step=k-1)
    end

    @test sys.rebuild_calls == 0
end

mutable struct GeometryUpdater <: AbstractUpdater
    fired::Vector{Float64}
end
GeometryUpdater() = GeometryUpdater(Float64[])
function update!(upd::GeometryUpdater, sys::ToySystem, u, p, t)
    push!(upd.fired, t)
    return :geometry
end

@testset "Rebuild contract: geometry invalidation triggers rebuild once" begin
    sys = ToySystem()
    u = [0.0]
    gupd = GeometryUpdater()
    add_update!(sys, AtTimes([0.5]), gupd)

    apply_scheduled_updates!(sys, u, nothing, 0.5; step=0)
    @test sys.rebuild_calls == 1
    @test gupd.fired == [0.5]
end