using Test
using PenguinSolverCore
import PenguinSolverCore: rhs!, mass_matrix, rebuild!, update!

if Base.find_package("SciMLBase") === nothing || Base.find_package("OrdinaryDiffEq") === nothing
    @testset "SciML extension (guarded)" begin
        @test true
    end
else
    using LinearAlgebra
    using SciMLBase
    using OrdinaryDiffEq

    mutable struct SciToySystem <: AbstractSystem
        val::Float64
        cache::InvalidationCache
        updates::UpdateManager
        rebuild_calls::Int
    end

    SciToySystem(val::Float64 = 0.0) = SciToySystem(val, InvalidationCache(), UpdateManager(), 0)

    function rhs!(du, sys::SciToySystem, u, p, t)
        du[1] = -u[1] + sys.val
        return du
    end

    mass_matrix(sys::SciToySystem) = Diagonal([1.0])

    rebuild!(sys::SciToySystem, u, p, t) = (sys.rebuild_calls += 1)

    mutable struct StepChangeUpdater <: AbstractUpdater
        fired::Vector{Float64}
    end

    StepChangeUpdater() = StepChangeUpdater(Float64[])

    function update!(upd::StepChangeUpdater, sys::SciToySystem, u, p, t)
        push!(upd.fired, t)
        sys.val = 2.0
        return :matrix
    end

    mutable struct SciPulseUpdater <: AbstractUpdater
        fired::Vector{Float64}
    end

    SciPulseUpdater() = SciPulseUpdater(Float64[])

    function update!(upd::SciPulseUpdater, sys::SciToySystem, u, p, t)
        push!(upd.fired, t)
        sys.val += 1.0
        return :matrix
    end

    @testset "SciML ODE problem + callbacks" begin
        sys = SciToySystem(0.0)
        upd = StepChangeUpdater()
        add_update!(sys, AtTimes([0.5]), upd)

        prob = sciml_odeproblem(sys, [0.0], (0.0, 1.0); p = nothing)
        sol = SciMLBase.solve(prob, OrdinaryDiffEq.Tsit5(); reltol = 1e-9, abstol = 1e-9)

        expected = 2.0 * (1.0 - exp(-0.5))
        @test sys.val == 2.0
        @test sys.rebuild_calls == 1
        @test length(upd.fired) == 1
        @test isapprox(upd.fired[1], 0.5; atol=10eps(0.5))
        @test sol.u[end][1] ≈ expected atol = 5e-4
    end

    @testset "SciML: AtTimes enforces exact tstop + fires once" begin
        sys = SciToySystem(0.0)
        tstar = 0.37
        upd = StepChangeUpdater()
        add_update!(sys, AtTimes([tstar]), upd)

        prob = sciml_odeproblem(sys, [0.0], (0.0, 1.0); p = nothing)

        sol = SciMLBase.solve(prob, OrdinaryDiffEq.Tsit5();
            reltol = 1e-10, abstol = 1e-10,
            save_everystep = true,
        )

        @test length(upd.fired) == 1
        @test isapprox(upd.fired[1], tstar; atol = 1000eps(tstar))
        # "Lock" that solver actually visited tstar (tstop)
        @test any(t -> isapprox(t, tstar; atol = 1000eps(tstar)), sol.t)
        @test sys.rebuild_calls == 1
    end

    @testset "SciML: AtTimes unsorted times all fire exactly once" begin
        sys = SciToySystem(0.0)
        upd = StepChangeUpdater()

        ts = [0.8, 0.2, 0.5]  # unsorted
        add_update!(sys, AtTimes(ts), upd)

        prob = sciml_odeproblem(sys, [0.0], (0.0, 1.0); p = nothing)
        sol = SciMLBase.solve(prob, OrdinaryDiffEq.Tsit5();
            reltol = 1e-10, abstol = 1e-10,
            save_everystep = true,
        )

        @test length(upd.fired) == 3
        # compare as sets with tolerance (avoid exact float equality)
        fired_sorted = sort(upd.fired)
        ts_sorted = sort(ts)
        @test all(i -> isapprox(fired_sorted[i], ts_sorted[i]; atol = 1000eps(ts_sorted[i])), 1:3)

        @test all(ti -> any(t -> isapprox(t, ti; atol = 1000eps(ti)), sol.t), ts)
    end

    @testset "SciML: AtTimes outside tspan ignored" begin
        sys = SciToySystem(0.0)
        upd = StepChangeUpdater()
        add_update!(sys, AtTimes([-1.0, 2.0, 0.5]), upd)

        prob = sciml_odeproblem(sys, [0.0], (0.0, 1.0); p = nothing)
        sol = SciMLBase.solve(prob, OrdinaryDiffEq.Tsit5();
            reltol = 1e-10, abstol = 1e-10,
            save_everystep = true,
        )

        @test length(upd.fired) == 1
        @test isapprox(upd.fired[1], 0.5; atol = 1000eps(0.5))
    end

    @testset "SciML: Periodic schedule fires on cadence via tstops" begin
        sys = SciToySystem(0.0)
        per_upd = SciPulseUpdater()
        add_update!(sys, Periodic(0.3, 0.0), per_upd)  # 0.0, 0.3, 0.6, 0.9

        prob = sciml_odeproblem(sys, [0.0], (0.0, 1.0); p = nothing)
        sol = SciMLBase.solve(prob, OrdinaryDiffEq.Tsit5();
            reltol = 1e-10, abstol = 1e-10,
            save_everystep = true,
        )

        expected = [0.0, 0.3, 0.6, 0.9]
        fired_sorted = sort(per_upd.fired)
        tol = 1e-10
        @test length(per_upd.fired) == length(expected)
        @test all(i -> isapprox(fired_sorted[i], expected[i]; atol = tol, rtol = 0.0), eachindex(expected))

        # solver visited these times
        @test all(ti -> any(t -> isapprox(t, ti; atol = tol, rtol = 0.0), sol.t), expected)
    end

    @testset "SciML: core callback + user callback both execute" begin
        sys = SciToySystem(0.0)
        upd = StepChangeUpdater()
        add_update!(sys, AtTimes([0.5]), upd)

        user_hits = Ref(0)
        user_cb = SciMLBase.DiscreteCallback(
            (u,t,integrator) -> isapprox(t, 0.5; atol = 1000eps(0.5)),
            integrator -> (user_hits[] += 1)
        )

        prob = sciml_odeproblem(sys, [0.0], (0.0, 1.0); callback = user_cb, p = nothing)
        sol = SciMLBase.solve(prob, OrdinaryDiffEq.Tsit5(); reltol=1e-10, abstol=1e-10)

        @test length(upd.fired) == 1
        @test user_hits[] == 1
    end

    @testset "SciML: EveryStep not included unless requested" begin
        sys = SciToySystem(0.0)

        per = SciPulseUpdater()
        add_update!(sys, EveryStep(), per)

        prob = sciml_odeproblem(sys, [0.0], (0.0, 0.2); p=nothing)  # default include_every_step=false
        sol = SciMLBase.solve(prob, OrdinaryDiffEq.Tsit5(); reltol=1e-10, abstol=1e-10)

        @test isempty(per.fired)

        prob2 = sciml_odeproblem(sys, [0.0], (0.0, 0.2); p=nothing, include_every_step=true)
        sol2 = SciMLBase.solve(prob2, OrdinaryDiffEq.Tsit5(); reltol=1e-10, abstol=1e-10)

        @test !isempty(per.fired)
    end
end
