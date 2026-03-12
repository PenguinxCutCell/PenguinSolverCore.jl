using Test
using SparseArrays
using PenguinSolverCore

# -----------------------------------------------------------------------------
# Existing LinearSystem tests
# -----------------------------------------------------------------------------

struct ToyModel end

function PenguinSolverCore.assemble!(sys::LinearSystem, ::ToyModel, t, dt)
    n = length(sys.b)
    sys.A = spdiagm(0 => fill(2.0, n))
    sys.b .= 1.0
    return sys
end

function PenguinSolverCore.assemble_matrix!(sys::LinearSystem, ::ToyModel, t, dt, θ)
    n = length(sys.b)
    sys.A = spdiagm(0 => fill(2.0 + dt + θ, n))
    return sys
end

function PenguinSolverCore.assemble_rhs!(sys::LinearSystem, ::ToyModel, uⁿ, t, dt, θ)
    sys.b .= uⁿ .+ t .+ dt .+ θ
    return sys
end

@testset "LinearSystem solve!" begin
    A = spdiagm(0 => [2.0, 4.0])
    b = [2.0, 8.0]
    sys = LinearSystem(A, b)
    x = solve!(sys; method=:direct)
    @test x ≈ [1.0, 2.0]
end

@testset "theta_step! skeleton" begin
    sys = LinearSystem(spdiagm(0 => ones(2)), zeros(2))
    model = ToyModel()
    u0 = [1.0, 3.0]
    x = theta_step!(sys, model, u0, 0.0, 0.1, 0.5; method=:direct)
    @test length(x) == 2
    @test sys.last_dt == 0.1
end

# -----------------------------------------------------------------------------
# Coupling toy models
# -----------------------------------------------------------------------------

struct ProducerModel
    value::Float64
    log::Vector{Symbol}
end

mutable struct ProducerState
    velocity::Float64
    steady_calls::Int
    unsteady_calls::Int
end

struct FollowerModel
    gain::Float64
    log::Vector{Symbol}
end

mutable struct FollowerState
    velocity::Float64
    pressure::Float64
    steady_calls::Int
    unsteady_calls::Int
end

struct VelocityFromTemperatureModel
    α::Float64
    c::Float64
    log::Vector{Symbol}
end

struct TemperatureFromVelocityModel
    β::Float64
    c::Float64
    log::Vector{Symbol}
end

mutable struct ReciprocalState
    velocity::Float64
    temperature::Float64
    steady_calls::Int
    unsteady_calls::Int
end

struct NoGetModel end
struct NoSetModel end

mutable struct NoGetState
    foo::Float64
end

mutable struct NoSetState
    foo::Float64
end

# Snapshot hooks used by outer-coupling residual computations
PenguinSolverCore.copy_state(state::ProducerState) = ProducerState(state.velocity, state.steady_calls, state.unsteady_calls)
PenguinSolverCore.copy_state(state::FollowerState) = FollowerState(state.velocity, state.pressure, state.steady_calls, state.unsteady_calls)
PenguinSolverCore.copy_state(state::ReciprocalState) = ReciprocalState(state.velocity, state.temperature, state.steady_calls, state.unsteady_calls)

# Initialization hooks
PenguinSolverCore.initialize_state(::ProducerModel, init) = ProducerState(Float64(init === nothing ? 0.0 : init), 0, 0)
PenguinSolverCore.initialize_state(::FollowerModel, init) = FollowerState(Float64(init === nothing ? 0.0 : init), 0.0, 0, 0)
PenguinSolverCore.initialize_state(model::VelocityFromTemperatureModel, init) =
    ReciprocalState(0.0, Float64(init === nothing ? 0.0 : init), 0, 0)
PenguinSolverCore.initialize_state(model::TemperatureFromVelocityModel, init) =
    ReciprocalState(Float64(init === nothing ? 0.0 : init), 0.0, 0, 0)
PenguinSolverCore.initialize_state(::NoGetModel, init) = NoGetState(Float64(init === nothing ? 0.0 : init))
PenguinSolverCore.initialize_state(::NoSetModel, init) = NoSetState(Float64(init === nothing ? 0.0 : init))

# Solve hooks
function PenguinSolverCore.advance_steady!(block::CoupledBlock{M,S,C}; kwargs...) where {M<:ProducerModel,S<:ProducerState,C}
    block.state.velocity = block.model.value
    block.state.steady_calls += 1
    push!(block.model.log, block.name)
    return block
end

function PenguinSolverCore.advance_unsteady!(block::CoupledBlock{M,S,C}, t, dt; kwargs...) where {M<:ProducerModel,S<:ProducerState,C}
    block.state.velocity = block.model.value + t + dt
    block.state.unsteady_calls += 1
    push!(block.model.log, block.name)
    return block
end

function PenguinSolverCore.advance_steady!(block::CoupledBlock{M,S,C}; kwargs...) where {M<:FollowerModel,S<:FollowerState,C}
    block.state.pressure = block.model.gain * block.state.velocity
    block.state.steady_calls += 1
    push!(block.model.log, block.name)
    return block
end

function PenguinSolverCore.advance_unsteady!(block::CoupledBlock{M,S,C}, t, dt; kwargs...) where {M<:FollowerModel,S<:FollowerState,C}
    block.state.pressure = block.model.gain * block.state.velocity + dt
    block.state.unsteady_calls += 1
    push!(block.model.log, block.name)
    return block
end

function PenguinSolverCore.advance_steady!(block::CoupledBlock{M,S,C}; kwargs...) where {M<:VelocityFromTemperatureModel,S<:ReciprocalState,C}
    block.state.velocity = block.model.α * block.state.temperature + block.model.c
    block.state.steady_calls += 1
    push!(block.model.log, block.name)
    return block
end

function PenguinSolverCore.advance_unsteady!(block::CoupledBlock{M,S,C}, t, dt; kwargs...) where {M<:VelocityFromTemperatureModel,S<:ReciprocalState,C}
    PenguinSolverCore.advance_steady!(block; kwargs...)
    block.state.unsteady_calls += 1
    return block
end

function PenguinSolverCore.advance_steady!(block::CoupledBlock{M,S,C}; kwargs...) where {M<:TemperatureFromVelocityModel,S<:ReciprocalState,C}
    block.state.temperature = block.model.β * block.state.velocity + block.model.c
    block.state.steady_calls += 1
    push!(block.model.log, block.name)
    return block
end

function PenguinSolverCore.advance_unsteady!(block::CoupledBlock{M,S,C}, t, dt; kwargs...) where {M<:TemperatureFromVelocityModel,S<:ReciprocalState,C}
    PenguinSolverCore.advance_steady!(block; kwargs...)
    block.state.unsteady_calls += 1
    return block
end

PenguinSolverCore.advance_steady!(block::CoupledBlock{M,S,C}; kwargs...) where {M<:NoGetModel,S<:NoGetState,C} = block
PenguinSolverCore.advance_steady!(block::CoupledBlock{M,S,C}; kwargs...) where {M<:NoSetModel,S<:NoSetState,C} = block
PenguinSolverCore.advance_unsteady!(block::CoupledBlock{M,S,C}, t, dt; kwargs...) where {M<:NoGetModel,S<:NoGetState,C} = block
PenguinSolverCore.advance_unsteady!(block::CoupledBlock{M,S,C}, t, dt; kwargs...) where {M<:NoSetModel,S<:NoSetState,C} = block

# Coupling field hooks
function PenguinSolverCore.get_coupling_field(block::CoupledBlock{M,S,C}, ::Val{:velocity}) where {M<:ProducerModel,S<:ProducerState,C}
    return block.state.velocity
end

function PenguinSolverCore.set_coupling_field!(block::CoupledBlock{M,S,C}, ::Val{:velocity}, data) where {M<:ProducerModel,S<:ProducerState,C}
    block.state.velocity = Float64(data)
    return block
end

function PenguinSolverCore.get_coupling_field(block::CoupledBlock{M,S,C}, ::Val{:velocity}) where {M<:FollowerModel,S<:FollowerState,C}
    return block.state.velocity
end

function PenguinSolverCore.set_coupling_field!(block::CoupledBlock{M,S,C}, ::Val{:velocity}, data) where {M<:FollowerModel,S<:FollowerState,C}
    block.state.velocity = Float64(data)
    return block
end

function PenguinSolverCore.get_coupling_field(block::CoupledBlock{M,S,C}, ::Val{:pressure}) where {M<:FollowerModel,S<:FollowerState,C}
    return block.state.pressure
end

function PenguinSolverCore.get_coupling_field(block::CoupledBlock{M,S,C}, ::Val{:velocity}) where {M,S<:ReciprocalState,C}
    return block.state.velocity
end

function PenguinSolverCore.get_coupling_field(block::CoupledBlock{M,S,C}, ::Val{:temperature}) where {M,S<:ReciprocalState,C}
    return block.state.temperature
end

function PenguinSolverCore.set_coupling_field!(block::CoupledBlock{M,S,C}, ::Val{:velocity}, data) where {M,S<:ReciprocalState,C}
    block.state.velocity = Float64(data)
    return block
end

function PenguinSolverCore.set_coupling_field!(block::CoupledBlock{M,S,C}, ::Val{:temperature}, data) where {M,S<:ReciprocalState,C}
    block.state.temperature = Float64(data)
    return block
end

function PenguinSolverCore.set_coupling_field!(block::CoupledBlock{M,S,C}, ::Val{:foo}, data) where {M<:NoGetModel,S<:NoGetState,C}
    block.state.foo = Float64(data)
    return block
end

function PenguinSolverCore.get_coupling_field(block::CoupledBlock{M,S,C}, ::Val{:foo}) where {M<:NoSetModel,S<:NoSetState,C}
    return block.state.foo
end

# -----------------------------------------------------------------------------
# Coupling helpers
# -----------------------------------------------------------------------------

function make_reciprocal_problem(α, β, c1, c2;
                                 ω=1.0,
                                 maxiter=80,
                                 atol=1e-10,
                                 rtol=1e-8,
                                 sweep=:GaussSeidel)
    log = Symbol[]
    block_a = CoupledBlock(:A,
                           VelocityFromTemperatureModel(α, c1, log),
                           ReciprocalState(0.0, 0.0, 0, 0))
    block_b = CoupledBlock(:B,
                           TemperatureFromVelocityModel(β, c2, log),
                           ReciprocalState(0.0, 0.0, 0, 0))

    relax = ω == 1.0 ? Dict{Symbol,Float64}() : Dict(:A => ω, :B => ω)
    mode = TwoWayCoupling([:A, :B];
                          maxiter=maxiter,
                          atol=atol,
                          rtol=rtol,
                          relaxation=relax,
                          norm_fields=[:velocity, :temperature],
                          sweep=sweep,
                          verbose=false)

    maps = [
        CouplingMap(:A, :B, :velocity),
        CouplingMap(:B, :A, :temperature),
    ]

    return CoupledProblem([block_a, block_b], mode; maps=maps), log
end

# -----------------------------------------------------------------------------
# Coupling tests
# -----------------------------------------------------------------------------

@testset "OneWayCoupling sequential transfer" begin
    log = Symbol[]
    block_a = CoupledBlock(:A, ProducerModel(2.5, log), ProducerState(0.0, 0, 0))
    block_b = CoupledBlock(:B, FollowerModel(3.0, log), FollowerState(0.0, 0.0, 0, 0))

    problem = CoupledProblem(
        [block_a, block_b],
        OneWayCoupling([:A, :B]);
        maps=[CouplingMap(:A, :B, :velocity)],
    )

    result = solve_coupled!(problem)

    @test result === problem
    @test log == [:A, :B]
    @test block_a.state.steady_calls == 1
    @test block_b.state.steady_calls == 1
    @test block_b.state.velocity == 2.5
    @test block_b.state.pressure ≈ 7.5

    step_coupled!(problem, 1.0, 0.25)
    @test block_a.state.unsteady_calls == 1
    @test block_b.state.unsteady_calls == 1
    @test block_b.state.velocity == 3.75
    @test block_b.state.pressure ≈ (3.0 * 3.75 + 0.25)
end

@testset "TwoWayCoupling fixed-point convergence" begin
    α = 0.6
    β = 0.4
    c1 = 1.2
    c2 = -0.3

    problem, log = make_reciprocal_problem(α, β, c1, c2; maxiter=80, atol=1e-12, rtol=1e-10)
    solved, history = solve_coupled!(problem; return_history=true)

    @test solved === problem
    @test !isempty(history.iterations)
    @test history.residuals[end] < history.residuals[1]

    threshold = solved.coupling_mode.atol + solved.coupling_mode.rtol * max(history.residuals[1], eps())
    @test history.residuals[end] <= threshold

    texact = (β * c1 + c2) / (1 - α * β)
    vexact = α * texact + c1

    block_a = solved.blocks[1]
    block_b = solved.blocks[2]

    @test block_a.state.velocity ≈ vexact atol=1e-8
    @test block_b.state.temperature ≈ texact atol=1e-8
    @test block_a.state.temperature ≈ texact atol=1e-8
    @test block_b.state.velocity ≈ vexact atol=1e-8
    @test log[1:2] == [:A, :B]
end

@testset "TwoWayCoupling unsteady step" begin
    problem, _ = make_reciprocal_problem(0.6, 0.4, 1.2, -0.3; maxiter=40, atol=1e-10, rtol=1e-8)
    stepped, history = step_coupled!(problem, 0.0, 0.1; return_history=true)
    @test stepped === problem
    @test !isempty(history.iterations)
    @test stepped.blocks[1].state.unsteady_calls > 0
    @test stepped.blocks[2].state.unsteady_calls > 0
end

@testset "TwoWayCoupling under-relaxation impact" begin
    α = 2.0
    β = -0.49
    c1 = 0.0
    c2 = 1.0

    plain_problem, _ = make_reciprocal_problem(α, β, c1, c2; ω=1.0, maxiter=40, atol=1e-10, rtol=1e-10)
    _, plain_history = solve_coupled!(plain_problem; return_history=true, throw_on_nonconvergence=false)

    relaxed_problem, _ = make_reciprocal_problem(α, β, c1, c2; ω=0.5, maxiter=40, atol=1e-10, rtol=1e-10)
    _, relaxed_history = solve_coupled!(relaxed_problem; return_history=true)

    @test length(plain_history.iterations) == 40
    @test length(relaxed_history.iterations) < length(plain_history.iterations)
    @test relaxed_history.residuals[end] < plain_history.residuals[end]

    texact = (β * c1 + c2) / (1 - α * β)
    vexact = α * texact + c1

    @test relaxed_problem.blocks[1].state.velocity ≈ vexact atol=1e-6
    @test relaxed_problem.blocks[2].state.temperature ≈ texact atol=1e-6
end

@testset "Missing coupling hook errors" begin
    missing_get_problem = CoupledProblem(
        [
            CoupledBlock(:A, NoGetModel(), NoGetState(1.0)),
            CoupledBlock(:B, NoSetModel(), NoSetState(0.0)),
        ],
        OneWayCoupling([:A, :B]);
        maps=[CouplingMap(:A, :B, :foo)],
    )

    err_get = try
        PenguinSolverCore.apply_coupling_map!(missing_get_problem, missing_get_problem.maps[1])
        nothing
    catch err
        err
    end
    @test err_get isa ArgumentError
    @test occursin("Missing get_coupling_field", sprint(showerror, err_get))

    missing_set_problem = CoupledProblem(
        [
            CoupledBlock(:A, NoSetModel(), NoSetState(1.0)),
            CoupledBlock(:B, NoSetModel(), NoSetState(0.0)),
        ],
        OneWayCoupling([:A, :B]);
        maps=[CouplingMap(:A, :B, :foo)],
    )

    err_set = try
        PenguinSolverCore.apply_coupling_map!(missing_set_problem, missing_set_problem.maps[1])
        nothing
    catch err
        err
    end
    @test err_set isa ArgumentError
    @test occursin("Missing set_coupling_field!", sprint(showerror, err_set))
end

@testset "CouplingHistory recording" begin
    problem, _ = make_reciprocal_problem(0.5, 0.25, 0.8, -0.1; maxiter=30, atol=1e-12, rtol=1e-10)
    _, history = solve_coupled!(problem; return_history=true)

    n = length(history.iterations)
    @test n == length(history.residuals)
    @test n == length(history.block_residuals)
    @test history.iterations == collect(1:n)
    @test all(all(name in (:A, :B) for name in keys(entry)) for entry in history.block_residuals)
end
