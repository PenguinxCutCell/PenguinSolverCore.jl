using Test
using SparseArrays
using PenguinSolverCore

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
