"""
    CouplingHistory

Residual history for coupled outer iterations.
"""
mutable struct CouplingHistory
    iterations::Vector{Int}
    residuals::Vector{Float64}
    block_residuals::Vector{Dict{Symbol,Float64}}
end

CouplingHistory() = CouplingHistory(Int[], Float64[], Dict{Symbol,Float64}[])

function _record_history!(history::CouplingHistory,
                          iteration::Integer,
                          global_residual::Real,
                          block_residuals::Dict{Symbol,Float64})
    push!(history.iterations, Int(iteration))
    push!(history.residuals, Float64(global_residual))
    push!(history.block_residuals, Dict{Symbol,Float64}(pairs(block_residuals)))
    return history
end

function _store_history!(problem::CoupledProblem, history::CouplingHistory)
    for block in problem.blocks
        if block.cache isa AbstractDict
            block.cache[:coupling_history] = history
        elseif block.metadata isa AbstractDict
            block.metadata[:coupling_history] = history
        end
    end
    return history
end

function _finalize_coupled_result(problem::CoupledProblem,
                                  history::CouplingHistory;
                                  return_history::Bool,
                                  store_history::Bool)
    store_history && _store_history!(problem, history)
    return return_history ? (problem, history) : problem
end

function _block_label(block::CoupledBlock)
    try
        return block_summary(block)
    catch err
        if err isa MethodError && err.f === block_summary
            return String(block.name)
        end
        rethrow()
    end
end

function _log_one_way_step(step::Integer, block::CoupledBlock)
    println("[OneWayCoupling] step=$step block=$(_block_label(block))")
end

function _log_two_way_iteration(iter::Integer,
                                block_residuals::Dict{Symbol,Float64},
                                global_residual::Float64,
                                threshold::Float64)
    block_part = join(("$name=$(round(value, sigdigits=5))" for (name, value) in sort(collect(block_residuals), by=first)), ", ")
    println("[TwoWayCoupling] iter=$iter global=$(round(global_residual, sigdigits=6)) threshold=$(round(threshold, sigdigits=6)) blocks={$block_part}")
end

function _log_relaxation(block_name::Symbol, ω::Real)
    println("[TwoWayCoupling] relaxation block=$block_name ω=$ω")
end
