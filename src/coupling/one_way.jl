"""
    solve_coupled!(problem::CoupledProblem; return_history=false, store_history=false, kwargs...)

Run coupled steady orchestration.

- `OneWayCoupling`: one ordered pass.
- `TwoWayCoupling`: outer Picard iteration.
"""
function solve_coupled!(problem::CoupledProblem;
                        return_history::Bool=false,
                        store_history::Bool=false,
                        kwargs...)
    mode = problem.coupling_mode
    if mode isa OneWayCoupling
        return _solve_one_way!(problem, mode; return_history=return_history, store_history=store_history, kwargs...)
    elseif mode isa TwoWayCoupling
        return _solve_two_way!(problem, mode; return_history=return_history, store_history=store_history, kwargs...)
    end

    throw(ArgumentError("Unsupported coupling mode type: $(typeof(mode))."))
end

"""
    step_coupled!(problem::CoupledProblem, t, dt; return_history=false, store_history=false, kwargs...)

Run one coupled unsteady step.

- `OneWayCoupling`: one ordered pass.
- `TwoWayCoupling`: outer Picard iteration at `(t, dt)`.
"""
function step_coupled!(problem::CoupledProblem, t, dt;
                       return_history::Bool=false,
                       store_history::Bool=false,
                       kwargs...)
    mode = problem.coupling_mode
    if mode isa OneWayCoupling
        return _step_one_way!(problem, mode, t, dt;
                              return_history=return_history,
                              store_history=store_history,
                              kwargs...)
    elseif mode isa TwoWayCoupling
        return _step_two_way!(problem, mode, t, dt;
                              return_history=return_history,
                              store_history=store_history,
                              kwargs...)
    end

    throw(ArgumentError("Unsupported coupling mode type: $(typeof(mode))."))
end

function _solve_one_way!(problem::CoupledProblem,
                         mode::OneWayCoupling;
                         return_history::Bool=false,
                         store_history::Bool=false,
                         kwargs...)
    _validate_problem_common!(problem)
    _validate_order!(problem, mode.order; mode_name="OneWayCoupling")

    history = CouplingHistory()
    for name in mode.order
        block = _block_or_error(problem, name)
        advance_steady!(block; kwargs...)
        apply_coupling_maps!(problem; from=name)
    end
    _record_history!(history, 1, 0.0, Dict{Symbol,Float64}(name => 0.0 for name in mode.order))

    return _finalize_coupled_result(problem, history;
                                    return_history=return_history,
                                    store_history=store_history)
end

function _step_one_way!(problem::CoupledProblem,
                        mode::OneWayCoupling,
                        t,
                        dt;
                        return_history::Bool=false,
                        store_history::Bool=false,
                        kwargs...)
    _validate_problem_common!(problem)
    _validate_order!(problem, mode.order; mode_name="OneWayCoupling")

    history = CouplingHistory()
    for name in mode.order
        block = _block_or_error(problem, name)
        advance_unsteady!(block, t, dt; kwargs...)
        apply_coupling_maps!(problem; from=name)
    end
    _record_history!(history, 1, 0.0, Dict{Symbol,Float64}(name => 0.0 for name in mode.order))

    return _finalize_coupled_result(problem, history;
                                    return_history=return_history,
                                    store_history=store_history)
end
