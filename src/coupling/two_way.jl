function _default_fields_for_block(problem::CoupledProblem, block_name::Symbol)
    fields = Symbol[]
    for map in problem.maps
        if map.from === block_name || map.to === block_name
            push!(fields, map.field)
        end
    end
    return unique(fields)
end

function _residual_fields(problem::CoupledProblem, mode::TwoWayCoupling, block_name::Symbol)
    if !isempty(mode.norm_fields)
        return mode.norm_fields
    end
    return _default_fields_for_block(problem, block_name)
end

function _snapshot_block(block::CoupledBlock)
    return CoupledBlock(block.name, block.model, copy_state(block.state);
                        cache=block.cache,
                        metadata=block.metadata)
end

function _snapshot_blocks(problem::CoupledProblem, order::AbstractVector{Symbol})
    snapshots = Dict{Symbol,CoupledBlock}()
    for name in order
        block = _block_or_error(problem, name)
        snapshots[name] = _snapshot_block(block)
    end
    return snapshots
end

function _relax_block_fields!(problem::CoupledProblem,
                              mode::TwoWayCoupling,
                              old_block::CoupledBlock,
                              new_block::CoupledBlock)
    ω = get(mode.relaxation, new_block.name, one(mode.atol))
    ω == one(mode.atol) && return new_block

    mode.verbose && _log_relaxation(new_block.name, ω)
    fields = _residual_fields(problem, mode, new_block.name)
    for field in fields
        old_field = _get_coupling_field_or_error(old_block, field)
        new_field = _get_coupling_field_or_error(new_block, field)
        relaxed = underrelax!(new_field, old_field, ω)
        _set_coupling_field_or_error!(new_block, field, relaxed)
    end
    return new_block
end

function _two_way_threshold(mode::TwoWayCoupling, baseline::Float64)
    return Float64(mode.atol + mode.rtol * baseline)
end

function _run_two_way!(problem::CoupledProblem,
                       mode::TwoWayCoupling,
                       stepper;
                       return_history::Bool,
                       store_history::Bool,
                       throw_on_nonconvergence::Bool)
    _validate_problem_common!(problem)
    _validate_order!(problem, mode.order; mode_name="TwoWayCoupling")

    for name in keys(mode.relaxation)
        _block_index(problem, name) === nothing && throw(ArgumentError(
            "TwoWayCoupling relaxation entry `$name` does not match any block."))
    end

    history = CouplingHistory()
    apply_coupling_maps!(problem)

    baseline = 0.0
    converged = false
    for iter in 1:mode.maxiter
        snapshots = _snapshot_blocks(problem, mode.order)

        if mode.sweep === :GaussSeidel
            for name in mode.order
                apply_coupling_maps!(problem; to=name)
                block = _block_or_error(problem, name)
                stepper(block)
                _relax_block_fields!(problem, mode, snapshots[name], block)
                apply_coupling_maps!(problem; from=name)
            end
        else
            apply_coupling_maps!(problem)
            for name in mode.order
                block = _block_or_error(problem, name)
                stepper(block)
            end
            for name in mode.order
                block = _block_or_error(problem, name)
                _relax_block_fields!(problem, mode, snapshots[name], block)
            end
            apply_coupling_maps!(problem)
        end

        block_residuals = Dict{Symbol,Float64}()
        for name in mode.order
            old_block = snapshots[name]
            new_block = _block_or_error(problem, name)
            fields = _residual_fields(problem, mode, name)
            block_residuals[name] = compute_block_residual(old_block, new_block, fields)
        end

        global_residual = compute_global_residual(block_residuals)
        if iter == 1
            baseline = max(global_residual, eps(Float64))
        end

        threshold = _two_way_threshold(mode, baseline)
        _record_history!(history, iter, global_residual, block_residuals)

        mode.verbose && _log_two_way_iteration(iter, block_residuals, global_residual, threshold)

        if global_residual <= threshold
            converged = true
            break
        end
    end

    if !converged && throw_on_nonconvergence
        last_residual = isempty(history.residuals) ? Inf : history.residuals[end]
        throw(ErrorException(
            "TwoWayCoupling did not converge after $(mode.maxiter) iterations; " *
            "last residual=$last_residual, atol=$(mode.atol), rtol=$(mode.rtol)."))
    end

    return _finalize_coupled_result(problem, history;
                                    return_history=return_history,
                                    store_history=store_history)
end

function _solve_two_way!(problem::CoupledProblem,
                         mode::TwoWayCoupling;
                         return_history::Bool=false,
                         store_history::Bool=false,
                         throw_on_nonconvergence::Bool=true,
                         kwargs...)
    stepper = block -> advance_steady!(block; kwargs...)
    return _run_two_way!(problem, mode, stepper;
                         return_history=return_history,
                         store_history=store_history,
                         throw_on_nonconvergence=throw_on_nonconvergence)
end

function _step_two_way!(problem::CoupledProblem,
                        mode::TwoWayCoupling,
                        t,
                        dt;
                        return_history::Bool=false,
                        store_history::Bool=false,
                        throw_on_nonconvergence::Bool=true,
                        kwargs...)
    stepper = block -> advance_unsteady!(block, t, dt; kwargs...)
    return _run_two_way!(problem, mode, stepper;
                         return_history=return_history,
                         store_history=store_history,
                         throw_on_nonconvergence=throw_on_nonconvergence)
end
