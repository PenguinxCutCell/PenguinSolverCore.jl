"""
    initialize_state(model, init)

Initialize a block state from a model and optional `init` payload.
Physics packages should overload this when constructing `CoupledBlock` without
an explicit `state` argument.
"""
function initialize_state(model, init)
    throw(MethodError(initialize_state, (model, init)))
end

"""
    advance_steady!(block::CoupledBlock; kwargs...)

Advance a coupled block to steady state.
"""
function advance_steady!(block::CoupledBlock; kwargs...)
    throw(MethodError(advance_steady!, (block,)))
end

"""
    advance_unsteady!(block::CoupledBlock, t, dt; kwargs...)

Advance a coupled block by one unsteady step at `(t, dt)`.
"""
function advance_unsteady!(block::CoupledBlock, t, dt; kwargs...)
    throw(MethodError(advance_unsteady!, (block, t, dt)))
end

"""
    get_coupling_field(block::CoupledBlock, ::Val{field}) where field

Extract a named coupling field from `block`.
"""
function get_coupling_field(block::CoupledBlock, ::Val{field}) where {field}
    throw(MethodError(get_coupling_field, (block, Val{field}())))
end

"""
    set_coupling_field!(block::CoupledBlock, ::Val{field}, data) where field

Inject a named coupling field into `block`.
"""
function set_coupling_field!(block::CoupledBlock, ::Val{field}, data) where {field}
    throw(MethodError(set_coupling_field!, (block, Val{field}(), data)))
end

"""
    copy_state(state)

Return a snapshot copy suitable for outer coupling residual comparisons.
"""
function copy_state(state)
    throw(MethodError(copy_state, (state,)))
end

"""
    coupling_residual(old_block::CoupledBlock, new_block::CoupledBlock, fields)

Optional custom residual hook for outer coupling convergence checks.
"""
function coupling_residual(old_block::CoupledBlock, new_block::CoupledBlock, fields)
    throw(MethodError(coupling_residual, (old_block, new_block, fields)))
end

"""
    block_summary(block::CoupledBlock)

Return a concise textual summary for diagnostics/logging.
"""
function block_summary(block::CoupledBlock)
    throw(MethodError(block_summary, (block,)))
end
