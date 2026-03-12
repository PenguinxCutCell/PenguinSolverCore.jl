block_names(problem::CoupledProblem) = [block.name for block in problem.blocks]

function _block_index(problem::CoupledProblem, name::Symbol)
    for (idx, block) in pairs(problem.blocks)
        if block.name === name
            return idx
        end
    end
    return nothing
end

function _block_or_error(problem::CoupledProblem, name::Symbol)
    idx = _block_index(problem, name)
    idx === nothing && throw(ArgumentError(
        "Unknown block `$name`. Available blocks: $(block_names(problem))."))
    return problem.blocks[idx]
end

function _validate_problem_common!(problem::CoupledProblem)
    isempty(problem.blocks) && throw(ArgumentError("CoupledProblem must contain at least one block."))

    names = block_names(problem)
    length(names) == length(unique(names)) ||
        throw(ArgumentError("CoupledProblem block names must be unique, got $names."))

    for map in problem.maps
        _block_index(problem, map.from) === nothing && throw(ArgumentError(
            "CouplingMap source block `$(map.from)` not found. Available blocks: $names."))
        _block_index(problem, map.to) === nothing && throw(ArgumentError(
            "CouplingMap target block `$(map.to)` not found. Available blocks: $names."))
    end

    return problem
end

function _validate_order!(problem::CoupledProblem, order::AbstractVector{Symbol}; mode_name::AbstractString)
    isempty(order) && throw(ArgumentError("$mode_name order cannot be empty."))
    _check_unique_symbols(order, mode_name)

    names = block_names(problem)
    order_set = Set(order)
    name_set = Set(names)
    order_set == name_set || throw(ArgumentError(
        "$mode_name order must include each block exactly once. order=$(collect(order)), blocks=$names"))

    return order
end

function _get_coupling_field_or_error(block::CoupledBlock, field::Symbol)
    try
        return get_coupling_field(block, Val(field))
    catch err
        if err isa MethodError && err.f === get_coupling_field
            throw(ArgumentError(
                "Missing get_coupling_field implementation for block `$(block.name)` and field `:$field`."))
        end
        rethrow()
    end
end

function _set_coupling_field_or_error!(block::CoupledBlock, field::Symbol, data)
    try
        return set_coupling_field!(block, Val(field), data)
    catch err
        if err isa MethodError && err.f === set_coupling_field!
            throw(ArgumentError(
                "Missing set_coupling_field! implementation for block `$(block.name)` and field `:$field`."))
        end
        rethrow()
    end
end

function _apply_transfer(map::CouplingMap, data, from_block::CoupledBlock, to_block::CoupledBlock, problem::CoupledProblem)
    transfer = map.apply!
    if applicable(transfer, data, from_block, to_block, problem)
        return transfer(data, from_block, to_block, problem)
    elseif applicable(transfer, data, from_block, to_block)
        return transfer(data, from_block, to_block)
    elseif applicable(transfer, data)
        return transfer(data)
    end

    throw(ArgumentError(
        "CouplingMap apply! for $(map.from)->$(map.to):$(map.field) must accept one of " *
        "(data), (data, from_block, to_block), or (data, from_block, to_block, problem)."))
end

"""
    apply_coupling_map!(problem, map)

Apply one directed coupling transfer map.
"""
function apply_coupling_map!(problem::CoupledProblem, map::CouplingMap)
    from_block = _block_or_error(problem, map.from)
    to_block = _block_or_error(problem, map.to)

    data = _get_coupling_field_or_error(from_block, map.field)
    mapped_data = _apply_transfer(map, data, from_block, to_block, problem)
    _set_coupling_field_or_error!(to_block, map.field, mapped_data)

    return mapped_data
end

"""
    apply_coupling_maps!(problem; from=nothing, to=nothing)

Apply all coupling maps in `problem`, optionally filtered by source/target block name.
"""
function apply_coupling_maps!(problem::CoupledProblem;
                              from::Union{Nothing,Symbol}=nothing,
                              to::Union{Nothing,Symbol}=nothing)
    for map in problem.maps
        (from === nothing || map.from === from) || continue
        (to === nothing || map.to === to) || continue
        apply_coupling_map!(problem, map)
    end
    return problem
end
