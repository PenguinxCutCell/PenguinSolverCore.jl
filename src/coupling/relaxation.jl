copy_state(state::Nothing) = nothing
copy_state(state::Number) = state
copy_state(state::AbstractArray) = copy(state)

function copy_state(state::Tuple)
    return map(copy_state, state)
end

function copy_state(state::NamedTuple)
    values = map(copy_state, Tuple(state))
    return NamedTuple{keys(state)}(values)
end

function copy_state(state::AbstractDict)
    copied = Dict{Any,Any}()
    for (k, v) in pairs(state)
        copied[k] = copy_state(v)
    end
    return copied
end

"""
    state_distance(a, b)

Infinity-norm-like distance between state-like objects.
"""
function state_distance(a::Number, b::Number)
    return abs(float(a - b))
end

function state_distance(a::AbstractArray, b::AbstractArray)
    size(a) == size(b) || throw(DimensionMismatch("state_distance array size mismatch: $(size(a)) vs $(size(b))."))

    d = 0.0
    @inbounds for i in eachindex(a, b)
        δ = abs(float(a[i] - b[i]))
        d = max(d, δ)
    end
    return d
end

function state_distance(a::Tuple, b::Tuple)
    length(a) == length(b) || throw(DimensionMismatch("state_distance tuple length mismatch."))

    d = 0.0
    for i in eachindex(a, b)
        d = max(d, state_distance(a[i], b[i]))
    end
    return d
end

function state_distance(a::NamedTuple, b::NamedTuple)
    keys(a) == keys(b) || throw(ArgumentError("state_distance NamedTuple keys mismatch."))

    d = 0.0
    for name in keys(a)
        d = max(d, state_distance(getfield(a, name), getfield(b, name)))
    end
    return d
end

function state_distance(a::AbstractDict, b::AbstractDict)
    keys(a) == keys(b) || throw(ArgumentError("state_distance dict keys mismatch."))

    d = 0.0
    for key in keys(a)
        d = max(d, state_distance(a[key], b[key]))
    end
    return d
end

"""
    underrelax!(new, old, ω)

Apply under-relaxation: `new <- ω*new + (1-ω)*old`.
"""
function underrelax!(new::Number, old::Number, ω::Real)
    return float(ω) * new + (1 - float(ω)) * old
end

function underrelax!(new::AbstractArray, old::AbstractArray, ω::Real)
    size(new) == size(old) || throw(DimensionMismatch("underrelax! array size mismatch."))

    ωf = float(ω)
    α = 1 - ωf
    @inbounds for i in eachindex(new, old)
        new[i] = ωf * new[i] + α * old[i]
    end
    return new
end

function underrelax!(new::Tuple, old::Tuple, ω::Real)
    length(new) == length(old) || throw(DimensionMismatch("underrelax! tuple length mismatch."))
    return ntuple(i -> underrelax!(new[i], old[i], ω), length(new))
end

function underrelax!(new::NamedTuple, old::NamedTuple, ω::Real)
    keys(new) == keys(old) || throw(ArgumentError("underrelax! NamedTuple keys mismatch."))
    values = ntuple(i -> underrelax!(Tuple(new)[i], Tuple(old)[i], ω), length(new))
    return NamedTuple{keys(new)}(values)
end

function underrelax!(new::AbstractDict, old::AbstractDict, ω::Real)
    keys(new) == keys(old) || throw(ArgumentError("underrelax! dict keys mismatch."))
    for key in keys(new)
        new[key] = underrelax!(new[key], old[key], ω)
    end
    return new
end

"""
    compute_block_residual(old_block, new_block, fields)

Compute one block residual for outer coupling convergence.
"""
function compute_block_residual(old_block::CoupledBlock,
                                new_block::CoupledBlock,
                                fields::AbstractVector{Symbol})
    try
        return Float64(coupling_residual(old_block, new_block, fields))
    catch err
        if err isa MethodError && err.f === coupling_residual
            return _default_block_residual(old_block, new_block, fields)
        end
        rethrow()
    end
end

function _default_block_residual(old_block::CoupledBlock,
                                 new_block::CoupledBlock,
                                 fields::AbstractVector{Symbol})
    if isempty(fields)
        return Float64(state_distance(old_block.state, new_block.state))
    end

    rmax = 0.0
    for field in fields
        old_data = _get_coupling_field_or_error(old_block, field)
        new_data = _get_coupling_field_or_error(new_block, field)
        rmax = max(rmax, state_distance(old_data, new_data))
    end
    return rmax
end

"""
    compute_global_residual(residuals)

Reduce block residuals to a single global residual.
"""
function compute_global_residual(residuals::Dict{Symbol,Float64})
    return isempty(residuals) ? 0.0 : maximum(values(residuals))
end

function compute_global_residual(residuals::Dict)
    if isempty(residuals)
        return 0.0
    end

    rmax = 0.0
    for value in values(residuals)
        rmax = max(rmax, Float64(value))
    end
    return rmax
end
