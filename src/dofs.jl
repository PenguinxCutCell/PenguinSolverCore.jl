struct DofMap{I<:Integer}
    indices::Vector{I}
end

DofMap(indices::AbstractVector{I}) where {I<:Integer} = DofMap{I}(collect(indices))

restrict(u, map::DofMap) = u[map.indices]

function prolong!(u, reduced, map::DofMap)
    length(reduced) == length(map.indices) || throw(DimensionMismatch("reduced vector length does not match dof map size"))
    for (k, idx) in enumerate(map.indices)
        u[idx] = reduced[k]
    end
    return u
end
