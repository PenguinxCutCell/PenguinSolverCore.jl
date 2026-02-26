mutable struct InvalidationCache
    dirty_matrix::Bool
    dirty_geometry::Bool
end

InvalidationCache() = InvalidationCache(false, false)

invalidate!(cache::InvalidationCache, ::Val{:matrix}) = (cache.dirty_matrix = true; cache)
invalidate!(cache::InvalidationCache, ::Val{:geometry}) = (cache.dirty_geometry = true; cache)

function invalidate!(cache::InvalidationCache, change::Symbol)
    if change === :matrix
        return invalidate!(cache, Val(:matrix))
    elseif change === :geometry
        return invalidate!(cache, Val(:geometry))
    elseif change === :rhs_only || change === :nothing
        return cache
    end
    throw(ArgumentError("unknown change flag `$change`; expected :rhs_only, :matrix, :geometry, or :nothing"))
end

needs_rebuild(cache::InvalidationCache) = cache.dirty_matrix || cache.dirty_geometry

function clear_invalidations!(cache::InvalidationCache)
    cache.dirty_matrix = false
    cache.dirty_geometry = false
    return cache
end

function maybe_invalidation_cache(sys)
    hasproperty(sys, :cache) || return nothing
    cache = getproperty(sys, :cache)
    cache isa InvalidationCache || throw(ArgumentError("system `$(typeof(sys))` has `cache` field but it is not an InvalidationCache"))
    return cache
end
