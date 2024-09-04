function propensity(reactant_counts::Vector{T}, r::Reaction{T})::Float64 where {T<:Real}
    log_rate = 0.0
    for (a, c) in r.reactants
        log_rate += log(reactant_counts[a]) * c
    end
    return exp(log_rate) * r.rate_constant
end

function propensity(reactant_counts::Vector{T}, r::SimpleReaction)::Float64 where {T<:Real}
    rate = r.rate_constant
    for a in r.reactants
        rate *= reactant_counts[a]
    end
    return rate
end

function calculate_all_propensity(rs::ReactionSystem{S})::Vector{Float64} where {S<:Any}
    return [propensity(rs.reactant_counts, r) for r in rs.reactions]
end

function calculate_all_propensity!(
        p::Vector{Float64}, 
        reactions::Vector{AbstractReaction}, 
        reactant_counts::Vector{T}
) where {T<:Real}
    s = 0.0
    for (i,r) in enumerate(reactions)
        p[i] = propensity(reactant_counts, r)
        s += p[i]
    end
    return s
end