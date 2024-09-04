function update!(A::Simulator{S}, reaction_id::ReactantId, dt::Float64) where {S<:AbstractReactionSystem}
    @views A.reactant_counts .+= A.RS.stoichiometries[:, reaction_id]
    if any(x->x<0, A.reactant_counts)
        return false
    end
    A.reaction_counts[reaction_id] += 1
    A.t += dt
    return true
end
