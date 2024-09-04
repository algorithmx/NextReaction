function simulate!(
    A::Simulator{S},
    t_end::Float64,
    sample_func::Function
) where {S<:AbstractReactionSystem}
    all_propensities = zeros(Float64, A.RS.N_reactions)
    succ = true
    while succ
        tot = calculate_all_propensity!(all_propensities, A.RS.reactions, A.reactant_counts)
        tot < 1e-16 && break
        tau, reaction_id = sample_func(all_propensities)
        push!(A.reaction_history, (A.t, reaction_id))
        succ = update!(A, reaction_id, tau)
        (A.t >= t_end) && break
    end
    return
end

simulate_direct_method!(rs::Simulator{S}, t_end::Float64) where {S<:AbstractReactionSystem} = simulate!(rs, t_end, sample_tau_and_reaction_1)

simulate_first_reaction!(rs::Simulator{S}, t_end::Float64) where {S<:AbstractReactionSystem} = simulate!(rs, t_end, sample_tau_and_reaction_2)