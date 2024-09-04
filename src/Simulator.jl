mutable struct Simulator{R<:AbstractReactionSystem}
    RS::R
    initial_reactant_counts::Vector{Float64}
    reactant_counts::Vector{Float64}
    reaction_history::Vector{Tuple{Float64,Int}}
    reaction_counts::Vector{Int}
    Q::PriorityQueue{Int,Float64}
    t::Float64
    function Simulator(
        reactants::Dict{V,ReactantId},
        reactions::Vector{AbstractReaction},
        initial_state::Vector{Float64}
    ) where {V<:Any}
        RS = ReactionSystem(reactants, reactions)
        N_reactions = RS.N_reactions
        new{ReactionSystem{Symbol}}(
            RS, 
            copy(initial_state), copy(initial_state), 
            Tuple{Float64,Int}[], zeros(Int, N_reactions),
            PriorityQueue{Int,Float64}(), 0.0)
    end
end


function Base.show(io::IO, A::Simulator)
    print(io,
        ("Reaction system (t = $(round(A.t,digits=4))):\n    Reactants:\n        "
         * join(["$k ($(A.reactant_counts[v]))" for (k, v) in A.RS.reactants], "\n        ")
         * "\n    Reactions:\n        "
         * join([to_string(x, A.RS) for x in A.RS.reactions], "\n        ")))
end


function init!(A::Simulator{S}, initial_state::Vector{Float64}) where {S<:AbstractReactionSystem}
    A.t = 0
    copyto!(A.initial_reactant_counts, initial_state)
    copyto!(A.reactant_counts, initial_state)
    A.reaction_counts = zeros(Int, A.N_reactions)
    A.reaction_history = Tuple{Float64,Int}[]
end


function reconstruct_reactant_counts_from_history(A::Simulator{S}) where {S<:AbstractReactionSystem}
    current_state = copy(A.initial_reactant_counts)
    current_time = 0.0
    counts_history = Tuple{Float64,Vector{Float64}}[(current_time, current_state),]
    for (tau, rid) in A.reaction_history
        current_time -= tau
        current_state .-= A.RS.stoichiometries[:, rid]
        push!(counts_history, (current_time, current_state))
    end
    return counts_history
end
