abstract type AbstractReactionSystem end

mutable struct ReactionSystem{S} <: AbstractReactionSystem
    N_reactants::Int
    N_reactions::Int
    reactant_list::Vector{S}
    reactants::Dict{S,ReactantId}
    reactions::Vector{AbstractReaction}
    stoichiometries::Matrix{Float64}
    masks::Vector{BitVector}
    DG::SimpleDiGraph
    function ReactionSystem(
        reactants::Dict{V,ReactantId},
        reactions::Vector{AbstractReaction}
    ) where {V<:Any}
        N_reactants = length(reactants)
        N_reactions = length(reactions)
        r_list = Vector{V}(undef, N_reactants)
        for (k, v) in reactants
            r_list[v] = k
        end
        new{V}(
            N_reactants, N_reactions,
            r_list, reactants, reactions, 
            build_stoichiometries(reactants, reactions)...,
            build_dependency_graph(reactions))
    end
end


to_string(r::SimpleReaction, s::ReactionSystem) =
    (join([string(s.reactant_list[r]) for r in r.reactants], " + ")
     * " → "
     * join([string(s.reactant_list[r]) for r in r.products], " + ")
     * " rate $(round(r.rate_constant,digits=3))")


to_string(r::Reaction{T}, s::ReactionSystem) where {T<:Integer} =
    (join(["$v $(s.reactant_list[k])" for (k, v) in r.reactants], " + ")
     * " → "
     * join(["$v $(s.reactant_list[k])" for (k, v) in r.products], " + ")
     * " rate $(round(r.rate_constant,digits=8))")


function build_stoichiometries(
    reactants::Dict{V,ReactantId},
    reactions::Vector{AbstractReaction}
)::Tuple{Matrix{Float64},Vector{BitVector}} where {V<:Any}
    N_reactants = length(reactants)
    N_reactions = length(reactions)
    stoichiometries = zeros(Float64, N_reactants, N_reactions)
    masks = Vector{BitVector}(undef, N_reactions)
    for (j, r) in enumerate(reactions)
        if r isa SimpleReaction
            for i in r.reactants
                stoichiometries[i, j] = -1
            end
            for i in r.products
                stoichiometries[i, j] = 1
            end
        else
            for (i, v) in r.reactants
                stoichiometries[i, j] = -v
            end
            for (i, v) in r.products
                stoichiometries[i, j] = v
            end
        end
    end
    for (j, r) in enumerate(reactions)
        masks[j] = BitVector(stoichiometries[:, j] .!= 0)
    end
    return stoichiometries, masks
end


function build_dependency_graph(reactions::Vector{R})::SimpleDiGraph where {R<:AbstractReaction}
    num_reactions = length(reactions)
    graph = SimpleDiGraph(num_reactions)
    for i in 1:num_reactions
        for j in 1:num_reactions
            if i == j
                add_edge!(graph, i, j)
            else
                if !isempty(intersect(reactions[i].affects, reactions[j].depends_on))
                    add_edge!(graph, i, j)
                end
            end
        end
    end
    return graph
end
