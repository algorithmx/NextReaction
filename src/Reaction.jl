ReactantId = Int

abstract type AbstractReaction end

mutable struct SimpleReaction <: AbstractReaction
    # assume that all rate constants
    rate_constant::Float64
    reactants::Vector{ReactantId}
    products::Vector{ReactantId}
    affects::Set{ReactantId}
    depends_on::Set{ReactantId}
end

mutable struct Reaction{T<:Real} <: AbstractReaction
    # assume that all rate constants
    rate_constant::Float64 
    reactants::Dict{ReactantId,T}
    products::Dict{ReactantId,T}
    affects::Set{ReactantId}
    depends_on::Set{ReactantId}
end

@inline function disjoint_union(A, B)::Set
    I = intersect(Set{Int}(A), Set{Int}(B))
    U = union(Set{Int}(A), Set{Int}(B))
    return Set{Int}(setdiff(U, I))
end

@inline Affects(r::Reaction) = disjoint_union(keys(r.reactants), keys(r.products))
@inline Affects(r::SimpleReaction) = disjoint_union(r.reactants, r.products)

@inline DependsOn(r::Reaction) = Set(keys(r.reactants))
@inline DependsOn(r::SimpleReaction) = Set(r.reactants)

add_information!(r::AbstractReaction) = (r.affects = Affects(r); r.depends_on = DependsOn(r); r)
