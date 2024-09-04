function sample_tau_and_reaction_2(L::Vector{Float64})::Tuple{Float64,Int}
    # p is the propensity vector
    min_tau = 1e100
    min_idx = 0
    for (i,p) in enumerate(L)
        tau = (-1 / p) * log(rand())
        if tau < min_tau
            min_tau = tau
            min_idx = i
        end
    end
    return min_tau, min_idx
end
