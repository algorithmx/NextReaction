function sample_tau_and_reaction_1(p::Vector{Float64})::Tuple{Float64,Int}
    # p is the propensity vector
    a0 = sum(p)
    # putatiVe time
    tau = (-1 / a0) * log(rand())
    sum_prop = 0.0
    r1 = rand() * a0
    for (mu, a) in enumerate(p)
        sum_prop += a
        if sum_prop >= r1
            return tau, mu
        end
    end
    return tau, length(p)
end
