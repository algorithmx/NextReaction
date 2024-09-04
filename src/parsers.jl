function parse_input_lines(
        lines::Vector{String}
)::Tuple{Dict{Symbol,ReactantId},Vector{AbstractReaction}, Vector{Float64}}
    # Read the reactants from the first line
    reactant_symbols = split(lines[1], r"\s+")
    reactants = Dict{Symbol,ReactantId}(
        Symbol(symbol)=> i 
        for (i, symbol) in enumerate(reactant_symbols))
    
    initial_state = parse.(Float64, split(lines[2], r"\s+"))

    # Read each line as a reaction
    reactions = []
    for line in lines[3:end]
        # Parse the reaction line
        md, tp, rate_constant, reactant_symbols, product_symbols = parse_line(line)
        # Convert reactant and product symbols to reactant IDs
        if md == :simple
            r_id = ReactantId[reactants[Symbol(s)] for s in reactant_symbols]
            p_id = ReactantId[reactants[Symbol(s)] for s in product_symbols]
            reaction = SimpleReaction(
                rate_constant, r_id, p_id,
                Set{ReactantId}(), Set{ReactantId}())
            add_information!(reaction)
        else
            @assert md == :coeff
            reaction = Reaction(
                rate_constant,
                Dict{ReactantId,tp}(reactants[Symbol(s)] => f for (s,f) ∈ reactant_symbols),
                Dict{ReactantId,tp}(products[Symbol(s)] => f for (s,f) ∈ product_symbols),
                Set{ReactantId}(),
                Set{ReactantId}())
            add_information!(reaction)
        end

        # Add the reaction to the list
        push!(reactions, reaction)
    end
    return reactants, reactions, initial_state
end


function parse_line(line::AbstractString)
    reactant_str, product_str = split(line, " -> ")
    reactant_symbols = split(reactant_str, " + ")
    product_symbols = split(product_str, " + ")
    rate_constant_str = split(product_symbols[end], " [")[2]
    rate_constant = parse(Float64, rate_constant_str[1:end-1])
    product_symbols[end] = split(product_symbols[end], " [")[1]

    reactants = Dict{String,Any}()
    products = Dict{String,Any}()
    Typ = Int
    for x in reactant_symbols
        if occursin(r"^\d+.\d+ ", x)
            Typ = Float64
            c, r = split(x, r"\s+")
            reactants[r] = parse(Typ, c)
        elseif occursin(r"^\d+ ", x)
            c, r = split(x, r"\s+")
            reactants[r] = parse(Typ, c)
        else
            reactants[x] = 1
        end
    end
    for symbol in product_symbols
        if occursin(r"^\d+.\d+ ", symbol)
            Typ = Float64
            c, p = split(symbol, r"\s+")
            products[p] = parse(Typ, c)
        elseif occursin(r"^\d+ ", symbol)
            c, p = split(symbol, r"\s+")
            products[p] = parse(Typ, c)
        else
            products[symbol] = 1
        end
    end
    md = :simple
    tp = :int
    if Typ == Int
        tp = :int
        if all(x->x==1, values(reactants)) && all(x->x==1, values(products))
            md = :simple
            return md, tp, rate_constant, collect(keys(reactants)), collect(keys(products))
        else
            md = :coeff
            return md, tp, rate_constant, Dict{String,Int}(reactants), Dict{String,Int}(products)
        end
    else
        @assert Typ == Float64
        tp = :float
        md = :coeff
        return md, tp, rate_constant, Dict{String,Float64}(reactants), Dict{String,Float64}(products)
    end
end