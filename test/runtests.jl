using Test
using NextReaction

@testset "simulate!" begin
    # Test that simulate! runs without error
    rs = parse_input_lines([
        "A B C",
        "A + B -> C  0.1",
        "C -> A + B  0.2",
    ])
    simulate!(rs)
end