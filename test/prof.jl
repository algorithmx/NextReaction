using Pkg
Pkg.activate("/home/dabajabaza/Documents/Research/KMC/NextReaction/NextReaction")
using NextReaction

using Profile
#using ProfileView
using PProf

Profile.clear()
Profile.Allocs.clear()


#ProfileView.@profview

rs = Simulator(parse_input_lines(readlines("./test/circle.chem"))...)
@time simulate_first_reaction!(rs, 10000.0)

rs = Simulator(parse_input_lines(readlines("./test/circle.chem"))...)
@time Profile.Allocs.@profile sample_rate = 4 simulate_first_reaction!(rs, 100.0)

PProf.Allocs.pprof(from_c=false)