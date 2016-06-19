module FakeGeno

using ProgressMeter
using Roots
using Distributions
using StatsBase
using PyPlot
using GLM

export
cor, #from StatsBase
calc_mrate,
create_pop_random,
make_pop,
evolve_pop_once!,
evolve_pop!,
calc_mrate,
showmatrix,
manhattanplot,
zhist,
calc_dosage,
heterozygote,
GWAS,
nullpermuteGWAS,
z2p

include("initial_pop.jl")
include("evolve_pop.jl")
include("plots.jl")
include("utility.jl")
include("gwas.jl")

end
