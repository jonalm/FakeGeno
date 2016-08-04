module GWASTools

using ProgressMeter
using Roots
using Distributions
using StatsBase
using StatsFuns
using PyPlot
using GLM
using HDF5
using MAT

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
probhist,
uniform_quantile,
calc_dosage,
heterozygote,
GWAS,
nullpermuteGWAS,
z2p,
z2nlp,
p2z,
nlp2z,
cbp2abspos,
abspos2cbp,
findsorted,
readrow,
extractmat

include("initial_pop.jl")
include("evolve_pop.jl")
include("plots.jl")
include("utility.jl")
include("gwas.jl")

end
