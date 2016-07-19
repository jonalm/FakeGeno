function mate!(pop::AbstractMatrix{Bool}, selected_alleles::Vector{Int}, sortedunique::Vector{Int})
    temp = pop[:,sortedunique]
    for (i, j) in enumerate(selected_alleles)
        tempidx = findfirst(sortedunique, j)
        pop[:,i] = temp[:,tempidx]
    end
    pop
end

function crossover_chr!(pop::AbstractMatrix{Bool}, ind::Int, crossoverindex::Int)
    temp = pop[crossoverindex:end,ind]
    pop[crossoverindex:end,ind] = pop[crossoverindex:end,ind+1]
    pop[crossoverindex:end,ind+1] = temp
    return pop
end

function crossover!(pop::AbstractMatrix{Bool}, hotspots::Vector{Int}, sortedunique::Vector{Int})
    for i in unique(Int[iseven(j)? j-1:j for j in sortedunique])
        crossoverindex = sample(hotspots)
        crossover_chr!(pop, i, crossoverindex)
    end
end

function mutate_column!(pop::AbstractMatrix{Bool}, column::Int, mutationrate::Float64)
    Nsnp = size(pop)[1]
    selected = find(rand(Nsnp) .< mutationrate)
    for s in selected
        pop[s,column] = !pop[s,column]
    end
    return pop
end

function mutate!(pop::AbstractMatrix{Bool}, mutationrate::Float64, sortedunique::Vector{Int})
    Ncolumns = size(pop)[2]
    for i in sortedunique
        mutate_column!(pop, i, mutationrate)
    end
    return pop
end

function evolve_pop_once!(pop::AbstractMatrix{Bool};
                          hotspots::Vector{Int}=Int[], mutationrate::Float64=-1.0)
    Ncol = size(pop)[2]
    selected_alleles = sample(1:Ncol,Ncol)
    sortedunique = sort!(unique(selected_alleles))

    if length(hotspots)>0
        crossover!(pop, hotspots, sortedunique)
    end
    if mutationrate>0
        mutate!(pop, mutationrate, sortedunique)
    end
    mate!(pop, selected_alleles, sortedunique)
end

function evolve_pop!(pop::AbstractMatrix{Bool}, Niter::Int;
                     hotspots::Vector{Int}=Int[], mutationrate::Float64=-1)
    @showprogress 1 "evolving population ..." for i in 1:Niter
        evolve_pop_once!(pop, hotspots=hotspots, mutationrate=mutationrate)
    end
end

# asymptotic prob for individual having two different alleles
function calc_mrate(Nind::Int, heterozygozity::Float64)
    heterozygote_pred(Nind, mrate) =  2(mrate*(1-mrate)) / (1 - 4*(mrate-0.5)^2 * (1-1/(2Nind)))
    obj(x) = heterozygote_pred(Nind, x) - heterozygozity
    fzero(obj,0.0,0.9)
end

function make_all_snps_minor!(pop::Matrix)
    select = sum(pop,2) / size(pop)[2] .> 0.5
    pop[select,:] = ~pop[select,:]
    pop
end

function remove_snps_under_maf_threshold(pop::Matrix, threshold::Float64)
    select = (sum(pop,2) / size(pop)[2]) .> threshold
    pop[select,:]
end

function make_pop(Nsnp::Int, Nind::Int, Nhs::Int;
                  hzygosity::Float64=0.4,
                  maf_cutoff::Float64=0.05,
                  initials_sample::Function=create_pop_random)
    @assert Nhs < (Nsnp-1)
    mrate = calc_mrate(Nind, hzygosity)
    Niter = 4*2*Nind

    println("================================================")
    println("generate sample population ....")
    println()
    println("# inidviduals  : $Nind")
    println("# chromosomes  : $(2*Nind)")
    println("# SNPs         : $Nsnp")
    println("# generations  : $Niter   ($(Niter/(2Nind)) x char. decay)")
    println("# hotspots     : $Nhs")
    println("# maf cutoff   : $maf_cutoff")
    println("# mutationrate : $mrate")
    println("# approx hz    : $hzygosity")
    println()
    println("================================================")

    pop = initials_sample(Nsnp, Nind)
    hotspots = sample(2:Nsnp,Nhs)
    evolve_pop!(pop, Niter, hotspots=hotspots, mutationrate=mrate)
    make_all_snps_minor!(pop)
    pop = remove_snps_under_maf_threshold(pop, maf_cutoff)
    println("# individuals after cutoff  : $(size(pop)[2])")
    println("================================================")
    pop
end

