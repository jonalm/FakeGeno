function simplecorr(pop::AbstractMatrix{Bool}, i,j)
    Ncol = size(pop)[2]
    pi = 1.* sum(pop[i,:]) / Ncol
    pj = 1.* sum(pop[j,:]) / Ncol
    pij = 1.* sum(pop[j,:].*pop[i,:]) / Ncol
    (pij - pi*pj) / sqrt(pi*(1-pi) * pj *(1-pj))
end

function calc_dosage(pop::AbstractMatrix{Bool})
    Nsnp = size(pop)[1]
    Nind = Int(size(pop)[2]/2)
    dose = zeros(Int, Nsnp, Nind)
    for i in 1:Nind
        dose[:,i] = pop[:,2*i-1] + pop[:,2*i]
    end
    dose
end

function zygotefreq(dosage)
    Nsnp = size(dosage)[1]
    Nind = size(dosage)[2]
    res = zeros(Float64, 3, Nsnp)
    for i in 1:Nsnp
        res[:, i] = counts(dosage[i,:],0:2) / Nind
    end
    res'
end

function heterozygote(pop::AbstractMatrix{Bool})
    zf = zygotefreq(calc_dosage(pop))
    pHW = sqrt(zf[:,1])
    Float64[zf[:,2] 2*pHW.*(1-pHW)]
end

z2p(zscores) = [ ccdf(Normal(), abs(z)) for z in zscores]
removenan(vec) = vec[~isnan(vec)]
