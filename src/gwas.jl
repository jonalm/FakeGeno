function GWAS{T<:Real}(dosage::Matrix{T}, cases::AbstractVector{Bool})
    Nsnp, Nind = size(dosage)
    @assert length(cases)==Nind
    
    dosType = eltype(dosage)
    dosT = dosage'
    zscores = Vector{Float64}(Nsnp)
    one_vec = ones(dosType, Nind,1)
    @showprogress 1 "calculating GWAS z-scores ..." for i in 1:Nsnp
        try
            logit = glm(dosType[one_vec dosT[:,i]], cases, Binomial(), LogitLink())
            zscores[i] = coef(logit)[2] / stderr(logit)[2]
        catch err
            if isa(err, Base.LinAlg.PosDefException)
                zscores[i]=NaN
            elseif isa(err, LoadError)
                zscores[i]=NaN
            else
                throw(err)
            end
        end
    end
    zscores
end

function nullpermuteGWAS(pop::AbstractMatrix{Bool}, prevalence::Float64, Npermute::Int)
    Nind = Int(size(pop)[2]/2)
    Nsnp = size(pop)[1]
    dosage = calc_dosage(pop)

    zmat = Array(Float64, Nsnp, Npermute)
    @showprogress 1 "looping over case/control permutations ..." for i in 1:Npermute
        cases = rand(Nind).<prevalence;
        zscores = GWAS(dosage, cases)
        zmat[:, i] = zscores
    end
    zmat
end
