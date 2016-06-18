function GWAS(pop::AbstractMatrix{Bool}, casecontrol::AbstractVector{Bool})
    dosT = calc_dosage(pop)' # each column corresponds to snp
    (Nind, Nsnp) = size(dosT)
    one_vec = ones(Nind)

    zscores = Vector{Float64}(Nsnp)

    @showprogress 1 "calculating GWAS z-scores ..." for i in 1:Nsnp
        try
            logit = glm(Float64[one_vec dosT[:,i]], casecontrol, Binomial(), LogitLink())
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

    zmat = Array(Float64, Nsnp, Npermute)
    @showprogress 1 "looping over case/control permutations ..." for i in 1:Npermute
        casecontrol = rand(Nind).<prevalence;
        zscores = GWAS(pop, casecontrol)
        zmat[:, i] = zscores
    end
    zmat
end
