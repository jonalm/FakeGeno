function calc_dosage(pop::AbstractMatrix{Bool})
    Nsnp = size(pop)[1]
    Nind = Int(size(pop)[2]/2)
    dose = zeros(Int8, Nsnp, Nind)
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

z2p(zscore::Float64) = 2*ccdf(Normal(), abs(zscore))
z2p(zscore::Vector{Float64}) = Float64[z2p(z) for z in zscore]
z2nlp(zscore::Float64) = -log10(z2p(zscore))
z2nlp(zscore::Vector{Float64}) = -log10(z2p(zscore))

p2z(p::Float64, s::Float64=1.0) = -norminvcdf(p/2) * sign(s)
nlp2z(nlp::Float64, s::Float64=1.0) = p2z(10.^(-nlp), s)

# https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
CHRLENGHThg19 = cumsum(Int[0,
                           249250621, # chr1 lenght
                           243199373,
                           198022430,
                           191154276,
                           180915260, # chr5 length
                           171115067, 
                           159138663,
                           146364022,
                           141213431,
                           135534747, # chr10 lenght
                           135006516,
                           133851895,
                           115169878,
                           107349540,
                           102531392, # chr15 lenght
                           90354753,
                           81195210,
                           78077248,
                           59128983,
                           63025520,  # chr20 length
                           48129895,
                           51304566])

function cbp2abspos(c::Integer, bp::Integer, chrlength::Vector{Int}=CHRLENGHThg19)
    chrlength[c]+bp
end

function abspos2cbp(abspos::Integer, chrlength::Vector{Int}=CHRLENGHThg19)
    temp = findfirst(chrlength .>= abspos, true) - 1
    chr = temp < 0 ? 23 : temp
    (chr, abspos - chrlength[chr])
end

function findsorted{T<:Any}(A::Vector{T}, v::T)
    I = (0, length(A)) # interval,
    while true
        L = I[2] - I[1] # interval length
        if L<=1
            return A[I[2]] == v ? I[2] : -1
        end
        half = I[1] + L >>> 1
        I = v<=A[half] ? (I[1], half) : (half, I[2])
    end
end

function readrow(file::HDF5.HDF5File,ir_objname::String,jc_objname::String,colnum::Int)
    # read row indices from a CSC sparse matrix format (compressed sparse column matrix)
    # note that ir and jc in HDF5 file is assumed to be indexed from 0
    # rowindices = ir[jc[columnindex]:jc[columnindex+1]-1]
    jc1, jc2 = file[jc_objname][colnum:(colnum+1)]
    Int[i+1 for i in file[ir_objname][Int(jc1+1):Int(jc2)]]
end

function extractmat(matfn::String, regex::Regex=r""; flat::Bool=false)
    flatten(x) = flat ? reshape(x, prod(size(x))) : x
    matopen(matfn) do file
        hits = String[]
        for n in names(file)
            ismatch(regex, n) && push!(hits,n)
        end
        length(hits)==0 && error("$regex not found")
        length(hits)>1 && error("no unique field, specify by regex\nFound : $(join(hits,','))\n")
        read(file, hits[1]) |> flatten
    end
end
