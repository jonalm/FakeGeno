function create_bitpop_triange(Nsnp::Int, Nind::Int)
    res = bitpack(zeros(Bool, Nsnp, 2Nind))
    for i in 1:2Nind
        for j in 0:(i-1)
            res[j%Nsnp+1,i] = !res[j%Nsnp+1,i]
        end
    end
    res
end

function create_pop_banded(Nsnp::Int, Nind::Int)
    res = Array(Bool,Nsnp,2Nind)
    for i in 1:size(res)[2]
        if i%2 == 0
            res[:,i]=true
        else
            res[:,i]=false
        end
    end
    res
end

function create_pop_random(Nsnp::Int, Nind::Int)
    bitunpack(rand(Nsnp, 2Nind).<0.5)
end
