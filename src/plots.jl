include("mycolors.jl")

function showmatrix(matrix::AbstractMatrix; minmax=[0,1], cmap="binary")
    fig = plt[:figure]()
    ax = fig[:add_axes]([0,0,1,1])
    ax[:imshow](matrix,interpolation="nearest",cmap=cmap, vmin=minmax[1], vmax=minmax[2])
    ax[:axis]("off")
    fig
end

function zhist{T<:Real}(zscores::Vector{T}; lim=3.0)
    x = collect(linspace(-3.5,3.5,40))
    fig = plt[:figure](figsize=(5,3))
    ax = fig[:add_subplot](111)
    ax[:hist](zscores, bins=x, color=BLUE, ec="white",normed=1)
    ax[:plot](x, pdf(Normal(),x),"--",color=RED, lw=2)
    ax[:set_ylim]([0,0.5])
    ax[:set_xlim]([-lim,lim])
    fig
end

function manhattanplot{T<:Real}(zscores::Vector{T}, highlight::AbstractVector{Bool}=Bool[])
    N = length(zscores)
    if length(highlight) != N
        highlight=BitArray(N)
        highlight[:] = false
    end

    neglogp = -log10(z2p(zscores))
    snpnumber = collect(1:N)
        

    range = [0,N+1]
    Bonf1 = -log10(0.05/N)
    Bonf2 = -log10(0.01/N)
    
    fig = plt[:figure](figsize=(5,3))
    ax = fig[:add_subplot](111)
    ax[:plot](snpnumber[~highlight], neglogp[~highlight], "xk")
    ax[:plot](snpnumber[highlight], neglogp[highlight], "o", color=RED)

    ax[:plot](range, [Bonf2, Bonf2], "--", color=GREEN)
    ax[:plot](range, [Bonf1, Bonf1], ":", color=GREEN)
    ax[:set_xlim](range)
    fig
end

