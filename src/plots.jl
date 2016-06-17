function showmatrix(matrix::AbstractMatrix; minmax=[0,1], cmap="binary")
    fig = plt[:figure]()
    ax = fig[:add_axes]([0,0,1,1])
    ax[:imshow](matrix,interpolation="nearest",cmap=cmap, vmin=minmax[1], vmax=minmax[2])
    ax[:axis]("off")
    fig
end

function zhist(zscores::AbstractVector)
    x = collect(linspace(-3.5,3.5,40))
    BLUE = "#348ABD"
    RED = "#A60628"
    fig = plt[:figure](figsize=(5,3))
    ax = fig[:add_subplot](111)
    ax[:hist](zscores, bins=x, color=BLUE, ec="white",normed=1)
    ax[:plot](x, pdf(Normal(),x),"--",color=RED, lw=2)
    ax[:set_ylim]([0,0.5])
    ax[:set_xlim]([x[1],x[end]])
    fig
end

function manhattanplot(zscores::AbstractVector, snpnumber::AbstractVector)
    @assert length(zscores) == length(snpnumber)
    RED = "#A60628"
    log = log10

    neglogp = -log(calcpvals(zscores))
    Bonf1 = -log(0.05/length(snpnumber))
    Bonf2 = -log(0.01/length(snpnumber))
    range = [snpnumber[1],snpnumber[end]]

    fig = plt[:figure](figsize=(5,3))
    ax = fig[:add_subplot](111)
    ax[:plot](snpnumber, neglogp, "xk")
    ax[:plot](range, [Bonf2, Bonf2], "--", color=RED)
    ax[:plot](range, [Bonf1, Bonf1], ":", color=RED)

    ax[:set_ylim]([0,9])
    ax[:set_xlim](range)
    fig
end
