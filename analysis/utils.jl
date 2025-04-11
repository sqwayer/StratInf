using Distributions
""" Standard error of the mean """
sem(X) = std(X)/sqrt(length(X))

""" Stats on a sliding window """
function movstat(stat, X, win)
    """Computes stat of vector X on a moving window win"""
    R = similar(X)
    for i in eachindex(X)
        idx = max(1, i-win):min(i+win,length(X))
        R[i] = stat(X[idx])
    end
    return R
end

""" Mutual information """
function mutual_info(X::Vector{T},win::Int) where T
    N = length(X)
    MI = zeros(N-1)
    #PJ = Dict{String, Vector{Float64}}("00" => Float64[], "01" => Float64[], "10" => Float64[], "11" => Float64[])
    for n = 1:N-1
        idx1 = max(1,n-win):min(N-1, n+win)
        idx2 = idx1 .+ 1
        MI[n], pp = mutual_info(X[idx1], X[idx2])
        # for (k, v) in pp
        #     push!(PJ[k], v)
        # end
    end
    return MI
end

function mutual_info(X::Vector{T}, Y::Vector{T}; allX = unique(X), allY = unique(Y)) where T
    minval = sqrt(nextfloat(0.0))
    MI = 0.
    PJ = Dict{String, Float64}()#zeros(length(allX)*length(allY))
    base = 2#length(allX)
    for x in allX
        pX = max(mean(X .== x), minval)
        for y in allY
            pY = max(mean(Y .== y), minval)
            pJoint = max(mean( (X .== x) .& (Y .== y) ), minval)
            MI += pJoint * log(base, pJoint / pX / pY)
            #push!(PJ, pJoint)
            PJ[string(x, y)] = log(base, pJoint / pX / pY)
        end
    end
    return MI, PJ
end

""" Hodges Lehmann estimator """
function hodges_lehmann(X)
    n = length(X)
    Y = eltype(X)[]
    for i in 1:n-1
        for j in i+1:n
            push!(Y, (X[i] + X[j])/2)
        end
    end
    return Y
end

""" Cluster-based permutation test """
function cluster_perm_test(X; niter=1e6)
    Tstats, clusters, posclusters, negclusters = cluster_tstat(X)
    pvalues = zero(Tstats)

    # Random permutations
    flip = [-1, 1]
    shuffledX = copy(X)
    for i = 1:niter
        for j in eachindex(shuffledX)
            shuffledX[j] *= rand(flip)
        end
        shuffledTstats, _, _, _ = cluster_tstat(shuffledX)
        for j in eachindex(pvalues)
            pvalues[j] += maximum(abs, shuffledTstats) > abs(Tstats[j])
        end  
    end
    pvalues ./= niter

    # Attribute pvalues / significativity to each individual trials
    trialpval = ones(size(X, 2))
    trialsig = zeros(Int, size(X, 2))
    for i in eachindex(clusters)
        trialpval[clusters[i]] .= pvalues[i]
        trialsig[clusters[i]] .= pvalues[i] < 0.05
    end

    return (trialsig=trialsig, trialpval=trialpval, Tstats=Tstats, pvalues=pvalues, clusters=clusters)
end

function cluster_tstat(X)
    # X : Matrix df x trials
    df, trials = size(X)

    # Pre-allocate
    tstats = zeros(trials)
    posvals = zeros(trials)
    negvals = zeros(trials)

    for i in 1:trials
        actualdf = count(!isnan, view(X, :,i)) 
        tstats[i] = mean(filter(!isnan, view(X, :,i))) / sqrt(var(filter(!isnan, view(X, :,i))) / actualdf)
        posvals[i] = (pdf(TDist(actualdf), tstats[i]) < 0.05) && (tstats[i] > 0.0)
        negvals[i] = (pdf(TDist(actualdf), tstats[i]) < 0.05) && (tstats[i] < 0.0)
    end

    # Clustering
    posclusters = find_clusters(posvals)
    negclusters = find_clusters(negvals)

    clusters = vcat(posclusters, negclusters)

    # cluster mass 
    Tstats = [sum(tstats[ci]) for ci in clusters]

    return Tstats, clusters, posclusters, negclusters
end

function find_clusters(V)

    C = [Int[]]
    j = 1
    for i in eachindex(V)
        if V[i] == 1
            push!(C[j], i)
        elseif i > 1 && V[i-1] == 1
            j += 1
            push!(C, Int[])
        end
    end
    if isempty(C[end])
        deleteat!(C, length(C))
    end
    return C
end

function shuffle_data(X)
    flip = [-1, 1]
    shuffledX = copy(X)
    for i in eachindex(shuffledX)
        shuffledX[i] *= rand(flip)
    end
    return shuffledX
end