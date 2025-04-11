using Distributions, SpecialFunctions, LinearAlgebra, Printf, JLD2, FileIO
maxval = 700.0
function BMS(paths, models; criterion=:BIC)
    nModels=4
    BIC = [[] for k in 1:nModels]
    FFX = zeros(nModels)
    ## Load fits
    for path in paths
        ## BIC and FFX
        for k = 1:nModels
            ModelResult = load("$path/$(models[k]).jld2", "Result")
            append!(BIC[k], ModelResult[criterion])
            FFX[k] += sum(ModelResult[criterion])
        end
        
    end
    nSub = length(BIC[1])
    BIC = hcat(BIC...)
    XP, MAtt, MFq = exceedance_probability(-BIC)
    Result = Dict{Symbol, Any}(
    :nSub => nSub,
    :nModels => nModels,
    :labels => models,
    :BIC => BIC,
    :FFX => FFX,
    :MAttributions => MAtt,
    :MFq => MFq,
    :XP => XP
    )
    return Result
end

function BMS(path::String; criterion=:BIC)

    ## Load fits
    filenames = filter(x -> occursin(".jld2", x), readdir(path))
    nModels = length(filenames)

    ## BIC and FFX
    BIC = Array{Array{Float64, 1}, 1}(undef, nModels)
    FFX = zeros(nModels)
    for k = 1:nModels
        ModelResult = load("$path/$(filenames[k])", "Result")
        BIC[k] = ModelResult[criterion]
        FFX[k] = sum(ModelResult[criterion])
    end
    nSub = length(BIC[1])
    BIC = hcat(BIC...)

    ## Exceedance probability
    XP, MAtt, MFq = exceedance_probability(-BIC)

    ## Return result
    Result = Dict{Symbol, Any}(
    :nSub => nSub,
    :nModels => nModels,
    :labels => [f[1:end-5] for f in filenames],
    :BIC => BIC,
    :FFX => FFX,
    :MAttributions => MAtt,
    :MFq => MFq,
    :XP => XP
    )
    return Result
end

function exceedance_probability(logME; famillies = [])
    ## From Stephan et al. 2009
    alpha, gnk = variational_dirichlet(logME)
    nModels = size(logME, 2)
    if !isempty(famillies)
        tmp = zeros(length(famillies))
        for i in eachindex(famillies)
            tmp[i] = sum(alpha[famillies[i]])
        end
        alpha = tmp
        nModels = length(famillies)
    end
    
    EP = zeros(nModels)
    r = zeros(nModels, 2) # Mean and variance
    a0 = sum(alpha)
    Ea = alpha ./ a0
    r[:,1] = Ea
    r[:,2] = Ea .* (1.0 .- Ea) ./ (a0 + 1)
    if nModels == 2
        EP[1] = cdf(Beta(alpha[2], alpha[1]), 0.5)
        EP[2] = 1.0 - EP[1]
    else
        nSamples = 1e6
        for i = 1:nSamples
            q = rand(Dirichlet(alpha))
            for j = 1:nModels
                EP[j] += q[j] == maximum(q)
            end
        end
        EP ./= nSamples
    end
    return EP, gnk, r
end

function variational_dirichlet(logME)
    cc = 1e-3
    c = 1.0
    nSub, nModels = size(logME)
    logU = zeros(nSub, nModels)
    gnk = zeros(nSub,nModels)
    alpha0 = ones(nModels)
    alpha = copy(alpha0)
    while c > cc
        for n = 1:nSub
            for k = 1:nModels
                logU[n,k] = logME[n,k] + digamma(alpha[k]) - digamma(sum(alpha))
            end
            logU[n,:] .-= mean(logU[n,:])
            for k = 1:nModels
                logU[n,k] = sign(logU[n,k]) * min(maxval,abs(logU[n,k]))
            end
        end
        U = exp.(logU)
        gnk .= U ./ sum(U,dims=2)
        beta = zeros(nModels)
        for k = 1:nModels
            beta[k] = sum(gnk[:,k])
        end
        prev = copy(alpha)
        alpha = alpha0 + beta
        c = norm(alpha - prev)
    end
    return alpha, gnk
end

function compute_bor(logME, alpha, gnk)

    nSub, nModels = size(logME)
    alpha0 = sum(alpha)
    F0 = 0.0
    F1 = 0.0

    Elogr = digamma.(alpha) .- digamma(alpha0)
    Sqf = sum(loggamma.(alpha)) - loggamma(alpha0) - sum((alpha .- 1) .* Elogr)
    Sqm = 0.0
    ELJ = loggamma(nModels) - nModels * loggamma(1.0)# + (1 - 1) * sum(Elogr)#sum((alpha0 - 1) .* Elogr)

    for n = 1:nSub
        normL = logME[n,:] .- maximum(logME[n,:])
        softmax!(normL)
        Sqm -= sum(gnk[n,:] .* log.(gnk[n,:] .+ 0.0))

        for k = 1:nModels
            F0 += normL[k] * (logME[n,k] - log(normL[k]) - log(nModels))
            ELJ += gnk[n,k] * (Elogr[k] + logME[n,k])
        end
    end
    F1 = ELJ + Sqf + Sqm
    bor = 1/(1 + exp(F1 - F0))
    return bor, F0, F1
end