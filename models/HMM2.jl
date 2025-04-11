using Turing, StatsFuns, LinearAlgebra
const hidden_states = [
    [1,2,3],
    [1,3,2],
    [2,1,3],
    [2,3,1],
    [3,1,2],
    [3,2,1],
    [1,1,2],
    [1,2,1],
    [2,1,1],
    [1,1,3],
    [1,3,1],
    [3,1,1],
    [2,2,1],
    [2,1,2],
    [1,2,2],
    [3,3,1],
    [3,1,3],
    [1,3,3],
    [2,2,3],
    [2,3,2],
    [3,2,2],
    [3,3,2],
    [3,2,3],
    [2,3,3],
    [1,1,1],
    [2,2,2],
    [3,3,3]
]
nh = length(hidden_states)
tmp = fill(4, nh+1, nh+1)
for i = 1:nh, j = 1:nh
    tmp[i,j] = sum(hidden_states[i] .≠ hidden_states[j])
end
const Htrans = tmp
const Hidx = [findall(Htrans .== i) for i = 0:4]
const Hn = [sum(Htrans[:,1] .== i) for i = 0:4] 
##
function build_trans_mat(θ, τ)
    nh = length(hidden_states)
    p = logistic.(θ)
    p[2:5] .*=  (0.5 - 0.5*p[1]) / sum(p[2:5])
    p[1] = 0.5 * (1 + p[1])

    for i = 2:5
        p[i] /= Hn[i]
    end
    l = log.(p)
    Ltrans = getindex.([l], Htrans .+ 1)
    Ltrans[end,:] .= log((1 - τ)/nh) # Probability of leaving exploration
    Ltrans[end,end] = log(τ) # Probability of staying in exploration
    return Ltrans
end

function get_probas(θ)
    p = logistic.(θ)
    p[2:5] .*=  (0.5 - 0.5*p[1]) / sum(p[2:5])
    p[1] = 0.5 * (1 + p[1])

    return p[2:5]
end

function get_llh(::Val{:obs}, ρ, hs, s, a, r)
    return xor(hs[s] == a, r==1) ? log(ρ) : log(1 - ρ)
end

function get_llh(::Val{:act}, ρ, hs, s, a, r)
    return hs[s] ≠ a ? log(ρ/2) : log(1 - ρ)
end

@model HMM2(actobs, S,A,R) = begin
    θ ~ MvNormal(5, 1.0)
    τ ~ Beta()
    ρ ~ truncated(Beta(), 0.0, 0.5)

    nh = length(hidden_states)
    Lprior = Vector{typeof(ρ)}(fill(log(1/nh+1), nh+1))
    LLH = zeros(typeof(ρ),nh+1)
    Lpost = Vector{typeof(ρ)}(fill(log(1/nh+1), nh+1))
    bestPrior = Array{Int, 2}(undef, length(S), nh+1)

    Ltrans = build_trans_mat(θ, τ)
    
    L = zeros(typeof(ρ),length(S),nh+1)

    tmp_prior = zeros(typeof(ρ), size(Ltrans))
    tmp_post = zeros(typeof(ρ), size(Ltrans))
    for t in eachindex(S)
        s = S[t]
        a = A[t]
        r = R[t]
        tmp_prior .= Ltrans 
        tmp_prior .+= Lprior
        tmp_post .= Ltrans
        tmp_post .+= Lpost
  
        # Forward pass
        for hi in eachindex(hidden_states)
           Lprior[hi] = logsumexp(@views(tmp_prior[:,hi]))
            LLH[hi] = get_llh(actobs, ρ, hidden_states[hi], s, a, r)
            bv, bestPrior[t,hi] = findmax(@views(tmp_post[:,hi]))
            Lpost[hi] = bv + LLH[hi]
        end
        # Exploration state
        Lprior[end] = logsumexp(@views(tmp_prior[:,end]))
        if actobs == Val(:obs)
            LLH[end] = r == 1 ? log(1/3) : log(2/3) 
        else
            LLH[end] = log(1/3)
        end
        bv, bestPrior[t,end] = findmax(@views(tmp_post[:,end]))
        Lpost[end] = bv + LLH[end]

        # Update prior
        Lprior .+= LLH
        margL = logsumexp(Lprior)
        Lprior .-= margL
        Turing.@addlogprob! margL
    end

    # Backtracking 
    H = Vector{Int}(undef, length(S))
    H[end] = argmax(Lpost) # The last state is the one with the full trajectory's MAP
    for t = length(S)-1:-1:1 # For each previous trial
        H[t] = bestPrior[t+1, H[t+1]] # the MAP state is the best prior state for the next trial
    end
    return H
end