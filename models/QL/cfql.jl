""" Counterfactual Q-learning """

function cfql(nactions, nstates)
    @model Counterfactual_QL(na, ns) = begin
        logitα ~ Normal()
        logitκ ~ Normal()
        β ~  truncated(LogNormal(), 0.0, 50.0)#
        ϵ ~ Beta() 
        
        return (α = logistic(logitα), β = β, ϵ = ϵ, κ = logistic(logitκ),
                na = na,
                ns = ns,
                Values = Matrix{typeof(logitα)}(fill(1/na, na, ns)),
                Qactions = Vector{typeof(logitα)}(fill(1/na, na))) # pre-allocates the Qvalues to speed up the observation step
    end

    mdl = Counterfactual_QL(nactions, nstates)

    @evolution mdl begin
        for i = 1:na
            if i == a
                delta_rule!(Values, s, a, r, α)
            else
                delta_rule!(Values, s, i, 1.0 - r, κ*α)
            end
        end
    end

    @observation mdl begin
        Qactions .= @views(Values[:,s])
        Qactions .*= β
        Categorical(epsilon_greedy!(softmax!(Qactions), ϵ))
    end

    return mdl
end



## Hierarchical 
function generate_latent(::Val{:Counterfactual_QL}, na, ns ; α, β, ϵ, κ)
    return  (α = α, β = β, ϵ = ϵ, κ = κ,
            na = na,
            ns = ns,
            Values = Matrix{typeof(α)}(fill(1/na, na, ns)),
            Qactions = Vector{typeof(α)}(fill(1/na, na)))
end
@model function cfql_hierarchical(A, data)
    nd = length(data)
    μ ~ MvNormal(4, 1.0)
    #τ ~ MvNormal(4, 1.0)
    τ = zeros(4)
    logitα ~ MvNormal(nd, 1.0)
    logitκ ~ MvNormal(nd, 1.0)
    logβ ~ MvNormal(nd, 1.0)
    logitϵ ~ MvNormal(nd, 1.0)

    mdl = cfql(3,3)
    for n = 1:nd
        
        α=logistic(μ[1] + exp(τ[1]) * logitα[n])
        κ=logistic(μ[2] + exp(τ[2]) * logitκ[n])
        β=exp(μ[3] + exp(τ[3]) * logβ[n])
        ϵ=logistic(μ[4] + exp(τ[4]) * logitϵ[n])
         
        θ = generate_latent(Val(mdl.name), 3, 3; α=α, β=β, ϵ=ϵ, κ=κ)
        P = [AnimalBehavior.cycle!(θ, mdl, obs) for obs in data[n]]
        A[n] ~ arraydist(P)
    end    
end
