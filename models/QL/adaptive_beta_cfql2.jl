""" Counterfactual Q-learning with adaptive inverse temperature following Khamassi et al. 2018 """

function adaptive_beta_cfql2(nactions, nstates)
    @model DynamicBeta_Counterfactual_QL2(na, ns) = begin
        logitα ~ Normal()
        λ ~ Beta()
        logitκ ~ Normal()
        R₀ ~ truncated(Beta(), 0.0, 1/na)
        η ~ truncated(LogNormal(), 0.0, 50.0)
        ϵ ~ Beta() 
        
        return (α=logistic(logitα), λ = λ, η=η,R₀=R₀, ϵ = ϵ, κ = logistic(logitκ), # Hyperparameters
                R = Vector{typeof(λ)}([1.0/na, R₀]), # Initial latent values
                β = Vector{typeof(λ)}([max(0.0, η * (1/na - R₀))]),
                na = na,
                ns = ns,
                Values = Matrix{typeof(logitα)}(fill(1/na, na, ns)))
    end

    mdl = DynamicBeta_Counterfactual_QL2(nactions, nstates)

    @evolution mdl begin
        
        # Update Values
        for i = 1:na
            if i == a
                delta_rule!(Values, s, a, r, α)
            else
                delta_rule!(Values, s, i, 1.0 - r, κ*α)
            end
        end

        # Update R  
        R[1] += λ*(r - R[1])
        R[2] += λ*(R[1] - R[2])

        β[1] = min(500.0, max(0.0, β[1] + η * (R[1] - R[2])))

    end

    @observation mdl begin
        P = epsilon_greedy!(softmax!(β[1] * @views(Values[:,s])), ϵ)
        Categorical(P)

    end
    return mdl
end