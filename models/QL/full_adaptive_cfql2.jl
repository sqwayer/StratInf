function full_adaptive_cfql2(nactions, nstates)
    @model Dynamic_Counterfactual_QL2(na, ns) = begin
        μ₀ ~ Beta()
        λₐ ~ Beta()
        λᵦ ~ Beta()
        logitκ ~ Normal()
        R₀ ~ truncated(Beta(), 0.0, 1/na)
        η ~ truncated(LogNormal(), 0.0, 50.0)
        ϵ ~ Beta() 
        
        return (μ₀ = μ₀, λₐ = λₐ, λᵦ = λᵦ, η=η,R₀=R₀, ϵ = ϵ, κ = logistic(logitκ), # Hyperparameters
                R = Vector{typeof(λₐ)}([1.0/na, R₀]), # Initial latent values
                β = Vector{typeof(λₐ)}([max(0.0, η * (1/na - R₀))]),
                μ = [μ₀],
                na = na,
                ns = ns,
                Values = Matrix{typeof(μ₀)}(fill(1/na, na, ns)))
    end

    mdl = Dynamic_Counterfactual_QL2(nactions, nstates)

    @evolution mdl begin
        α = μ[1]
        δ = abs(r - Values[a,s]) 
        # Update Values
        for i = 1:na
            if i == a
                delta_rule!(Values, s, a, r, α)
            else
                delta_rule!(Values, s, i, 1.0 - r, κ*α)
            end
        end
        μ[1] *= 1 - λₐ
        μ[1] += λₐ * abs(δ)
        # Update R  
        R[1] += λᵦ*(r - R[1])
        R[2] += λᵦ*(R[1] - R[2])

        β[1] = min(500.0, max(0.0, β[1] + η * (R[1] - R[2])))

    end

    @observation mdl begin
        P = epsilon_greedy!(softmax!(β[1] * @views(Values[:,s])), ϵ)
        Categorical(P)

    end
    return mdl
end