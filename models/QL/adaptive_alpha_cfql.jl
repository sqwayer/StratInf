""" Counterfactual Q-learning with adaptive learning rate following the Pearce-Hall model (P & H 1982) """

function adaptive_alpha_cfql(nactions, nstates)
    @model PearceHall_Counterfactual_QL(na, ns) = begin
        μ₀ ~ Beta()
        λ ~ Beta()
        logitκ ~ Normal()
        β ~ truncated(LogNormal(), 0.0, 50.0)
        ϵ ~ Beta(1, 10) 
        
        return (μ₀=μ₀, λ = λ, β = β, ϵ = ϵ, κ = logistic(logitκ),
                na = na,
                ns = ns,
                Values = Matrix{typeof(μ₀)}(fill(1/na, na, ns)),
                μ = [μ₀])
    end

    mdl = PearceHall_Counterfactual_QL(nactions, nstates)

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
        # Update V
        μ[1] *= 1 - λ
        μ[1] += λ * abs(δ)

    end

    @observation mdl begin
        Categorical(epsilon_greedy!(softmax!(β * @views(Values[:,s])), ϵ))
    end
    return mdl
end