function probe_PH(nactions, nstates)
    @model PROBEmodelPH(na, ns) = begin
        μ₀ ~ Beta()
        λ ~ Beta()
        κ ~ Beta()
        β ~ truncated(LogNormal(), 0.0, 50.0)
        ϵ ~ Beta(1, 10) 
        Θ ~ MvNormal(3, 1.0)
        # τ ~ Beta(10, 1)
        # θ ~ truncated(Beta(), 0.0, 0.5)#truncated(Normal(), -2, 3)
        η ~ Beta(1, 1)
        #m ~ Categorical(5)
        τ = logistic(Θ[1])
        θ = 0.5 * logistic(Θ[2])
        η = logistic(Θ[3])
        return generate_latent(Val(:PROBEmodelPH), na, ns, μ₀ = μ₀,λ=λ, κ = κ, β = β, ϵ = ϵ, τ = τ, η = η, θ = θ, m = 1) 
    end

    mdl = PROBEmodelPH(nactions, nstates)
    @evolution mdl begin
        # Update reliability 
        nmonitored = findall(Reliabilities .> 0.0)
        for k in nmonitored
            Reliabilities[k] = τ * Reliabilities[k] + (1 - τ) * (1 - Reliabilities[k]) / (length(nmonitored) - 1) # Prior
            lh = Probas[a,s,k] / sum(Probas[:,s,k]) # Likelihood
            Reliabilities[k] *= r == 1 ? lh : 1 - lh
        end
        Reliabilities ./= sum(Reliabilities)

        # Learn values and update probas
        α = μ[1]
        δ = abs(r - Values[a,s,Actor[1]]) 
        for i = 1:na
            if i == a
                Values[i,s,Actor[1]] += α * (r - Values[i,s,Actor[1]])
            else
                Values[i,s,Actor[1]] += κ * α * (1 - r - Values[i,s,Actor[1]])
            end
        end

        Probas[:,s,Actor[1]] .+= (1 - r)/(na-1)
        Probas[a,s,Actor[1]] += r*na/(na-1) - 1/(na-1)

        # Update learning rate
        μ[1] *= 1 - λ
        μ[1] += λ * abs(δ)
        
        # Set selection
        set_selection!(Reliabilities, Values, Probas, Actor, Order, Probing, LTM, Pointer, na, θ, η, m)
    end

    @observation mdl begin
        Qactions .= @views(Values[:,s,Actor[1]])
        Qactions .*= β
        return Categorical(epsilon_greedy!(softmax!(Qactions), ϵ))
    end

    return mdl
end

function generate_latent(::Val{:PROBEmodelPH}, na, ns ; μ₀ = μ₀,λ=λ, κ, β, ϵ, τ, θ, η, m)
    initF = ( (τ - 1.0) * 0.3 + τ * 0.7 ) / (2 * τ - 1.0)
    return (μ₀ = μ₀,λ=λ, κ = κ, β = β, ϵ = ϵ, τ = τ, θ = θ, η = η,  m = m,
        na = na,
        ns = ns,
        Values = Array{typeof(λ),3}(fill(1/na, na, ns, m+2)),
        Probas = Array{typeof(λ),3}(fill(1/na, na, ns, m+2)),
        LTM = Array{typeof(λ),3}(zeros(na, ns, 2)),
        Reliabilities = Vector{typeof(λ)}(vcat(initF, zeros(m), 1 - initF)),
        Actor = [1],
        Order = vcat(1, zeros(Int, m-1)), 
        Pointer = [1], # LTM pointer
        Probing = [false], # Probing state
        μ = [μ₀],
        Qactions = Vector{typeof(β)}(fill(1/na, na)))
end