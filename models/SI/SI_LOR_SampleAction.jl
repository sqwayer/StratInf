"""
Strategic inference model, with reliabilities defined as log odds ratios (compared to the probability of a random strategy)
"""

const all_sets = [
    [1 1 1; 0 0 0; 0 0 0],
    [1 1 0; 0 0 1; 0 0 0],
    [1 1 0; 0 0 0; 0 0 1],
    [1 0 1; 0 1 0; 0 0 0],
    [1 0 0; 0 1 1; 0 0 0],
    [1 0 0; 0 1 0; 0 0 1],
    [1 0 1; 0 0 0; 0 1 0],
    [1 0 0; 0 0 1; 0 1 0],
    [1 0 0; 0 0 0; 0 1 1],
    [0 1 1; 1 0 0; 0 0 0],
    [0 1 0; 1 0 1; 0 0 0],
    [0 1 0; 1 0 0; 0 0 1],
    [0 0 1; 1 1 0; 0 0 0],
    [0 0 0; 1 1 1; 0 0 0],
    [0 0 0; 1 1 0; 0 0 1],
    [0 0 1; 1 0 0; 0 1 0],
    [0 0 0; 1 0 1; 0 1 0],
    [0 0 0; 1 0 0; 0 1 1],
    [0 1 1; 0 0 0; 1 0 0],
    [0 1 0; 0 0 1; 1 0 0],
    [0 1 0; 0 0 0; 1 0 1],
    [0 0 1; 0 1 0; 1 0 0],
    [0 0 0; 0 1 1; 1 0 0],
    [0 0 0; 0 1 0; 1 0 1],
    [0 0 1; 0 0 0; 1 1 0],
    [0 0 0; 0 0 1; 1 1 0],
    [0 0 0; 0 0 0; 1 1 1]
]

const pmin = 1e-40

function SI_LOR_SampleAction(nactions, nstates)

    @model StrategicInferenceLORSA(na, ns) = begin
        ρ ~ truncated(Beta(), 0.5, 1.0 - 1e-8)
        ω ~ truncated(Beta(), 1e-8, 1.0 - 1e-8)
        ω_a ~ Beta()
        ϵ ~ Beta()
        β ~ LogNormal()

        Nsets = length(all_sets)
        ϕ₀ = 1/(Nsets+1) # All sets + the random one

        return (ρ = ρ, ω = ω, ω_a = ω_a,ϵ=ϵ, β=β, na = na, ns = ns,
            all_sets = all_sets,
            Nsets = Nsets,
            ϕ₀ = ϕ₀,
            Φ = Vector{typeof(ω_a)}(zeros(Nsets)), # Initial reliabilities
            Φ_a = Vector{typeof(ω_a)}(zeros(Nsets)), # Initial long term priors
            C = [1],
            Q = Matrix{typeof(ω_a)}(fill(1/na, na, ns)),
            P = zeros(typeof(ϵ), na),
            S = Vector{typeof(ω_a)}(zeros(Nsets)) # pre-allocate
            )
    end

    mdl = StrategicInferenceLORSA(nactions, nstates)

    @evolution mdl begin
        C[1] += 1
        α = 1 / C[1]

        # Update long term priors and reliabilities of all sets
        ϕ_ar = -logsumexp(Φ_a) # Long term prior for random strat
        ϕ_r = -logsumexp(Φ) # Reliability of random strat
        for k = 1:Nsets
            # Transform LORs into probabilities
            eϕ_a = exp(Φ_a[k]) * logistic(ϕ_ar)
            eϕ  = exp(Φ[k]) * logistic(ϕ_r)

            # Update longterm prior
            Φ_a[k] = log( (1 - α) * eϕ_a + α * (ω_a * eϕ + (1 - ω_a) * ϕ₀)) - log((1 - α) * logistic(ϕ_ar) + α * (ω_a * eϕ_a + (1 - ω_a) * ϕ₀))

            # Update reliability
            Φ[k] = log(ω * eϕ + (1-ω) * eϕ_a) - log(ω * logistic(ϕ_r) + (1-ω) * logistic(ϕ_ar))

            Φ[k] += all_sets[k][a,s] == r ? log(ρ) : log(1-ρ)
            Φ[k] -= r == 1 ? log(1/na) : log(1-1/na)
            
        end
        ϕ_r = -logsumexp(Φ)
        
        S .= exp.(Φ)
        S .*= logistic(ϕ_r)

        Q .= 0.0
        for k = 1:Nsets
            Q .+= S[k] .* all_sets[k]
        end
        Q .+= logistic(ϕ_r) * fill(1/na, na, ns) 
        
    end

    @observation mdl begin
        P .= softmax!(β .* @views(Q[:,s]))
        return Categorical(epsilon_greedy!(P, ϵ))
    end

    return mdl
end
