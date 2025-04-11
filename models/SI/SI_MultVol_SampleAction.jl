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

const typeIdx = hcat([findall([sum(any(all_sets[i] .≠ all_sets[j], dims=1)) == 1 for j = 1:27]) for i = 1:27], [findall([sum(any(all_sets[i] .≠ all_sets[j], dims=1)) == 2 for j = 1:27]) for i = 1:27], [findall([sum(any(all_sets[i] .≠ all_sets[j], dims=1)) == 3 for j = 1:27]) for i = 1:27])


function SI_MultVol_SampleAction(nactions, nstates)

    @model StrategicInferenceMultVolSA(na, ns) = begin
        ρ ~ truncated(Beta(), 0.5, 1)
        ω₀ ~ Beta() # Long-term memory
        ω ~ MvNormal(5, 1.0) # Volatility for each type of switches + random + stability
        ωᵣ ~ Beta() # Proba switching out of random
        ϵ ~ Beta()
        β ~ LogNormal()

        Nsets = length(all_sets)

        return (ρ = ρ, ω₀ = ω₀,ω = softmax(ω),ωᵣ= ωᵣ,ϵ=ϵ, β=β, na = na, ns = ns,
            all_sets = all_sets, typeIdx = typeIdx,
            Nsets = Nsets,
            Γ = Vector{typeof(ω₀)}(fill(1/(Nsets+1), Nsets+1)),
            tmp_gamma = Vector{typeof(ω₀)}(fill(1/(Nsets+1), Nsets+1)),
            Π = Vector{typeof(ω₀)}(fill(1/(Nsets+1), Nsets+1)),
            Φ = Vector{typeof(ρ)}(fill(1/(Nsets+1), Nsets+1)),
            C = [1.0],
            ϕ₀ = 1/(Nsets+1),
            Q = Matrix{typeof(ρ)}(fill(1/na, na, ns)),
            P = zeros(typeof(ϵ), na), 
            )
    end

    mdl = StrategicInferenceMultVolSA(nactions, nstates)

    @evolution mdl begin
        # Update long-term prior
        C[1] += 1
        α = 1 / C[1]
        Γ .*= 1 - α
        #Γ .+= α .* (ω₀ * ϕ₀ .+ (1-ω₀) .* Φ ./ sum(Φ))
        tmp_gamma .= Φ
        tmp_gamma .*= 1-ω₀
        tmp_gamma .+= ω₀ * ϕ₀
        tmp_gamma .*= α
        Γ .+= tmp_gamma
        
        # Update prior
        Π .= Φ
        Π[1:end-1] .*= ω[end]
        Π[end] *= 1 - ωᵣ
        for i = 1:Nsets
            for st = 1:3 # Switch type 
                for j in typeIdx[i,st]
                    Z = 0.0
                    for k in typeIdx[j,st]
                        Z += Γ[k]
                    end
                    Π[i] += ω[st] * Φ[j] * Γ[i] / Z
                end
            end
            Π[i] += ωᵣ * Φ[end] * Γ[i] / (1 - Γ[end]) # From random set 
            Π[end] += ω[end-1] * Φ[i] * Γ[end] / (1 - Γ[i]) # From set i to random
            # Likelihood
            Π[i] *= all_sets[i][a,s] == r ? ρ : 1-ρ  
        end
        Π[end] *= r == 1 ? 1/na : 1-1/na # likelihood for random set
        # Normalize posterior
        Φ .= Π
        Z = sum(Φ)
        Φ ./= Z
        Q .= 0.0
        for i in 1:Nsets
            Q .+= Φ[i] .* all_sets[i]
        end
        Q .+= Φ[end]/na
    end

    @observation mdl begin

        P .= softmax!(β .* @views(Q[:,s]))
        if !isprobvec(P)
            @show P
            @show Φ
        end
        return Categorical(epsilon_greedy!(P, ϵ))
        

    end

    return mdl
end