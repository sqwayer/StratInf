function set_selection!(Reliabilities, Values, Probas, Actor, Order, Probing, LTM, Pointer, na, θ, η, m)

    bestF, bestTS = findmax(Reliabilities[1:end-1]) 
    if bestF > 0.5 && bestTS < m+1 # One TS is selected and it's not the Probe
        Actor[1] = bestTS
    elseif bestF > 0.5 && bestTS == m+1 # The winning TS is the Probe
        # Find a slot
        Actor[1] = argmin(Order)

        if Reliabilities[Actor[1]] > 0 # The slot is occupied -> Transfer the set to the LTM
            LTM[:,:,1] += Values[:,:,Actor[1]]
            LTM[:,:,2] += Probas[:,:,Actor[1]] ./ sum(Probas[:,:,Actor[1]], dims=1)
            Pointer[1] += 1
        end

        # Transfer the probe
        Reliabilities[Actor[1]] = Reliabilities[m+1]
        Values[:,:,Actor[1]] = Values[:,:,m+1]
        Probas[:,:,Actor[1]] = Probas[:,:,m+1]

    else # There is no winning TS --> Probe
        probe_creation!(Reliabilities, Values, Probas, LTM, Pointer, na, θ, η, m)
        Probing[1] = true
        Actor[1] = m+1
        return
    end

    # Remove the probe (always if there is a winning set)
    Probing[1] = false
    Reliabilities[m+1] = 0.0
    Reliabilities ./= sum(Reliabilities)

    # Update order (always if there is a winning set)
    Order[Actor[1]] = maximum(Order) + 1
end

function probe_creation!(Reliabilities, Values, Probas, LTM, Pointer, na, θ, η, m)

    # Create a probe from the LTM
    Values[:,:,m+1] .= LTM[:,:,1]
    Probas[:,:,m+1] .= LTM[:,:,2]

    # Add the currently monitored sets
    monitored = findall(Reliabilities .> 0.0)
    for k in monitored
        Values[:,:,m+1] .+= Values[:,:,k]
        Probas[:,:,m+1] .+= Probas[:,:,k] ./ sum(Probas[:,:,k], dims=1)
    end

    # Normalize
    Values[:,:,m+1] ./= Pointer[1] + length(monitored)
    Probas[:,:,m+1] ./= Pointer[1] + length(monitored)

    # Add noise
    Values[:,:,m+1] .*= η
    Values[:,:,m+1] .+= (1 - η)/na
    Probas[:,:,m+1] .*= η
    Probas[:,:,m+1] .+= (1 - η)/na

    # Initialize the reliability
    # H = 0.0 # Buffer Entropy
    # for k in monitored
    #     H += log(max(1e-20, Reliabilities[k])) * Reliabilities[k]
    # end
    # Fopt = 1/(1+exp(-H))
    # Fopt = min(1, max(0, θ*0.5 + (1-θ)*Fopt)) # Black magic that breaks the sampler

    Reliabilities .*= (1-θ)#(1 - Fopt)
    Reliabilities[m+1] = θ# Fopt

    return
end
