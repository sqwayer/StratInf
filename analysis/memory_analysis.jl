## Supplementary analysis to add a test on the effect of recurrence
# Look at the performance in a block depending on its similarity with the recurrent rule 

function get_recurrent_rules(df)
    rules = []

    gdf = groupby(df[df.condition .== 2,:], :blockNum)
    for block in gdf
        rule = zeros(Int, 3)
        for stim = 1:3
            idx = findfirst(block.stimulus .== stim)
            rule[stim] = block[idx, :corResp]
        end

        if !in(rule, rules)
            push!(rules, rule)
        end
    end
    return rules
end

function fb_as_recurrent(stim, choice, fb, recRule)
    N = length(stim)
    concordant = 0
    discordant = 0
    for i in eachindex(stim)
        if choice[i] == recRule[stim[i]]
            concordant += fb[i]
            discordant += !fb[i]
        else
            discordant += fb[i]
        end
    end
    return concordant, discordant

    # if concordant == 0 && discordant == 0
    #     return 0, 0
    # end

    # ncat = 7
    # cutoffs = 1/ncat:1/ncat:1
    # concordCat = findfirst(concordant / (concordant + discordant) .<= cutoffs)
    # discordCat = findfirst(discordant / (concordant + discordant) .<= cutoffs)
    # return cutoffs[concordCat], cutoffs[discordCat]
end

function distance_from_recurrent!(df, nrec=1)

    df[!, :choiceAsRec] .= 0
    df[!, :correctAsRec] .= 0

    for k in 1:nrec
        df[!, Symbol("distToRec_$k")] .= 3
        df[!, Symbol("concordEvRec_$k")] .= NaN
        df[!, Symbol("discordEvRec_$k")] .= NaN
    end

    sessdf = groupby(df, :sessNum)
    for sess in sessdf
        recRules = get_recurrent_rules(sess)
        blockdf = groupby(sess, :blockNum)
        for block in blockdf
            dist = fill(3, nrec)
            for stim = 1:3
                # Find the current rule for the given stim and distance to recurrent rule
                idx = findfirst(block.stimulus .== stim)
                if !isnothing(idx)
                    for k in eachindex(recRules)
                        if block[idx, :corResp] == recRules[k][stim]
                            dist[k] -= 1
                        end
                    end
                    

                    # Code the choice compared to the recurrent rules 
                    ind = findall(block.stimulus .== stim)
                    for k in eachindex(recRules)
                        block[ind, :choiceAsRec] .+= k .* (block[ind, :choice] .== recRules[k][stim])
                        block[ind, :correctAsRec] .+= k .* (block[ind, :corResp] .== recRules[k][stim])
                    end

                    # Compute the number/frequency of anterior trials with concordant/discordant evidence with the recurrent rule
                    if idx > 1
                        stimVec = block[1:idx-1, :stimulus]
                        choiceVec = block[1:idx-1, :choice]
                        fbVec = block[1:idx-1, :fb]

                        for k in eachindex(recRules)
                            concordant, discordant = fb_as_recurrent(stimVec, choiceVec, fbVec, recRules[k])
                            block[idx, Symbol("concordEvRec_$k")] = concordant 
                            block[idx, Symbol("discordEvRec_$k")] = discordant 
                        end

                    end
                end
            end
            for k in 1:nrec
                block[:,Symbol("distToRec_$k")] .= dist[k]
            end
        end
    end
end

function evidence_effect!(df)
    df[!,:congruentChoice] .= NaN
    df[!,:congruentDiff] .= NaN
    df[!,:evForRec] .= NaN
    df[!, :incongruentDiff] .= NaN
    df[!, :consistent] .= NaN
    blockdf = groupby(df, [:sessNum, :blockNum])
    for block in blockdf

        for consist = 0:1
            if consist == 1
                # Find all choices consistent with a recurrent rule
                fbIdx = findall(block.choiceAsRec .> 0)
            elseif consist == 0
                # Find all choices inconsistent with a recurrent rule
                fbIdx = findall(block.choiceAsRec .== 0)
            end
        
            for idx in fbIdx
                
                evForRec = block.choiceAsRec[idx] # choice was evidence for rec rule 
                
                posEv = block.fb[idx]

                evStim = block.stimulus[idx]

                congruentChoice = 1
                prevCongruentChoice = 1
                distFromEv = Inf

                for stim = 1:3
                    if stim ≠ evStim
                        stimIdx = findfirst(block.stimulus[idx+1:end] .== stim) # Find the next presentation
                    
                        absIdx = block[idx,:trialNum]
                        prevIdx = findlast(df[1:absIdx-1, :stimulus] .== stim)
                        if !isnothing(stimIdx) && !isnothing(prevIdx)
                            if stimIdx < distFromEv 
                                distFromEv = stimIdx 
                            end
                            stimIdx += idx

                            if consist == 1
                                # Congruent => consistent with the initial recurrent rule
                                congruentChoice *= block.choiceAsRec[stimIdx] == evForRec

                                prevCongruentChoice *= df.choiceAsRec[prevIdx] == evForRec
                            elseif consist == 0
                                # Congruent => consistent with any recurrent rule
                                congruentChoice *= block.choiceAsRec[stimIdx] ≠ 0

                                prevCongruentChoice *= df.choiceAsRec[prevIdx] ≠ 0
                            end
                        end

                        
                    end
                end
                block.congruentChoice[idx] = congruentChoice
                block.congruentDiff[idx] = congruentChoice - prevCongruentChoice
                block.evForRec[idx] = posEv
                block.incongruentDiff[idx] = prevCongruentChoice
                block.consistent[idx] = consist
            end
        end
    end
end

function group_level_mem!(groupdf, nrec=1)
    gdf = groupby(groupdf, :subject)
    for ddf in gdf
        ddf[!,:trialNum] = 1:nrow(ddf)
        distance_from_recurrent!(ddf, nrec)
        evidence_effect!(ddf)
    end
end

