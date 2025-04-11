function preprocess(df)
    # with pilot data :
    # intVars = [:blockNum, :condition_1, :stimulus, :choice, :corResp, :persevResp]
    # floatVars = [:rt]
    # boolVars = [:stableAS, :newBlock, :fb, :correct, :deterBlock, :persev, :explo, :trap]
    #
    # trialsIdx = findall((df.part .== "stim") .*  (df.training .== "false") .* (df.answered .== "true")) # only with pilot data
   trialsIdx = findall(df.answered)
    # df = df[trialsIdx, [:subject, :task, :sessNum, :blockNum, :condition_1, :deterBlock, :stableAS, :newBlock, :stimulus, :choice, :fb, :correct, :persev, :explo, :trap, :corResp, :persevResp, :rt]] # with pilot data
    df = df[trialsIdx, [:subject, :task, :sessNum, :blockNum, :trialsInBlock, :condition, :stableAS, :isStable_1, :isStable_2, :isStable_3, :newBlock, :stimulus, :choice, :fb, :correct, :persev, :explo, :trap, :corResp, :persevResp, :rt]]

    # only for pilot data
    # df[findall(df[:,:stableAS] .== "\""),:stableAS] .= "false"
    # df[!,intVars] = parse.(Int, df[!,intVars])
    # df[!,floatVars] = parse.(Float64, df[!,floatVars])
    # df[!,boolVars] = parse.(Bool, df[!,boolVars])
    # rename!(df, :condition_1 => :condition) # With pilot data
    df.stableAS = convert.(Bool, df.stableAS)
    df.isStable_1 = convert.(Bool, df.isStable_1)
    df.isStable_2 = convert.(Bool, df.isStable_2)
    df.isStable_3 = convert.(Bool, df.isStable_3)
    df.rt = convert.(Float64, df.rt)
    # Recode everything from 1
    # df.condition = df.condition .+ 1 .+ df.deterBlock * 2 # 0 -> 1 if proba, 2 if deter, 1 -> 3 if proba, 4 if deter / only for pilot data
    df.choice .+= 1
	df.corResp .+= 1
    df.persevResp .+= 1
    df.stimulus .+= 1

    # Correct some mislabeling of stable stims in task 1
    # if df.task[1] == "WMM1"
    #     df = correct_stableAS(df)
    # end

    # Add mutual information of successive choices
   # df.MI = fill(0.0, nrow(df))
    # for s = 1:3
    #     sIdx = findall(df.stimulus .== s)
    #     df.MI[sIdx[2:end]] .= mutual_info(df.choice[sIdx], 2)
    # end

    # Add a correct or persev variable and nexCondition
    df[!,:nextCondition] .= 0
    perblock = groupby(df, :blockNum)
    cop = []
    for gi = 1:length(perblock)
        gpb = perblock[gi]
        append!(cop, correctOrPersev(gpb))
        if gi < length(perblock)
            gpb.nextCondition .= perblock[gi+1].condition[1]
        end
    end
    df[:,:correctOrPersev] = cop

    # Trials/Presentations counts : from the reversal (positive) to the reversal (negative) 
    trials_in_block!(df)
    pres_in_block!(df)

    return df
end

function correct_stableAS(df)
    """ To correct some mislabeling of stable stims """
    stableIdx = findall(df.stableAS)
    for i in stableIdx
        if df[i+10, :condition] != 2
            df[i,:stableAS] = false
        end
    end
    return df
end

function correctOrPersev(perblock)
    """ Adds a variable to count correct at the end of an episode and persev the rest of the time (except for stable stims) """
    stableAS = perblock.stableAS
    correct = perblock.correct
    persev = perblock.persev

    cop = copy(correct)
	
    if perblock.task[1] ≠ "WMM2" || perblock.condition[1] ≠ 4
        cop[1:end-10] = stableAS[1:end-10] .* correct[1:end-10] .+ (1 .- stableAS[1:end-10]) .* persev[1:end-10]
    end
	
    return cop
end

function recompute_latency(df)
    # Column newSwitch = boolean indicating a switch
    df[!,:newSwitch] = vcat(true, diff(df.blockId) .≠ 0)

    # Column T0latency_1 = latency (in presentations) between the actual change and the 1st transition
    # Column T0latency_2 = latency (in presentations) between the 1st transition and the 2nd one  
    # Column T0latency_3 = latency (in presentations) between the 1st transition and the 3rd one, if any   
    df[!,:T0latency_1] = fill(Inf, nrow(df)) 
    df[!,:T0latency_2] = fill(Inf, nrow(df)) 
    df[!,:T0latency_3] = fill(Inf, nrow(df)) 

    # Column T0latency_trials = latency between the change and the 1st swith in trials
    df[!,:T0latency_trials] = fill(Inf, nrow(df)) 

    perblock = groupby(df, [:sessNum, :blockNum]) # Spilt by episodes
    
    for pbi in 2:length(perblock)
        pb = perblock[pbi]
        T0s = [pb[1,Symbol("T0_$s")] for s = 1:3] # Vector of transition latencies
        # - 1st transition = latency between newBlock and change of blockId, in presentations
        epChange = DataFrames.row(pb[1,:]) # First trial of the episode
        switchIdx = findfirst(df.newSwitch[epChange:end])
        if !isnothing(switchIdx) && !all(isinf.(T0s)) # If a switch exists
            switchIdx += epChange - 1 # First transition index in the whole session referential
            stim1 = df.stimulus[switchIdx]
            pb.T0latency_1 .= sum(df.stimulus[epChange:switchIdx] .== stim1) - 1 # #of presentations of the stim in switchIdx during the interval - 1 (1 presentation = 0 additional presentation from the episode change)
            pb.T0latency_trials .= switchIdx - epChange
            
            # - 2nd transition = smallest non zero T0_$s
            if any(T0s .> 0)
                lat = minimum(filter(!isequal(0), T0s))
                stim2 = findfirst(isequal(lat), T0s)
                pb.T0latency_2 .= lat

            # - 3rd transition = highest non zero T0_$s that is not stim2
                lat = maximum(filter(!isequal(0), T0s))
                stim3 = findfirst(isequal(lat), T0s)
                if stim3 ≠ stim2
                    pb.T0latency_3 .= lat
                end
            end
        end
    end
    return df
end
function trials_in_block!(df; blockId = :blockNum)
    df[!,:trialsInBlock] = zeros(Int, nrow(df))
    df[!,:negtrialsInBlock] = zeros(Int, nrow(df))
    gdf = groupby(df, [blockId, :subject, :task, :sessNum])
    for g in gdf
        g.trialsInBlock .= 1:nrow(g)
        g.negtrialsInBlock .= -nrow(g):-1
    end
end

function pres_in_block!(df; blockId = :blockNum, condition = :condition)
    df[!,:presInBlock] = zeros(Int, nrow(df))
    df[!,:negpresInBlock] = zeros(Int, nrow(df))
    #df[!,:nextCondition] = zeros(Int, nrow(df))
    gdf = groupby(df, [blockId, :subject, :task, :sessNum, :stimulus])
    # for gi in 1:length(gdf)
    #     g = gdf[gi]
    #     g.presInBlock .= 1:nrow(g)
    #     g.negpresInBlock .= -nrow(g):-1
    #     # if gi < length(gdf)
    #     #     g.nextCondition .= gdf[gi+1][1,condition]
    #     # end
    # end
    for g in gdf
        for i in 1:nrow(g)
            g.presInBlock[i] = i
            g.negpresInBlock[i] = i-nrow(g)-1
        end
    end
end

function nextCondition!(df)
    maxblocknum = maximum(df.blockNum)
    gdf = groupby(df, [:subject, :sessNum, :blockNum])
    df[!,:nextCondition] .= 0
    for gi = 1:length(gdf)
        if gdf[gi].blockNum[1] < maxblocknum 
            for i in 1:nrow(gdf[gi])
                gdf[gi].nextCondition[i] = gdf[gi+1].condition[1]
            end
        end
    end
end