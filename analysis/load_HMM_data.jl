function load_HMM_data(folder)
    flist = filter(x -> occursin(".csv", x), readdir(folder))
    datadf = DataFrame()

    for f in flist
        tmp = CSV.read(string(folder, "/", f), DataFrame)
        tmp[!, :Subject] .= f
        tmp.sessNum .= parse(Int, f[7])
        tmp[!,:zrt] = (tmp.rt .- mean(tmp.rt)) ./ std(tmp.rt)
        
        # Add model predictions
        for mdl in ["cfql", "full_adaptive_cfql2", "SI_MultVol_SampleAction", "probe_LRRel"]
            mdlDf = CSV.read(string("/Users/sami/PhD/Model_Tasks_Data/Data/WMM/Fits/Raw2/all_fits/Fits/", mdl, "_", tmp.task[1], "/choicesProba/", f), DataFrame)
            
            tmp[!, Symbol("choiceProb_$mdl")] = mdlDf[:,Symbol("choiceProb_$mdl")]
            tmp[!, Symbol("persevProb_$mdl")] = mdlDf[:,Symbol("persevProb_$mdl")]
            tmp[!, Symbol("correctOrPersevProb_$mdl")] = mdlDf[:,Symbol("correctOrPersevProb_$mdl")]
            tmp[!, Symbol("H_$mdl")] = mdlDf[:,Symbol("entropy_$mdl")]
            tmp[!, Symbol("L_$mdl")] = mdlDf[:,Symbol("latent_$mdl")]

        end
        
        append!(datadf, tmp)
        
    end 
return datadf
end

function process_HMM_data!(datadf)

    # Recompute pre-swtich perseveration
    datadf[!,:correct_or_persev] = copy(datadf.persev)
    pres_in_block!(datadf, blockId = :blockNum) 
    datadf[datadf.negpresInBlock .>= -3,:correct_or_persev] = datadf[datadf.negpresInBlock .>= -3,:correct]

    # Recode for stable stims 
    datadf[datadf.stableAS, :correct_or_persev] = datadf[datadf.stableAS, :correct]

    # Block count locked on HMM switch
    datadf[!,:switchNum] = zeros(Int, nrow(datadf)) 
    gdf = groupby(datadf, [:Subject, :task, :sessNum])
    for g in gdf
        idx = diff(vcat(0, findall(g.HMMSwitch .== 1) .- 1, nrow(g)))

        t = 0
        for b in eachindex(idx)
            g[t+1:t+idx[b], :switchNum] .= b
            t += idx[b]
        end 
    end
    pres_in_block!(datadf, blockId = :switchNum)
    trials_in_block!(datadf, blockId = :switchNum)
    gdf = groupby(datadf, [:subject, :sessNum, :blockNum])
    datadf[!,:nextCondition] .= 0
    for gi = 1:length(gdf)
        if gdf[gi].blockNum[1] < 39 
            gdf[gi].nextCondition .= gdf[gi+1].condition[1]
        end
    end
    # For task 1
    df1 = datadf[datadf.task .== "WMM1",:]
    df1.condition[df1.condition .== 4] .= 1 # Merge conditions 1 and 4

    for t = 1:nrow(df1)
        if df1[t, Symbol("isStable_$(df1.stimulus[t])")] && df1.condition[t] > 0 
            df1.condition[t] = 4 # Make stable stims a special condition
        end
    end
    df1[(df1.condition .== 4) .* (df1.presInBlock .<= 10), :persev] .= df1[(df1.condition .== 4) .* (df1.presInBlock .<= 10), :correct] # persev = correct for the first presentations in condition 4 (stable associations). Otherwise perseveration would be = 0

    return df1

end