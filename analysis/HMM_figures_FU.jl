using CSV, DataFrames, StatsPlots, Statistics
## Hyperparameters second paradigm
flist = readdir("WMM/FU_HMM/Hyperparameters")
df = DataFrame()

for f in flist
    tmp = CSV.read(string("WMM/FU_HMM/Hyperparameters/", f), DataFrame)
    tmp[!, :Subject] .= f
    append!(df, tmp)
end 

volsA = [[],[],[]]
stosA = [[],[],[]]
distA = [[],[],[]]
fulldistA = [[],[],[]]

volsO = [[],[],[]]
stosO = [[],[],[]]
distO = [[],[],[]]
fulldistO = [[],[],[]]

gdf = groupby(df, [:Subject, :Agent, :Task])
for g in gdf
    t = g[in.(g.parameters, [["p[1]", "p[2]", "p[3]"]]),:mean]
    v = sum(t)
    s = g[g.parameters .== "ρ",:mean][1]
    
    if g.Task[1] == "task1" && g.Agent[1] == "act"
        push!(volsA[1], v)
        push!(stosA[1], s)
        push!(distA[1], t ./ sum(t))
        push!(fulldistA[1], vcat(1 - v, t))
    elseif g.Task[1] == "task2" && g.Agent[1] == "act"
        push!(volsA[2], v)
        push!(stosA[2], s)
        push!(distA[2], t ./ sum(t))
        push!(fulldistA[2], vcat(1 - v, t))
    elseif g.Task[1] == "task3" && g.Agent[1] == "act"
        push!(volsA[3], v)
        push!(stosA[3], s)
        push!(distA[3], t ./ sum(t))
        push!(fulldistA[3], vcat(1 - v, t))
    elseif g.Task[1] == "task1" && g.Agent[1] == "obs"
        push!(volsO[1], v)
        push!(stosO[1], s)
        push!(distO[1], t ./ sum(t))
        push!(fulldistO[1], vcat(1 - v, t))
    elseif g.Task[1] == "task2" && g.Agent[1] == "obs"
        push!(volsO[2], v)
        push!(stosO[2], s)
        push!(distO[2], t ./ sum(t))
        push!(fulldistO[2], vcat(1 - v, t))
    elseif g.Task[1] == "task3" && g.Agent[1] == "obs"
        push!(volsO[3], v)
        push!(stosO[3], s)
        push!(distO[3], t ./ sum(t))
        push!(fulldistO[3], vcat(1 - v, t))
    else
        println("that shouldn't happen")
    end
end


## Fig. 1 : Volatility
scatter(volsO[1], volsA[1], color=:blue, alpha=0.4, label="Task 1")
#scatter!([mean(volsO[1])], [mean(volsA[1])], xerror = std(volsO[1])/sqrt(length(volsO[1])), yerror = std(volsA[1])/sqrt(length(volsA[1])),color=:black, markerstrokecolor=:blue, markerstrokewidth=2, label="Task 1")
scatter!(volsO[2], volsA[2], color=:pink, alpha=0.6, label="Task 2")
#scatter!([mean(volsO[2])], [mean(volsA[2])], xerror = std(volsO[2])/sqrt(length(volsO[2])), yerror = std(volsA[2])/sqrt(length(volsA[2])),color=:red, markerstrokecolor=:pink, markerstrokewidth=2, label="Task 2")
scatter!(volsO[3], volsA[3], color=:green, alpha=0.4, label="Task 3")
#scatter!([mean(volsO[3])], [mean(volsA[3])], xerror = std(volsO[3])/sqrt(length(volsO[3])), yerror = std(volsA[3])/sqrt(length(volsA[3])),color=:black, markerstrokecolor=:green, markerstrokewidth=2, label="Task 3")
plot!([0, 0.4], [0, 0.4], label="", color=:black, linestyle=:dash)
xlabel!("Ideal Observer")
ylabel!("Subjects")
title!("Volatility")
xlims!(0, 0.15)


## Fig. 2 : Stochasticity
scatter(stosO[1], stosA[1], label="Task1", color=:blue, alpha=0.4)
#scatter!([mean(stosO[1])], [mean(stosA[1])], xerror = std(stosO[1])/sqrt(101), yerror = std(stosA[1])/sqrt(101),color=:black, markerstrokecolor=:blue, markerstrokewidth=2, label="Task 1")
scatter!(stosO[2], stosA[2], label="Task 2", color=:pink, alpha=0.5)
#scatter!([mean(stosO[2])], [mean(stosA[2])], xerror = std(stosO[2])/sqrt(103), yerror = std(stosA[2])/sqrt(103),color=:red, markerstrokecolor=:pink, markerstrokewidth=2, label="Task 2")
scatter!(stosO[3], stosA[3], label="Task 3", color=:green, alpha=0.5)
plot!([0, 0.4], [0, 0.4], label="", color=:black, linestyle=:dash)
xlabel!("Ideal Observer")
ylabel!("Subjects")
title!("Stochasticity")
xlims!(0.05, 0.3)

## Fig. 3 : Distributions
plot(distO[1], color = :grey, alpha = 0.3, label="")
plot!(mean(distO[1]), color=:black, linewidth=3,label = "Ideal Observer")
plot!(distA[1], color = :blue, alpha = 0.2, label="")
plot!(mean(distA[1]), color=:blue, linewidth=3,label = "Subjects")
xticks!(1:3)
xlabel!("Distance to current set")
ylabel!("Normalized transition probability")
title!("Task 1")
##
plot(distO[2], color = :grey, alpha = 0.3, label="")
plot!(mean(distO[2]), color=:black, linewidth=3,label = "Ideal Observer")
plot!(distA[2], color = :pink, alpha = 0.5, label="")
plot!(mean(distA[2]), color=:red, linewidth=3,label = "Subjects")
xticks!(1:3)
xlabel!("Distance to current set")
ylabel!("Normalized transition probability")
title!("Task 2")

##
plot(distO[3], color = :grey, alpha = 0.3, label="")
plot!(mean(distO[3]), color=:black, linewidth=3,label = "Ideal Observer")
plot!(distA[3], color = :green, alpha = 0.4, label="")
plot!(mean(distA[3]), color=:green, linewidth=3,label = "Subjects")
xticks!(1:3)
xlabel!("Distance to current set")
ylabel!("Normalized transition probability")
title!("Task 3")

## Data
flist = readdir("WMM/FU_HMM/Data")
datadf = DataFrame()

for f in flist
    tmp = CSV.read(string("WMM/FU_HMM/Data/", f), DataFrame)
    tmp[!, :Subject] .= f
    append!(datadf, tmp)
end 

# Recompute pre-swtich perseveration
datadf[!,:correct_or_persev] = copy(datadf.persev)
pres_in_block!(datadf, blockId = :blockNum) 
datadf[datadf.negpresInBlock .> -3,:correct_or_persev] = datadf[datadf.negpresInBlock .> -3,:correct]

# Recode for stable stims 
datadf[datadf.stableAS, :correct_or_persev] = datadf[datadf.stableAS, :correct]

# Block count locked on HMM switch
datadf[!,:switchNum] = zeros(Int, nrow(datadf)) 
gdf = groupby(datadf, [:Subject, :task, :sessNum])
for g in gdf
    idx = diff(vcat(0, findall(g.HMMSwitch .== 1) .-  1, nrow(g)))
    t = 0
    for b in eachindex(idx)
        g[t+1:t+idx[b], :switchNum] .= b
        t += idx[b]
    end 
end
pres_in_block!(datadf, blockId = :switchNum)

## Env 2
df2 = datadf[datadf.task .== "task2",:]
df2.condition[df2.stableAS] .= 3

## Global switches
compdf = df2[(df2.switchType .== 3) .* (0 .< df2.condition .<= 2),:]
grpstats = grp_stats_hmm(compdf)
grp_plot(grpstats, "persev", [1,2,3]; xlims=(-3, 8), ylims=(0, 1), yticks=0:0.2:1, label=["New rule" "Partial switch - changing stim"], xlabel="Stimulus presentations", ylabel="Perseveration")
savefig("WMM/Figures/HMM/expe2_env2_global_complete.pdf")
## Partial switches global VS non global
partdf = df2[df2.condition .== 3,:]
partdf.persev .= partdf.correct_or_persev
partdf.condition[partdf.switchType .< 3] .= 1 # Recode conditions as switch types
partdf.condition[partdf.switchType .== 3] .= 2
grpstats2 = grp_stats_hmm(partdf)
grp_plot(grpstats2, "persev", [1,2,3]; xlims=(-3, 8), ylims=(0, 1), yticks=0:0.2:1, label=["Non global switches" "Global switches"], xlabel="Stimulus presentations", ylabel="Perseveration")

savefig("WMM/Figures/HMM/expe2_env2_globalVSnonglobal_partial.pdf")
## Env 3
df3 = datadf[datadf.task .== "task3",:]
df3.condition[df3.stableAS] .+= 2

## Global switches
compdf = df3[(df3.switchType .== 3) .* (0 .< df3.condition .<= 2),:]
grpstats = grp_stats_hmm(compdf)
grp_plot(grpstats, "persev", [1,2,3]; xlims=(-3, 8), ylims=(0, 1), yticks=0:0.2:1, label=["New rule" "Partial switch - changing stim"], xlabel="Stimulus presentations", ylabel="Perseveration")

## Partial switches global VS non global
partdf = df3[df3.condition .>= 3,:]
partdf.persev .= partdf.correct_or_persev
partdf.condition[partdf.switchType .< 3] .-= 2 # Recode conditions as switch types
#partdf.condition[partdf.switchType .== 3] .= 2
grpstats2 = grp_stats_hmm(partdf)
grp_plot(grpstats2, "persev", [1,2,3]; xlims=(-3, 8), ylims=(0, 1), yticks=0:0.2:1, label=["Non global switches" "Global switches"], xlabel="Stimulus presentations", ylabel="Perseveration")

## Trap sensitivity

sdf = groupby(df2, :subject)
trapmat = zeros(length(sdf), 4, 2) # P(error) at -1, +1, +2, +3, before and after the switch
for si in 1:length(sdf)
    g = sdf[si]
    trapidx = findall(g.trap .== 1)
    perror = zeros(4, 2)
    count = zeros(Int, 4, 2)
    for t in trapidx
        if g.stableAS[t] # Stable association only
            pre_post_switch = (g.presInBlock[t] < abs(g.negpresInBlock[t])) + 1 # 1 : pre-switch, 2 : post switch
            if (pre_post_switch == 1) || ((pre_post_switch == 2) && (g.switchType[t] < 3))  # if post switch : non global switch
                stim = g.stimulus[t]
                # Get the last and 3 next presentations of the same stimulus
                tprev = findlast(g.stimulus[1:t-1] .== stim)
                tnext = findall(g.stimulus[t+1:end] .== stim)[1:min(end, 3)]

                # Compute the error rate before and after the trap
                if !isnothing(tprev)
                    perror[1,pre_post_switch] += 1 - g.correct[tprev]
                    count[1,pre_post_switch] += 1
                end
                for i in eachindex(tnext)
                    perror[i+1,pre_post_switch] += 1 - g.correct[tnext[i]]
                    count[i+1,pre_post_switch] += 1
                end
            end
        end
    end
    @show count
    perror ./= count # averaging
    trapmat[si,:,:] .= perror
end

##
cpal = palette(:Dark2)
boxplot([0 2 4 6], trapmat[:,:,1], color=cpal[1], alpha=0.5, label=["Before switch" "" "" ""])
dotplot!([0 2 4 6], trapmat[:,:,1], color=cpal[1], alpha=0.5, label="")
boxplot!([1 3 5 7], trapmat[:,:,2], color=cpal[2], alpha=0.5, label=["After switch" "" "" ""])
dotplot!([1 3 5 7], trapmat[:,:,2], color=cpal[2], alpha=0.5, label="")
xticks!(0.5:2:6.5, ["-1", "+1", "+2", "+3"])
xlabel!("Presentations from trap")
ylabel!("P(error)")
title!("Trap sensitivity before and after a non global switch for stable associations", titlefontsize=10)

## Latency Env1
df1 = datadf[datadf.task .== "task1",:]
gdf = groupby(df1, [:subject, :blockNum])
latdf = DataFrame(subject = zeros(Int, length(gdf)), blockNum = zeros(Int, length(gdf)), condition = zeros(Int, length(gdf)), switchType = zeros(Int, length(gdf)), latency_trials = zeros(Int, length(gdf)), latency_pres = zeros(length(gdf)))

for gi in 1:length(gdf)
    g = gdf[gi]
    latdf[gi, :subject] = g.subject[1]
    latdf[gi, :blockNum] = g.blockNum[1]
    latdf[gi, :condition] = g.condition[1] < 4 ? g.condition[1] : 3 # Reclassify stable associations in condition 3
    latdf[gi, :switchType] = g.switchType[1]
    if g.blockNum[1] > 1
        idx = findfirst(g.HMMSwitch .== 1)
        latdf[gi, :latency_trials] = idx
        mean([g.stimulus[1:idx] .== stim for stim = 1:3])
        latdf[gi, :latency_pres] = mean([count(g.stimulus[1:idx] .== stim) for stim = 1:3]) # Averagre presentations from the begining of the block to the switch
    end
end

##Recode conditions as condtiion X switchType
latdf = latdf[latdf.condition .> 0, :]
latdf = latdf[latdf.switchType .== 3,:] # keep only global switches

sdf = groupby(latdf, :condition)
latsum = combine(sdf, :latency_trials => mean => :latency_mean, :latency_trials => sem => :latency_sem)

@df latsum bar(:latency_mean, yerror=:latency_sem, group=:condition, bar_width=0.7, label="", palette=:Dark2, legend_position=:topleft)
xticks!([1,2], ["New rule", "Recurent rule"])
ylabel!("Rule change - Internal switch latency (# trials)")
title!("Switch latency for global switches")
ylims!(0, 10)


## Latency Env2
df2 = datadf[datadf.task .== "task2",:]
gdf = groupby(df2, [:subject, :blockNum])
latdf = DataFrame(subject = zeros(Int, length(gdf)), blockNum = zeros(Int, length(gdf)), condition = zeros(Int, length(gdf)), switchType = zeros(Int, length(gdf)), latency_trials = zeros(Int, length(gdf)), latency_pres = zeros(length(gdf)))

for gi in 1:length(gdf)
    g = gdf[gi]
    latdf[gi, :subject] = g.subject[1]
    latdf[gi, :blockNum] = g.blockNum[1]
    latdf[gi, :condition] = g.condition[1] < 4 ? g.condition[1] : 3 # Reclassify stable associations in condition 3
    latdf[gi, :switchType] = g.switchType[1]
    if g.blockNum[1] > 1
        idx = findfirst(g.HMMSwitch .== 1)
        latdf[gi, :latency_trials] = idx
        mean([g.stimulus[1:idx] .== stim for stim = 1:3])
        latdf[gi, :latency_pres] = mean([count(g.stimulus[1:idx] .== stim) for stim = 1:3]) # Averagre presentations from the begining of the block to the switch
    end
end

##Recode conditions as condtiion X switchType
latdf = latdf[latdf.condition .> 0, :]
latdf = latdf[latdf.switchType .== 3,:] # keep only global switches

sdf = groupby(latdf, :condition)
latsum = combine(sdf, :latency_trials => mean => :latency_mean, :latency_trials => sem => :latency_sem)

cpal = palette(:Dark2)[[1,3]]
@df latsum bar(:latency_mean, yerror=:latency_sem, group=:condition, bar_width=0.7, label="", palette=cpal, legend_position=:topleft)
xticks!([1,2], ["New rule", "Partial change"])
ylabel!("Rule change - Internal switch latency (# trials)")
title!("Switch latency for global switches")
ylims!(0,10)