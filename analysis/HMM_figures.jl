using CSV, DataFrames, StatsPlots, Statistics

## Data
folder = "data/HMM_data"
datadf = load_HMM_data(folder)
df1 = process_HMM_data!(datadf)


## Show RTs locked on strategic changes
alldf = copy(df1)
grpstats = grp_stats_hmm(alldf)
grpstats = grpstats[0 .< grpstats.condition,: ]

# Mean and sem of Z-RT for correct and incorrect responses
summaryRT_sub = combine(groupby(df1, [:subject, :correct]), :zrt => mean => :zrt)
summaryRT_group = combine(groupby(summaryRT_sub, :correct), :zrt => mean, :zrt => sem)

grp_plot_hmm(grpstats, "rt", [1,2,3,4]; xlims=(-2, 9), xticks=-2:2:9, label="", ylims = (-0.5, 0.5), xlabel="Stimulus presentations", ylabel="RT (z-scored)", legend_position = :topright, linestyle=[:solid :solid :solid :solid], background_color=:transparent, foreground_color=:black, labelfontsize=20, tickfontsize=14, dpi=300)


plot!([-2, 9], repeat(summaryRT_group.zrt_mean', 2, 1), ribbon=repeat(summaryRT_group.zrt_sem', 2, 1), linewidth=3, linestyle=:dash, label="", color=:grey)

## Exploration bouts

df1[!,:inExplo] .= zeros(Int, nrow(df1))
df1[!,:beforeExplo] .= fill(-1e5, nrow(df1))
df1[!,:afterExplo] .= fill(1e5, nrow(df1))
df1[!,:outExplo] .= fill(-1e5, nrow(df1))
before_after = 5
maxInExplo = 25 
gdf = groupby(df1, [:subject, :sessNum])
for gd in gdf
    exploIn = findall((gd.HMMSwitch .== 1) .* (gd.switchType .== -3)) # Trials indices of explo in
    exploOut = findall((gd.ExploOut .== 1) .* (gd.switchType .== -3)) # Trials indices of explo out
    if length(exploIn) ≠ length(exploOut)
        error()
    end
    for t in eachindex(exploIn)
        tin = exploIn[t]
        tout = exploOut[t]
        for s = 1:3 # Find all presentations of each stim
            # Find all previous presentation before explo in and all further presentation after explo out
            prevIdx = findall(gd.stimulus[1:tin-1] .== s)
            postIdx = findall(gd.stimulus[tout+1:end] .== s) .+ tout
            for z = 1:before_after
                gd.beforeExplo[prevIdx[end-z+1]] = -z+1
                gd.afterExplo[postIdx[min(z, length(postIdx))]] = z
            end

            # Sort all trials in explo by presentation from explo out
            prevIdx = findall(gd.stimulus[1:tout-1] .== s) 
            for z = 1:maxInExplo
                zidx = prevIdx[max(length(prevIdx)-z+1, 1)]
                if zidx < tin # if this presentation is before entering explo
                    break 
                end
                gd.outExplo[zidx] = -z+1
            end
        end

        # Sort all trials in explo by presentation from explo in
        for z = 1:maxInExplo
            if tin + z - 1 > tout#gd.ExploOut[tin+z-1] == 1#
                break 
            end
            gd.inExplo[tin+z-1] = gd.presInBlock[tin+z-1]
        end

    end
end

exploWin = 7

# Before switching to exploration
gdf = groupby(df1[df1.condition .<= 4,:], [:subject, :sessNum, :beforeExplo])
cc = combine(gdf, :correct_or_persev => mean => :persev, :correct => mean => :cor, :zrt => mean => :rt)
before_df = combine(groupby(cc, :beforeExplo), :cor => mean, :cor => sem, :persev => mean, :persev => sem, :rt => mean, :rt => sem)

# During exploration
gdf = groupby(df1[df1.condition .<= 4,:], [:subject, :sessNum,:inExplo])
cc_in = combine(gdf, :persev => mean => :persev, :correct => mean => :cor, :zrt => mean => :rt)
in_df = combine(groupby(cc_in, :inExplo), :persev => mean, :persev => sem, :cor => mean, :cor => sem, :rt => mean, :rt => sem)

# After switching out of exploration
gdf = groupby(df1[df1.condition .<= 4,:], [:subject, :sessNum,:afterExplo])
cc = combine(gdf, :correct => mean => :cor, :persev => mean => :persev, :zrt => mean => :rt)
after_df = combine(groupby(cc, :afterExplo), :cor => mean, :cor => sem, :persev => mean, :persev => sem, :rt => mean, :rt => sem)

# During exploration (locked by the last exploration trial)
gdf = groupby(df1[df1.condition .<= 4,:], [:subject, :sessNum,:outExplo])
cc_out = combine(gdf, :correct => mean => :cor, :persev => mean => :persev, :zrt => mean => :rt)
out_df = combine(groupby(cc_out, :outExplo), :cor => mean, :cor => sem, :persev => mean, :persev => sem, :rt => mean, :rt => sem)

# All exploration trials
dur_df = combine(groupby(df1[df1.switchType .== -3,:], [:presInBlock, :condition]), :ExploOut => (x -> Float64(sum(x))) => :ExploOut)
for ci = 0:4
    dur_df.ExploOut[dur_df.condition .== ci] ./= sum(dur_df.ExploOut[dur_df.condition .== ci])
end


## RTs when switching in exploration
plot([-3, 8], repeat(summaryRT_group.zrt_mean', 2, 1), ribbon=repeat(summaryRT_group.zrt_sem', 2, 1), linewidth=3, linestyle=:dash, label="", color=:grey)

@df before_df[1 .> before_df.beforeExplo .> -4,:] plot!(:beforeExplo, :rt_mean, ribbon=:rt_sem, label="", linewidth = 3, color=:black, xticks=(-2:2:exploWin), size=(500, 500), background_color=:transparent, foreground_color=:black, labelfontsize=20, tickfontsize = 14, dpi=300, ylabel="RT (z-scored)", xlabel="Stimulus presentations")
@df in_df[0 .< in_df.inExplo .<= 8,:] plot!(:inExplo,  :rt_mean, ribbon=:rt_sem, label="", linewidth = 3, color=:black, xticks=(-2:2:exploWin))

## Statistical significance (using cluster-based permutation test)

X = fill(NaN, length(unique(cc_in.subject)), 8)
gp = groupby(cc_in, :subject)
for i in 1:length(gp)
    # find baseline for the subject 
    sidx = findfirst((summaryRT_sub.subject .== gp[i].subject[1]) .* (summaryRT_sub.correct .== false))
    baseline = summaryRT_sub.zrt[sidx]
    for j = 1:maximum(gp[i].inExplo)
        jidx = gp[i].inExplo[j]
        if 8 >= jidx > 0
            X[i,jidx] = gp[i].rt[j] - baseline
        end
    end
end
res = cluster_perm_test(X; niter=1e5)
# cluster : [1] p = 0.0078

plot!(res.clusters[1] .+ [-0.3, 0.3], 0.55 .* ones(length(res.clusters[1]) + 1), linewidth=5, label="", color =:black, alpha=0.5)


## RTs when switching out of exploration
plot([-4, 5], repeat(summaryRT_group.zrt_mean', 2, 1), ribbon=repeat(summaryRT_group.zrt_sem', 2, 1), linewidth=3, linestyle=:dash, label="", color=:grey)

tmp = out_df[out_df.outExplo .>= -4,:]
sort!(tmp, :outExplo)
@df tmp plot!(:outExplo, :rt_mean, ribbon=:rt_sem, label="", linewidth = 3, color=[:blue :black], xticks=(-6:2:6))
@df after_df[after_df.afterExplo .< 6,:] plot!(:afterExplo,:rt_mean, ribbon = :rt_sem, label="", linewidth = 3, color=[:blue :black], size=(500, 500), background_color=:transparent, foreground_color=:black, labelfontsize=20, tickfontsize = 14, dpi=300, ylabel="RT (z-scored)", xlabel="Stimulus presentations", legendfontsize=12)

X = fill(NaN, length(unique(cc_out.subject)), 5)
gp = groupby(cc_out, :subject)
for i in 1:length(gp)
    # find baseline for the subject 
    sidx = findfirst((summaryRT_sub.subject .== gp[i].subject[1]) .* (summaryRT_sub.correct .== false))
    baseline = summaryRT_sub.zrt[sidx]
    for j = 1:nrow(gp[i])
        jidx = gp[i].outExplo[j]
        if 0 >= jidx >= -4
            X[i,Int(jidx+5)] = gp[i].rt[j] - baseline
        end
    end
end
res = cluster_perm_test(X; niter=1e5)
# cluster : [2, 3, 4] (= (-3, -2, -1 from switch)), p < 0.001
plot!(res.clusters[1] .- 5, 0.55 .* ones(length(res.clusters[1])), linewidth=5, label="", color =:blue, alpha=0.5)


## Perseveration : global VS Overlapping switches  (complete)
compdf = df1[0 .< df1.condition .<= 3,:]
compdf.condition .+= 3 .* (abs.(compdf.switchType) .== 3) 
compdf.nextCondition[compdf.nextCondition .> 3] .= 0
grpstats = grp_stats_hmm(compdf)
grpstats = grpstats[0 .< grpstats.condition,: ] # Remove other condition for the pre-switch trials
grp_plot_hmm(grpstats, "persev", [1,2,3]; xlims=(-2, 9), xticks=-2:2:9,ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. perseverative choice", legend_position = :topright, linestyle=[:solid :solid :solid :dot :dot :dot ], background_color=:transparent, foreground_color=:black, labelfontsize=20, tickfontsize=14, dpi=300)
hline!([1/3], color=:black, linewidth=3, linestyle=:dash, label="", size=(500, 500), legendfontsize=14)

## Statistical Significance (cluster based permutation test)
tmp = df1[(0 .< df1.condition .<= 3) .* (df1.presInBlock .<= 10),:]
gp = groupby(tmp, [:presInBlock, :subject, :condition])
cc = combine(gp, [:persev, :switchType] => ((x, y) -> mean(x[abs.(y) .< 3]) - mean( x[abs.(y) .== 3])) => :diff1_2)

pl = Plots.current()
for ci = 1:3
    cci = cc[cc.condition .== ci,:]
    X = zeros(length(unique(cci.subject)), 10)
    gp = groupby(cci, :subject)
    for i in 1:length(gp)
        for j = 1:10
            jidx = gp[i].presInBlock[j]
            X[i,jidx] = gp[i].diff1_2[j]
        end
    end
    res = cluster_perm_test(X; niter=1e5)
    plot!(res.clusters[1], (0.95 - 0.05 * (ci-1)) .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[ci], alpha=0.5)
end
plot(pl)


## Global switches in condition 3 (partial) for stable stims
alldf = copy(df1)
grpstats = grp_stats_hmm(alldf)
grpstats = grpstats[0 .< grpstats.condition,: ] 
grp_plot_hmm(grpstats, "persev", [1,2,3,4]; xlims=(-2, 9), xticks=-2:2:9,ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. perseverative choice", legend_position = :topright, linestyle=[:solid :solid :solid :solid], background_color=:transparent, foreground_color=:black, labelfontsize=20, tickfontsize=14, dpi=300)
hline!([1/3], color=:black, linewidth=3, linestyle=:dash, label="", size=(500, 500), legendfontsize=14)

## Global switches in condition 3 (partial) for stable stims : split between switch types
partdf = df1[df1.condition .== 4,:]
partdf.persev .= partdf.correct#partdf.correct_or_persev
partdf.condition[(partdf.switchType .== 2) .| (partdf.switchType .== 1)] .= 1 # Recode conditions as switch types
partdf.condition[abs.(partdf.switchType) .== 3] .= 2
grpstats_part = grp_stats_hmm(partdf)
grpstats_part = grpstats_part[0 .< grpstats_part.condition .<= 2,: ]
grp_plot_hmm(grpstats_part, "persev", [4,4]; linestyle=[:solid :dot],xlims=(-2, 9), xticks=-2:2:9, ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. perseverative choice", size=(500, 500), legend_position = :topright, background_color=:transparent, foreground_color=:black, labelfontsize=20, tickfontsize = 14, dpi=300)

## Statistical Significance (cluster based permutation test)
tmp = partdf[partdf.presInBlock .<= 10,:]
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, [:correct, :condition] => ((x, y) -> mean(x[y .== 1]) - mean( x[y .== 2])) => :diff1_2)

X = zeros(length(unique(cc.subject)), 10)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:10
        jidx = gp[i].presInBlock[j]
        X[i,jidx] = gp[i].diff1_2[j]
    end
end
res = cluster_perm_test(X; niter=1e5)
plot!(res.clusters[1], 0.95 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[4], alpha=0.5)


hline!([1/3], color=:black, linewidth=3, linestyle=:dash, label="")

## Perseveration : switch in random exploration

@df before_df[1 .> before_df.beforeExplo .> -4,:] plot(:beforeExplo, :persev_mean, ribbon=:persev_sem, label="", linewidth = 3, color=:black, xticks=(-2:2:exploWin), size=(500, 500), background_color=:transparent, foreground_color=:black, labelfontsize=20, tickfontsize = 14, dpi=300, ylabel="Prop. perseverative choice", xlabel="Stimulus presentations")
@df in_df[0 .< in_df.inExplo .<= exploWin,:] plot!(:inExplo,  :persev_mean, ribbon=:persev_sem, label="", linewidth = 3, color=:black, xticks=(-2:2:exploWin), ylims=(0, 1.0))
hline!([1/3], color=:black, linewidth=3, linestyle=:dash, label="")

##
@df dur_df[(dur_df.presInBlock .<= exploWin) .* (0 .< dur_df.condition .<= 2),:] plot(:presInBlock, :ExploOut, group=:condition, linewidth = 5, palette=StatsPlots.palette(:Dark2)[1:2], alpha=1.0, label="", xticks=(2:2:exploWin), xlims=(0, 9), yticks=(0:0.1:0.3), ylims=(0, 0.45),labelfontsize=20, tickfontsize = 14, xlabel="Stimulus presentations", ylabel="Probability of switching out\nof the random strategy", background_color=:transparent, foreground_color=:black, size=(500, 500), dpi=300)


##
tmp = out_df[out_df.outExplo .>= -4,:]
sort!(tmp, :outExplo)
@df tmp plot(:outExplo, :cor_mean, ribbon=:cor_sem, label="", linewidth = 3, color=[:blue :black], ylims=(0, 1.0), xticks=(-6:2:6))
@df after_df[after_df.afterExplo .< 6,:] plot!(:afterExplo,:cor_mean, ribbon = :cor_sem, label="", linewidth = 3, color=[:blue :black], ylabel="Prop. correct choice", xlabel="Stimulus presentations", legendfontsize=12)
hline!([1/3], color=:black, linewidth=3, linestyle=:dash, label="", size=(500, 500), background_color=:transparent, foreground_color=:black, labelfontsize=20, tickfontsize = 14, dpi=300)

## Latency
gdf = groupby(df1, [:subject, :blockNum, :sessNum])
latdf = DataFrame(subject = zeros(Int, length(gdf)), blockNum = zeros(Int, length(gdf)), condition = zeros(Int, length(gdf)), switchType = zeros(Int, length(gdf)), latency_trials = zeros(Int, length(gdf)), latency_pres = zeros(length(gdf)))

for gi in 1:length(gdf)
    g = gdf[gi]
    latdf[gi, :subject] = g.subject[1]
    latdf[gi, :blockNum] = g.blockNum[1]
    latdf[gi, :condition] = g.condition[1] < 4 ? g.condition[1] : 3 # Reclassify stable associations in condition 3 g.condition[1]#
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


sdf = groupby(latdf, [:subject, :condition, :switchType])
latsum = combine(sdf, :latency_pres => mean => :latency_mean)
latdf1 = latsum[latsum.condition .== 3,:] # for later
latdf1.condition .= 1

latpop = combine(groupby(latsum, [:condition, :switchType]), :latency_mean => median => :latency_median, :latency_mean =>( x -> quantile(x, 0.25)) => :latency_q1, :latency_mean =>( x -> quantile(x, 0.75)) => :latency_q3)
latpop.err = (latpop.latency_q3 .- latpop.latency_q1) ./ 2
latpop.cerr =  (latpop.latency_q3 .+ latpop.latency_q1) ./ 2

lpglobal = latpop[latpop.switchType .== 3, :]
lsglobal = latsum[latsum.switchType .== 3,:]
@df lpglobal bar(:condition, :latency_median, group=:condition, label = "", ylabel = "Rule change → Strategy change\nlatency (# presentations)", palette=StatsPlots.palette(:Dark2)[1:3],alpha=0.6, linewidth=0, size=(500, 500),  dpi=300, background_color = :transparent, foreground_color=:black, labelfontsize=16, tickfontsize=12)

bar!(latpop[latpop.switchType .== 2, :condition] .- 4, latpop[latpop.switchType .== 2, :latency_median], group=latpop[latpop.switchType .== 2, :condition],label = "", alpha=0.6, linewidth=0)

bar!(latpop[latpop.switchType .== 1, :condition] .- 8, latpop[latpop.switchType .== 1, :latency_median], group=latpop[latpop.switchType .== 1, :condition],label = "", alpha=0.6, linewidth=0)

bar!(latpop[latpop.switchType .== -3, :condition] .+ 4, latpop[latpop.switchType .== -3, :latency_median], group=latpop[latpop.switchType .== -3, :condition],label = "", alpha=0.6, linewidth=0)


@df lsglobal dotplot!(:condition, :latency_mean, color=:grey, msw=0, markersize=3,label="", bar_width=0.3)
@df lpglobal scatter!(:condition, :cerr, label="", yerror=:err, color=:black, markersize=0, linewidth=5)

lpoverlap = latpop[latpop.switchType .== 2, :]
lsoverlap = latsum[latsum.switchType .== 2,:]
lpoverlap.condition .-= 4
lsoverlap.condition .-= 4
@df lsoverlap dotplot!(:condition, :latency_mean, color=:grey, msw=0, markersize=3,label="", bar_width=0.3)
@df lpoverlap scatter!(:condition, :cerr, label="", yerror=:err, color=:black, markersize=0, linewidth=5)

lpoverlap = latpop[latpop.switchType .== 1, :]
lsoverlap = latsum[latsum.switchType .== 1,:]
lpoverlap.condition .-= 8
lsoverlap.condition .-= 8
@df lsoverlap dotplot!(:condition, :latency_mean, color=:grey, msw=0, markersize=3,label="", bar_width=0.3)
@df lpoverlap scatter!(:condition, :cerr, label="", yerror=:err, color=:black, markersize=0, linewidth=5)

lprandom = latpop[latpop.switchType .== -3, :]
lsrandom = latsum[latsum.switchType .== -3,:]
lprandom.condition .+= 4
lsrandom.condition .+= 4
@df lsrandom dotplot!(:condition, :latency_mean, color=:grey, msw=0, markersize=3,label="", bar_width=0.3)
@df lprandom scatter!(:condition, :cerr, label="", yerror=:err, color=:black, markersize=0, linewidth=5)

xticks!([-6, -2, 2, 6], ["2 simil.", "1 simil.", "0 simil.", "Random"])
xlabel!("New behavioral strategy")


## Stats 
# Latency per condition 
latsub = combine(groupby(latsum, [:subject, :condition]), :latency_mean => mean => :latency_mean)
X1 = latsub.latency_mean[latsub.condition .== 1]
X2 = latsub.latency_mean[latsub.condition .== 2]
X3 = latsub.latency_mean[latsub.condition .== 3]
SignedRankTest(X1, X2)

# Per condition x switchType
gdf  = groupby(latsum, [:subject, :switchType])
subidx = unique(latsum.subject)
sw = [-3,1,2,3]
conds = [1,2,3]
L = zeros(length(subidx), length(sw), length(conds))
for i in eachindex(subidx)
    for j in eachindex(sw)
        for k in eachindex(conds)
            id = findfirst((latsum.subject .== subidx[i]) .* (latsum.condition .== conds[k]) .* (latsum.switchType .== sw[j]))
            if !isnothing(id)
                L[i,j,k] = latsum.latency_mean[id]
            else
                L[i,j,k] = NaN
            end
        end
    end
end

## Random 
X1 = L[:,1,1]
X2 = L[:,1,2]
X3 = L[:,1,3]
## Random : New vs recurrent
X = X1 .- X2 
SignedRankTest(filter(!isnan, X))
## Random : New vs partial
X = X1 .- X3
SignedRankTest(filter(!isnan, X))

## Global 
X1 = L[:,4,1]
X2 = L[:,4,2]
X3 = L[:,4,3]
## Global : New vs recurrent
X = X1 .- X2 
SignedRankTest(filter(!isnan, X))
## Global : New vs partial
X = X1 .- X3
SignedRankTest(filter(!isnan, X))

## 1 simil 
X1 = L[:,3,1]
X2 = L[:,3,2]
X3 = L[:,3,3]
## 1 simil : New vs recurrent
X = X1 .- X2 
SignedRankTest(filter(!isnan, X))
## 1 simil : New vs partial
X = X1 .- X3
SignedRankTest(filter(!isnan, X))

## 2 simil 
X1 = L[:,2,1]
X2 = L[:,2,2]
X3 = L[:,2,3]
## 2 simil : New vs recurrent
X = X1 .- X2 
SignedRankTest(filter(!isnan, X))
## 2 smil : New vs partial
X = X1 .- X3
SignedRankTest(filter(!isnan, X))


## Random exploration duration per condition
randdf = df1[(df1.switchType .== -3) .* (3 .> df1.condition .> 0),:]
gdf = groupby(randdf, [:subject, :switchNum, :sessNum])
latdf = DataFrame(subject = zeros(Int, length(gdf)), blockNum = zeros(Int, length(gdf)), condition = zeros(Int, length(gdf)), switchType = zeros(Int, length(gdf)), latency_trials = zeros(Int, length(gdf)), latency_pres = zeros(length(gdf)))

for gi in 1:length(gdf)
    g = gdf[gi]
    latdf[gi, :subject] = g.subject[1]
    latdf[gi, :blockNum] = g.blockNum[1]
    latdf[gi, :condition] = g.condition[1]#g.condition[1] < 4 ? g.condition[1] : 3 # Reclassify stable associations in condition 3
    latdf[gi, :switchType] = g.switchType[1]
    idx = findfirst(g.ExploOut .== 1)
    if !isnothing(idx)
        latdf[gi, :latency_trials] = idx
        mean([g.stimulus[1:idx] .== stim for stim = 1:3])
        latdf[gi, :latency_pres] = mean([count(g.stimulus[1:idx] .== stim) for stim = 1:3]) # Averagre presentations from the begining of the block to the switch
    else
        latdf[gi, :latency_trials] = -1
        latdf[gi, :latency_pres] = -1
    end
end
latdf = latdf[latdf.latency_pres .> 0,:]
sdf = groupby(latdf, [:subject, :condition])
latsum = combine(sdf, :latency_pres => mean => :latency_mean)
ss = combine(groupby(latsum, :condition), :latency_mean => mean => :latency_mean, :latency_mean => sem => :latency_sem)
bar(ss.condition, ss.latency_mean, color=StatsPlots.palette(:Dark2)[1:2], alpha=0.6,  linewidth=0, bar_width=0.8, label="")
dotplot!(latsum.condition, latsum.latency_mean, group=latsum.condition, label = "", color=:grey, msw=0, markersize=5/2, bar_width=0.3, xticks=[], yticks=(2:2:8, []), ylims=(0, 9), ymirror=true, yrotation=90, tickfontsize=14, thickness_scaling=2, size=(200, 700), background_color=:transparent, foreground_color=:black, dpi=300)
scatter!( ss.condition, ss.latency_mean,yerror=ss.latency_sem, color=:black, markerstrokewidth=10/2, markersize=0, label="")


## Switch type per condition
## Env 1
df1[df1.condition .== 4, :condition] .= 3 # Re-merge stable and changing stims in condition 3
stdf = df1[df1.newBlock .== 1,:]
stdf = stdf[stdf.condition .> 0,:]
gdf = groupby(stdf, [:subject, :condition])

stdf1 = combine(gdf, :switchType => (x -> mean(x .== -3)) => :random, :switchType => (x -> mean(x .== 3)) => :global, :switchType => (x -> mean(x .== 2)) => :paired, :switchType => (x -> mean(x .== 1)) => :local, :switchType => (x -> mean(0 .< x .<= 2)) => :overlapping, :switchType => (x -> mean(abs.(x) .== 3)) => :nonoverlapping)
stdf1 = stdf1[stdf1.condition .> 0,:]
##
# @df stdf1 violin(:condition, :overlapping, group=:condition, label="", ylims=(0, 1), ylabel="Proportion of strategic changes", palette=StatsPlots.palette(:Dark2)[1:3], alpha=0.6, linewidth=0, size=(500, 500), tickfontsize=14, labelfontsize=18, dpi=300)
bar_sum=combine(groupby(stdf1, :condition), :random => mean, :random => sem, :global => mean, :global => sem, :paired => mean, :paired => sem, :local => mean, :local => sem)

bar(bar_sum.condition, bar_sum.local_mean, yerror=bar_sum.local_sem, label="", ylims=(0, 1), ylabel="Proportion of strategic changes", color=StatsPlots.palette(:Dark2)[1:3], alpha=0.6, linewidth=0, msw=5, msc=:black, size=(500, 500), tickfontsize=12, labelfontsize=18, background_color=:transparent, foreground_color=:black, dpi=300)
dotplot!(stdf1.condition, stdf1.local, group=stdf1.condition, color=:grey, msw=0, markersize=3,label="", ylims=(0, 1), bar_width=0.3, ylabel="Proportion of strategic changes")
scatter!(bar_sum.condition, bar_sum.local_mean, yerror=bar_sum.local_sem, label="", markersize=0, msc=:black, msw=5)


bar!(bar_sum.condition .+ 4, bar_sum.paired_mean, yerror=bar_sum.paired_sem, label="", ylims=(0, 1), ylabel="Proportion of strategic changes", color=StatsPlots.palette(:Dark2)[1:3], alpha=0.6, linewidth=0, msw=5, msc=:black, size=(500, 500), tickfontsize=12, labelfontsize=18, background_color=:transparent, foreground_color=:black, dpi=300)
dotplot!(stdf1.condition.+4, stdf1[:,:paired], group=stdf1.condition, color=:grey, msw=0, markersize=3, label="", ylims=(0, 1), bar_width=0.3,ylabel="Proportion of strategic changes")
scatter!(bar_sum.condition.+4, bar_sum.paired_mean, yerror=bar_sum.paired_sem, label="", markersize=0, msc=:black, msw=5)


bar!(bar_sum.condition .+ 8, bar_sum.global_mean, yerror=bar_sum.global_sem, label="", ylims=(0, 1), ylabel="Proportion", color=StatsPlots.palette(:Dark2)[1:3], alpha=0.6, linewidth=0, msw=5, msc=:black, size=(500, 500), tickfontsize=12, labelfontsize=18, background_color=:transparent, foreground_color=:black, dpi=300)
dotplot!(stdf1.condition .+8, stdf1[:,:global], group=stdf1.condition, color=:grey, msw=0, markersize=3, label="", ylims=(0, 1), bar_width=0.3,ylabel="Proportion")
scatter!(bar_sum.condition .+8, bar_sum.global_mean, yerror=bar_sum.global_sem, label="", markersize=0, msc=:black, msw=5)


bar!(bar_sum.condition .+ 12, bar_sum.random_mean, yerror=bar_sum.random_sem, label="", ylims=(0, 1), ylabel="Proportion", color=StatsPlots.palette(:Dark2)[1:3], alpha=0.6, linewidth=0, msw=5, msc=:black, size=(500, 500), tickfontsize=12, labelfontsize=18, background_color=:transparent, foreground_color=:black, dpi=300)
dotplot!(stdf1.condition .+12, stdf1[:,:random], group=stdf1.condition, color=:grey, msw=0, markersize=3, label="", ylims=(0, 1), bar_width=0.3,ylabel="Proportion")
scatter!(bar_sum.condition .+12, bar_sum.random_mean, yerror=bar_sum.random_sem, label="", markersize=0, msc=:black, msw=5)

xticks!([2, 6, 10, 14], ["2 simil.", "1 simil.", "0 simil.", "Random"])
xlabel!("New behavioral strategy")

## Function adding the HMM profile of the data to plot over simulations 
summary_data_hmm = copy(bar_sum)
function add_data_plot(summary_data_hmm)
    for i = 1:3
        scatter!([i - 0.3], [summary_data_hmm.local_mean[i]], yerror =  [summary_data_hmm.local_sem[i]], label="",  color =  StatsPlots.palette(:Dark2)[i], msc = StatsPlots.palette(:Dark2)[i], msw=5)

        scatter!([4 + i - 0.3], [summary_data_hmm.paired_mean[i]], yerror =  [summary_data_hmm.paired_sem[i]], label="",  color =  StatsPlots.palette(:Dark2)[i], msc = StatsPlots.palette(:Dark2)[i], msw=5)

        scatter!([8 + i - 0.3], [summary_data_hmm.global_mean[i]], yerror =  [summary_data_hmm.global_sem[i]], label="",  color =  StatsPlots.palette(:Dark2)[i], msc = StatsPlots.palette(:Dark2)[i], msw=5)

        scatter!([12 + i - 0.3], [summary_data_hmm.random_mean[i]], yerror =  [summary_data_hmm.random_sem[i]], label="",  color =  StatsPlots.palette(:Dark2)[i], msc = StatsPlots.palette(:Dark2)[i], msw=5)
    end
    plot!()
end

## Only partial switches
std_p = stdf1[stdf1.condition .== 3, [:overlapping, :nonoverlapping]]


hline([0.5], color=:black, linewidth=3, linestyle= :dash, label="", size=(500, 500))
@df std_p violin!([1], :overlapping, label="", xticks=(), ylims=(0, 1), ylabel="", palette=StatsPlots.palette(:Dark2)[4:4], alpha=0.6, linewidth=0, size=(500, 500), tickfontsize=14, labelfontsize=20, dpi=300)
@df std_p dotplot!([1], :overlapping,label="", xticks=[], xlims = (0.2, 1.8), ylims=(0, 1), ylabel="", palette=StatsPlots.palette(:Dark2)[4:4], alpha=0.6, linewidth=0, size=(500, 500), dpi=300, tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)

