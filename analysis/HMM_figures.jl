using CSV, DataFrames, StatsPlots, Statistics

## Hyperparameters first paradigm
flist = readdir("WMM/HMM/Hyperparameters")
df = DataFrame()

for f in flist
    tmp = CSV.read(string("WMM/HMM/Hyperparameters/", f), DataFrame)
    tmp[!, :Subject] .= f
    append!(df, tmp)
end 

volsA = [[],[]]
stosA = [[],[]]
distA = [[],[]]
fulldistA = [[],[]]

volsO = [[],[]]
stosO = [[],[]]
distO = [[],[]]
fulldistO = [[],[]]

gdf = groupby(df, [:Subject, :Agent, :Task])
for g in gdf
    t = g[in.(g.parameters, [["p[1]", "p[2]", "p[3]"]]),:mean]
    v = sum(t)
    s = g[g.parameters .== "ρ",:mean]
    
    if g.Task[1] == "WMM1" && g.Agent[1] == "act"
        push!(volsA[1], v)
        push!(stosA[1], s)
        push!(distA[1], t ./ sum(t))
        push!(fulldistA[1], vcat(1 - v, t))
    elseif g.Task[1] == "WMM2" && g.Agent[1] == "act"
        push!(volsA[2], v)
        push!(stosA[2], s)
        push!(distA[2], t ./ sum(t))
        push!(fulldistA[2], vcat(1 - v, t))
    elseif g.Task[1] == "WMM1" && g.Agent[1] == "obs"
        push!(volsO[1], v)
        push!(stosO[1], s)
        push!(distO[1], t ./ sum(t))
        push!(fulldistO[1], vcat(1 - v, t))
    elseif g.Task[1] == "WMM2" && g.Agent[1] == "obs"
        push!(volsO[2], v)
        push!(stosO[2], s)
        push!(distO[2], t ./ sum(t))
        push!(fulldistO[2], vcat(1 - v, t))
    else
        println("that shouldn't happen")
    end
end



## Fig. 1 : Volatility
scatter(volsO[1], volsA[1], label="", color=:blue, alpha=0.4)
scatter!([mean(volsO[1])], [mean(volsA[1])], xerror = std(volsO[1])/sqrt(101), yerror = std(volsA[1])/sqrt(101),color=:black, markerstrokecolor=:blue, markerstrokewidth=2, label="Task 1")
scatter!(volsO[2], volsA[2], label="", color=:pink, alpha=0.5)
scatter!([mean(volsO[2])], [mean(volsA[2])], xerror = std(volsO[2])/sqrt(103), yerror = std(volsA[2])/sqrt(103),color=:red, markerstrokecolor=:pink, markerstrokewidth=2, label="Task 2")
plot!([0, 0.4], [0, 0.4], label="", color=:black, linestyle=:dash)
xlabel!("Ideal Observer")
ylabel!("Subjects")
title!("Volatility")
xlims!(0, 0.15)


## Fig. 2 : Stochasticity
scatter(stosO[1], stosA[1], label="", color=:blue, alpha=0.4)
scatter!([mean(stosO[1])], [mean(stosA[1])], xerror = std(stosO[1])/sqrt(101), yerror = std(stosA[1])/sqrt(101),color=:black, markerstrokecolor=:blue, markerstrokewidth=2, label="Task 1")
scatter!(stosO[2], stosA[2], label="", color=:pink, alpha=0.5)
scatter!([mean(stosO[2])], [mean(stosA[2])], xerror = std(stosO[2])/sqrt(103), yerror = std(stosA[2])/sqrt(103),color=:red, markerstrokecolor=:pink, markerstrokewidth=2, label="Task 2")
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

## Figure 4 : Divergence VS performance
## Data
flist = readdir("WMM/HMM2/Data")
datadf = DataFrame()

for f in flist
    tmp = CSV.read(string("WMM/HMM2/Data/", f), DataFrame)
    tmp[!, :Subject] .= f
    tmp.sessNum .= parse(Int, f[7])
    tmp[!,:zrt] = (tmp.rt .- mean(tmp.rt)) ./ std(tmp.rt)
    append!(datadf, tmp)
end 
##
gdfd = groupby(datadf, [:Subject, :task])
perfdf = combine(gdfd, :correct => mean)

perf1 = perfdf[perfdf.task .== "WMM1",:correct_mean]
perf2 = perfdf[perfdf.task .== "WMM2",:correct_mean]

kldiv1 = [sum(fulldistA[1][i] .* log.(fulldistA[1][i] ./ fulldistO[1][i])) for i in eachindex(fulldistA[1])] 
kldiv2 = [sum(fulldistA[2][i] .* log.(fulldistA[2][i] ./ fulldistO[2][i])) for i in eachindex(fulldistA[2])] 

plot([0, 0.5], 0.703 .- 0.952 .* [0, 0.5], color = :blue, linewidth = 2, linestyle = :dash, label="Task 1")
plot!([0, 0.5], 0.745 .- 1.864 .* [0, 0.5], color = :red, linewidth = 2, linestyle = :dash, label="Task 2")

scatter!(kldiv1, perf1, color=:blue, alpha=0.4, label="")
scatter!(kldiv2, perf2, color=:pink, alpha=0.6, label="")

scatter!([mean(kldiv1)], [mean(perf1)], xerror = std(kldiv1)/sqrt(101), yerror = std(perf1)/sqrt(101), color=:black, markerstrokecolor=:blue, markerstrokewidth=3, label="")
scatter!([mean(kldiv2)], [mean(perf2)], xerror = std(kldiv2)/sqrt(103), yerror = std(perf2)/sqrt(103), color=:pink, markerstrokecolor=:red, markerstrokewidth=3, label="")

xlabel!("DKL[P(subject)|| P(ideal observer)]")
ylabel!("Performance")
ylims!(0.31, 0.86)
# xlims!(-0.01, 0.3)

## Figure 5 : Performance locked on internal switch
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
    # eidx = g.ExploOut .== 1 # Find all exploOut 
    # for i in eachindex(eidx)
    #     if eidx[i] && g.switchType[i] ≠ -3 # Remove all non random strategies
    #         eidx[i] = false
    #     end
    # end
    # idx = diff(vcat(0, findall(eidx) .- 1, nrow(g)))

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
## For task 1
df1 = datadf[datadf.task .== "WMM1",:]
df1.condition[df1.condition .== 4] .= 1 # Merge conditions 1 and 4
#df1.condition[df1.stableAS] .= 4
##
for t = 1:nrow(df1)
    if df1[t, Symbol("isStable_$(df1.stimulus[t])")] && df1.condition[t] > 0 
        df1.condition[t] = 4 # Make stable stims a special condition
    end
end
#df1[(df1.condition .== 4) .* (df1.presInBlock .<= 10), :persev] .= df1[(df1.condition .== 4) .* (df1.presInBlock .<= 10), :correct] # persev = correct for the first presentations in condition 4 (stable associations). Otherwise perseveration would be = 0
## Show perseveration locked on strategic changes
alldf = copy(df1)
grpstats = grp_stats_hmm(alldf)
grpstats = grpstats[0 .< grpstats.condition,: ] 
grp_plot_hmm(grpstats, "persev", [1,2,3,4]; xlims=(-2, 9), xticks=-2:2:9,ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. perseverative choice", legend_position = :topright, linestyle=[:solid :solid :solid :solid], background_color=:transparent, foreground_color=:black, labelfontsize=20, tickfontsize=14, dpi=300)
hline!([1/3], color=:black, linewidth=3, linestyle=:dash, label="", size=(500, 500), legendfontsize=14)

## Statistical Significance (cluster based permutation test)
# Stable effect (compared to pre switch performance)
M = combine(groupby(alldf[(alldf.negpresInBlock .>= -2) .* (alldf.condition .> 0),:], :subject), :correct_or_persev => mean => :persev_mean)
tmp = alldf[(alldf.condition .== 4) .* (alldf.presInBlock .<= 10),:]
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, :persev => mean)
X = zeros(length(unique(cc.subject)), 10)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:10
        jidx = gp[i].presInBlock[j]
        X[i,jidx] = gp[i].persev_mean[j] - M.persev_mean[i]
    end
end
res = cluster_perm_test(X; niter=1e5)
plot!(0:9, fill(mean(M.persev_mean), 10), color = :grey, linewidth=3, linestyle=:dash, label="", alpha=0.5)
plot!(res.clusters[1], 0.95 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[4], alpha=0.5)


## Global VS Overlapping switches  (complete)
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

#savefig("WMM/Figures/HMM/expe1_correct_exploOut.png")
#savefig("WMM/Figures/HMM/expe1_paired_complete.pdf")
#savefig("WMM/Figures/HMM/expe1_global_complete.pdf")
#savefig("WMM/Figures/HMM/expe1_random_complete.png")
#savefig("WMM/Figures/HMM/expe1_globalOrRandom_complete.pdf")

## Partial paired switches in condition 3
compdf = df1[df1.switchType .== 2,:]
# compdf.condition[compdf.switchType .== 2] .= 1
# compdf.nextCondition[compdf.switchType .== 2] .= 1
# compdf.condition[compdf.switchType .== 3] .= 2
# compdf.nextCondition[compdf.switchType .== 3] .= 2
grpstats = grp_stats_hmm(compdf)
grpstats = grpstats[grpstats.condition .== 3,: ] 
# Remove other condition for the pre-switch trials

grp_plot_hmm(grpstats, "persev", [3, 3]; xlims=(-2.1, 8.1), xticks=([-2.1, 0.1, 2.1, 4.1, 6.1, 8.1], [-3, 1, 3, 5, 7, 9]),ylims=(0, 1), yticks=0:0.2:1, label="Partial rule change\n(changing associations)", xlabel="Stimulus presentations", ylabel="Perseveration")
hline!([1/3], color=:black, linewidth=3, linestyle=:dash, label="Chance level", size=(500, 500))
#savefig("WMM/Figures/HMM/expe1_paired_partial.pdf")

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
#savefig("WMM/Figures/HMM/expe1_globalVSnonglobal_partial.png")

## Random exploration patterns
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
            if gd.ExploOut[tin+z-1] == 1
                break 
            end
            gd.inExplo[tin+z-1] = gd.presInBlock[tin+z-1]
        end

    end
end

exploWin = 9

gdf = groupby(df1[df1.condition .<= 4,:], [:subject, :sessNum, :beforeExplo])
cc = combine(gdf, :correct_or_persev => mean => :persev, :correct => mean => :cor)
before_df = combine(groupby(cc, :beforeExplo), :cor => mean, :cor => sem, :persev => mean, :persev => sem)

gdf = groupby(df1[df1.condition .<= 4,:], [:subject, :sessNum,:inExplo])
cc = combine(gdf, :persev => mean => :persev, :correct => mean => :cor)
in_df = combine(groupby(cc, :inExplo), :persev => mean, :persev => sem, :cor => mean, :cor => sem)

gdf = groupby(df1[df1.condition .<= 4,:], [:subject, :sessNum,:afterExplo])
cc = combine(gdf, :correct => mean => :cor, :persev => mean => :persev)
after_df = combine(groupby(cc, :afterExplo), :cor => mean, :cor => sem, :persev => mean, :persev => sem)

gdf = groupby(df1[df1.condition .<= 4,:], [:subject, :sessNum,:outExplo])
cc = combine(gdf, :correct => mean => :cor, :persev => mean => :persev)
out_df = combine(groupby(cc, :outExplo), :cor => mean, :cor => sem, :persev => mean, :persev => sem)


dur_df = combine(groupby(df1[df1.switchType .== -3,:], [:presInBlock, :condition]), :ExploOut => (x -> Float64(sum(x))) => :ExploOut)
for ci = 0:4
    dur_df.ExploOut[dur_df.condition .== ci] ./= sum(dur_df.ExploOut[dur_df.condition .== ci])
end
##
#@df before_df[1 .> before_df.beforeExplo .> -4,:] plot(:beforeExplo, :persev_mean, ribbon=:persev_sem, label="", linewidth = 3, color=:black, xticks=(-2:2:exploWin), size=(500, 500), background_color=:transparent, foreground_color=:black, labelfontsize=20, tickfontsize = 14, dpi=300, ylabel="Prop. perseverative choice", xlabel="Stimulus presentations")
#@df in_df[0 .< in_df.inExplo .<= exploWin,:] plot!(:inExplo,  :persev_mean, ribbon=:persev_sem, label="", linewidth = 3, color=:black, xticks=(-2:2:exploWin), ylims=(0, 1.0))
#@df dur_df[(dur_df.presInBlock .<= exploWin) .* (0 .< dur_df.condition .<= 2),:] plot(twinx(), :presInBlock, :ExploOut, group=:condition, linewidth = 5, palette=StatsPlots.palette(:Dark2)[1:2], alpha=1.0, label="", xticks=(-2:2:exploWin), ylims=(0, 0.35),labelfontsize=20, tickfontsize = 14, y_foreground_color_text=StatsPlots.palette(:Dark2)[2], y_foreground_color_border=StatsPlots.palette(:Dark2)[2], y_guidefontcolor=StatsPlots.palette(:Dark2)[2], ylabel="Probability of switching\nout of the random strategy")

@df dur_df[(dur_df.presInBlock .<= exploWin) .* (0 .< dur_df.condition .<= 2),:] plot(:presInBlock, :ExploOut, group=:condition, linewidth = 5, palette=StatsPlots.palette(:Dark2)[1:2], alpha=1.0, label="", xticks=(2:2:exploWin), xlims=(0, 9), yticks=(0:0.1:0.3), ylims=(0, 0.45),labelfontsize=20, tickfontsize = 14, xlabel="Stimulus presentations", ylabel="Probability of switching out\nof the random strategy", background_color=:transparent, foreground_color=:black, size=(500, 500), dpi=300)
#hline!([1/3], color=:black, linewidth=3, linestyle=:dash, label="")

##
tmp = out_df[out_df.outExplo .>= -4,:]
sort!(tmp, :outExplo)
@df tmp plot(:outExplo, :cor_mean, ribbon=:cor_sem, label="", linewidth = 3, color=[:blue :black], ylims=(0, 1.0), xticks=(-6:2:6))
@df after_df[after_df.afterExplo .< 6,:] plot!(:afterExplo,:cor_mean, ribbon = :cor_sem, label="", linewidth = 3, color=[:blue :black], ylabel="Prop. correct choice", xlabel="Stimulus presentations", legendfontsize=12)
hline!([1/3], color=:black, linewidth=3, linestyle=:dash, label="", size=(500, 500), background_color=:transparent, foreground_color=:black, labelfontsize=20, tickfontsize = 14, dpi=300)


## Task 2
df2 = datadf[datadf.task .== "WMM2",:]
##
for t = 1:nrow(df1)
    if df2[t, Symbol("isStable_$(df2.stimulus[t])")] && df2.condition[t] > 0
        df2.condition[t] += 5 # Make stable stims a special condition
    end
end
## Global switches  (complete)
compdf = df2[(df2.switchType .== -3) .* (0 .< df2.condition .<= 3),:]
grpstats = grp_stats_hmm(compdf)
grpstats = grpstats[0 .< grpstats.condition .<= 3,: ] # Remove other condition for the pre-switch trials
grp_plot_hmm(grpstats, "persev", [1,2,3]; xlims=(-2.1, 8.1), xticks=([-2.1, 0.1, 2.1, 4.1, 6.1, 8.1], [-3, 1, 3, 5, 7, 9]),ylims=(0, 1), yticks=0:0.2:1, label=["2/3 rule change" "1/3 rule change" "1/3 rule change + noise"], xlabel="Stimulus presentations", ylabel="Perseveration")
hline!([1/3], color=:black, linewidth=3, linestyle=:dash, label="Chance level", size=(500, 500))

#savefig("WMM/Figures/HMM/expe1_env2_global_complete.pdf")

## Global switches in condition 3 (partial) for stable stims : split between switch types
partdf = df2[df2.condition .== 6,:]
partdf.persev .= partdf.correct_or_persev
partdf.condition[(partdf.switchType .== 2) .| (partdf.switchType .== 1)] .= 1 # Recode conditions as switch types
partdf.condition[abs.(partdf.switchType) .== 3] .= 2
grpstats_part2 = grp_stats_hmm(partdf)
grpstats_part2 = grpstats_part2[0 .< grpstats_part2.condition .<= 2,: ]
grp_plot_hmm(grpstats_part2, "persev", [4,4]; linestyle=[:solid :dot],xlims=(-2, 9), xticks=-2:2:9, ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), legend_position = :topright, background_color=:transparent, foreground_color=:black, labelfontsize=20, tickfontsize = 14, dpi=300)
#savefig("WMM/Figures/HMM/expe1_env2_globalVSnonglobal_partial.png")

## Interference env1 VS env2
tmp = copy(grpstats_part2)
tmp.condition .+= 2
gg = vcat(grpstats_part, tmp)

grp_plot_hmm(gg, "persev", [3,3,5,5]; linestyle=[:solid :dot :solid :dot],xlims=(-2, 9), xticks=-2:2:9, ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), legend_position = :topright, background_color=:transparent, foreground_color=:black, labelfontsize=20, tickfontsize = 14, dpi=300)

## Trap sensitivity
sdf = groupby(df1, :subject)
trapmat = zeros(length(sdf), 4, 2) # P(error) at -1, +1, +2, +3, before and after the switch
for si in 1:length(sdf)
    g = sdf[si]
    trapidx = findall(g.trap .== 1)
    perror = zeros(4, 2)
    count = zeros(Int, 4, 2)
    for t in trapidx
        if g.stableAS[t] # Stable association only
            pre_post_switch = (g.presInBlock[t] < abs(g.negpresInBlock[t])) + 1 # 1 : pre-switch, 2 : post switch
            if (pre_post_switch == 1) || ((pre_post_switch == 2) && (abs(g.switchType[t]) < 3))  # if post switch : non global switch
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
cpal = StatsPlots.palette(:Dark2)
boxplot([0 2 4 6], trapmat[:,:,1], color=cpal[1], alpha=0.5, label=["Before switch" "" "" ""])
dotplot!([0 2 4 6], trapmat[:,:,1], color=cpal[1], alpha=0.5, label="")
boxplot!([1 3 5 7], trapmat[:,:,2], color=cpal[2], alpha=0.5, label=["After switch" "" "" ""])
dotplot!([1 3 5 7], trapmat[:,:,2], color=cpal[2], alpha=0.5, label="")
xticks!(0.5:2:6.5, ["-1", "+1", "+2", "+3"])
xlabel!("Presentations from trap")
ylabel!("P(error)")
title!("Trap sensitivity before and after a non global switch for stable associations", titlefontsize=10)

## Latency
## Env 1 
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
#latdf = latdf[latdf.switchType .== -3,:] # keep only global switches

sdf = groupby(latdf, [:subject, :condition, :switchType])
latsum = combine(sdf, :latency_pres => mean => :latency_mean)
latdf1 = latsum[latsum.condition .== 3,:] # for later
latdf1.condition .= 1

#latsum[0 .<= latsum.switchType .<= 2,:switchType] .= 2

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
savefig("WMM/Figures/HMM/expe1_latency2.png")

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


## Env 2

gdf = groupby(df2, [:subject, :blockNum, :sessNum])
latdf = DataFrame(subject = zeros(Int, length(gdf)), blockNum = zeros(Int, length(gdf)), condition = zeros(Int, length(gdf)), switchType = zeros(Int, length(gdf)), latency_trials = zeros(Int, length(gdf)), latency_pres = zeros(length(gdf)))

for gi in 1:length(gdf)
    g = gdf[gi]
    latdf[gi, :subject] = g.subject[1]
    latdf[gi, :blockNum] = g.blockNum[1]
    latdf[gi, :condition] = g.condition[1]
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
latdf2 = latsum[latsum.condition .== 3,:]
latdf2.condition .= 2
@df latsum[latsum.switchType .== -3,:] boxplot(:condition, :latency_mean, group=:condition, label = ["2/3 rule change" "1/3 rule change" "1/3 rule change + noise in stable" "No rule change + noise"], ylabel = "Rule change - Strategic switch latency\n(# presentations)", palette=StatsPlots.palette(:Dark2)[1:4], size=(500, 500))
boxplot!(latsum[latsum.switchType .== 3,:condition] .+ 5, latsum[latsum.switchType .== 3,:latency_mean], group=latsum[latsum.switchType .== 3,:condition], label = "")
boxplot!(latsum[latsum.switchType .== 2,:condition] .+ 10, latsum[latsum.switchType .== 2,:latency_mean], group=latsum[latsum.switchType .== 2,:condition], label = "")
xticks!([2.5, 7.5, 12.5], ["Random", "Global", "Overlapping"])
xlabel!("Switch type")
#savefig("WMM/Figures/HMM/expe1_env2_latency.pdf")

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


#savefig("WMM/Figures/HMM/expe1_env1_random_duration.pdf")
#savefig("WMM/Figures/HMM/expe1_env1_random_duration.png")

## Preswitch performance
# Task 1
idx = findall((datadf.blockNum .> 1) .* (datadf.newBlock .== 1) .* (datadf.task .== "WMM1"))
predf1 = DataFrame(subject = zeros(Int, length(idx)), correct = zeros(length(idx)), persev = zeros(length(idx)), explo = zeros(length(idx)), condition = zeros(Int, length(idx)), switchType = zeros(Int, length(idx)))

for i in eachindex(idx)
    predf1[i,:subject] = datadf.subject[idx[i]]
    predf1[i,:condition] = datadf.condition[idx[i]]
    predf1[i,:switchType] = datadf.switchType[idx[i]]
    predf1[i,:correct] = mean(datadf.correct[idx[i]-10:idx[i]-1])
    predf1[i,:persev] = mean(datadf.persev[idx[i]-10:idx[i]-1])
    predf1[i,:explo] = mean(datadf.explo[idx[i]-10:idx[i]-1])
end

predf1.condition[predf1.condition .== 4] .= 1
predf1 = predf1[predf1.condition .> 0,:]

gdf = groupby(predf1, [:subject, :condition, :switchType])
predf1 = combine(gdf, :correct => mean => :correct, :persev => mean => :persev, :explo => mean => :explo)

predf1[!,:allcondition] = ["C$(predf1.condition[i])S$(predf1.switchType[i])" for i = 1:nrow(predf1)] 

## Reaction times
# Z-scoring per session
df1[!,:zrt] = zeros(nrow(df1))
gdf = groupby(df1, [:subject, :sessNum])
for g in gdf
    mrt = movstat(mean, g.rt, 20)
    srt = movstat(std, g.rt, 20)
    g.zrt .= (g.rt .- mrt) ./ srt
end

# Summary per subject
#df1[df1.condition .== 4, :condition] .= 3 # Re-merge stable and changing stims in condition 3

## Local switches, only switching stims
swdf = copy(df1) # Only the switching stims
for t in 1:nrow(swdf)
    swdf[t,:switchType] = swdf[t,Symbol("switch_$(swdf.stimulus[t])")] == 1 ? swdf[t,:switchType] : 3 
end 
grpstats = grp_stats_hmm(swdf[1 .< swdf.switchType .< 3,:])

grpstats = grpstats[0 .< grpstats.condition .< 4,:]
recidx = findall(grpstats.condition .== 2)
grpstats.condition[grpstats.condition .== 3] .= 2
grpstats.condition[recidx] .= 3
grp_plot_hmm(grpstats, "rt", [1,3,2]; xlims=(-2.1, 8.1), xticks=([-2.1, 0.1, 2.1, 4.1, 6.1, 8.1], [-3, 1, 3, 5, 7, 9]),  ylims=(-0.25, 0.5), label=["Complete rule change" "Partial rule change" "Recurent rule"], xlabel="Stimulus presentations", ylabel="RT (z-score)", size=(500, 500))
#savefig("WMM/Figures/HMM/expe1_localswitching_RT.pdf")
## Local switches, only nonswitching stims
swdf = copy(df1) # Only the nonswitching stims
for t in 1:nrow(swdf)
    swdf[t,:switchType] = swdf[t,Symbol("switch_$(swdf.stimulus[t])")] == 1 ? 3 : swdf[t,:switchType]
end 
grpstats = grp_stats_hmm(swdf[1 .< swdf.switchType .< 3,:])

grpstats = grpstats[0 .< grpstats.condition .< 4,:]
recidx = findall(grpstats.condition .== 2)
grpstats.condition[grpstats.condition .== 3] .= 2
grpstats.condition[recidx] .= 3
grp_plot_hmm(grpstats, "rt", [1,3,2]; xlims=(-2.1, 8.1), xticks=([-2.1, 0.1, 2.1, 4.1, 6.1, 8.1], [-3, 1, 3, 5, 7, 9]),  ylims=(-0.25, 0.5), label=["Complete rule change" "Partial rule change" "Recurent rule"], xlabel="Stimulus presentations", ylabel="RT (z-score)", size=(500, 500))
#avefig("WMM/Figures/HMM/expe1_localnonswitching_RT.pdf")
## Global switches
grpstats = grp_stats_hmm(df1[df1.switchType .== 3,:])
recidx = findall(grpstats.condition .== 2)
grpstats.condition[grpstats.condition .== 3] .= 2
grpstats.condition[recidx] .= 3
grpstats = grpstats[0 .< grpstats.condition .< 4,:]
grp_plot_hmm(grpstats, "rt", [1,3,2]; xlims=(-2.1, 8.1), xticks=([-2.1, 0.1, 2.1, 4.1, 6.1, 8.1], [-3, 1, 3, 5, 7, 9]),  ylims=(-0.25, 0.5), label=["Complete rule change" "Partial rule change" "Recurent rule"], xlabel="Stimulus presentations", ylabel="RT (z-score)", size=(500, 500))
#savefig("WMM/Figures/HMM/expe1_global_RT.pdf")

## Transition ambiguity
## Env 1
idx = findall(df1.HMMSwitch .== 1)
ambdf = DataFrame(fb = zeros(length(idx)), st = zeros(Int, length(idx)), cond=zeros(Int, length(idx)), subject = df1.subject[idx])
for i in eachindex(idx)
    cp = findlast(df1.newBlock[1:idx[i]])
    ambdf.fb[i] = mean(df1.fb[cp:idx[i]])
    ambdf.cond[i] = df1.condition[idx[i]]
    if df1.switchType[idx[i]] == -3
        ambdf.st[i] = 0
    elseif df1.switchType[idx[i]] == 3
        ambdf.st[i] = 1
    elseif   df1.switchType[idx[i]] == 2
        ambdf.st[i] = 2
    else
        ambdf.st[i] = 3
    end
end
ambdf = ambdf[ambdf.cond .> 0,:]
gdf = groupby(ambdf, [:subject, :st, :cond])
cc = combine(gdf, :fb => mean)
cc1 = cc[cc.cond .== 3,:]
@df cc[cc.st .== 1,:] violin(:cond, :fb_mean, group=:cond, label= "",ylabel="Prop. of positive feedbacks\nsince rule change", palette=StatsPlots.palette(:Dark2)[1:3], size=(500, 600), alpha=0.6, linewidth=0, tickfontsize=14, labelfontsize=20)
@df cc[cc.st .== 1,:] dotplot!(:cond, :fb_mean, group=:cond, label= "", palette=StatsPlots.palette(:Dark2)[1:3], background_color=:transparent, foreground_color=:black)

violin!(cc[cc.st .== 2,:cond] .+ 4, cc[cc.st .== 2,:fb_mean], group = cc[cc.st .== 2,:cond],label="", alpha=0.6, linewidth=0)
dotplot!(cc[cc.st .== 2,:cond] .+ 4, cc[cc.st .== 2,:fb_mean], group = cc[cc.st .== 2,:cond],label="")
violin!(cc[cc.st .== 0,:cond] .+ 8, cc[cc.st .== 0,:fb_mean], group = cc[cc.st .== 0,:cond],label="", alpha=0.6, linewidth=0)
dotplot!(cc[cc.st .== 0,:cond] .+ 8, cc[cc.st .== 0,:fb_mean], group = cc[cc.st .== 0,:cond],label="")
xticks!([2, 6, 10], ["Global", "Overlapping", "Random"])
xlabel!("Switch type")
# 
hline!([1/3], color=:black, linestyle=:dash, linewidth=3, label="")
#savefig("WMM/Figures/HMM/expe1_fb_preswitch.pdf")

## Env 2
idx = findall(df2.HMMSwitch .== 1)
ambdf = DataFrame(fb = zeros(length(idx)), st = zeros(Int, length(idx)), cond=zeros(Int, length(idx)), subject = df2.subject[idx])
for i in eachindex(idx)
    cp = findlast(df2.newBlock[1:idx[i]])
    ambdf.fb[i] = mean(df2.fb[cp:idx[i]])
    ambdf.cond[i] = df2.condition[idx[i]]
    if df2.switchType[idx[i]] == -3
        ambdf.st[i] = 0
    elseif df2.switchType[idx[i]] == 3
        ambdf.st[i] = 1
    elseif   df2.switchType[idx[i]] == 2
        ambdf.st[i] = 2
    else
        ambdf.st[i] = 3
    end
end
ambdf = ambdf[ambdf.cond .> 0,:]
gdf = groupby(ambdf, [:subject, :st, :cond])
cc = combine(gdf, :fb => mean)
cc2 = cc[cc.cond .== 1,:]
@df cc[cc.st .== 0,:] boxplot(:cond, :fb_mean, group=:cond, label= ["2/3 rule change" "1/3 rule change" "1/3 rule change + noise in stable" "No rule change + noise"],ylabel="Proportion of positive feedbacks\nfrom rule change to strategic switch", palette=StatsPlots.palette(:Dark2)[1:4], size=(500, 500))

boxplot!(cc[cc.st .== 1,:cond] .+ 5, cc[cc.st .== 1,:fb_mean], group = cc[cc.st .== 1,:cond],label="")

boxplot!(cc[cc.st .== 2,:cond] .+ 10, cc[cc.st .== 2,:fb_mean], group = cc[cc.st .== 2,:cond],label="")
xticks!([2, 6, 10], ["Random", "Global", "Overlapping"])
xlabel!("Switch type")
# 
hline!([1/3], color=:black, linestyle=:dash, linewidth=3, label="Chance level", ylims=(0, 1.2), legend=:topleft)
#savefig("WMM/Figures/HMM/expe1_env2_fb_preswitch.pdf")

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

#violin!(stdf1.condition .+ 4, stdf1[:,:global], label="", group=stdf1.condition, alpha=0.6, linewidth=0)
#dotplot!(stdf1.condition .+ 4, stdf1[:,:global], label="", group=stdf1.condition)
#violin!(stdf1.condition .+ 8, stdf1[:,:random], label="", group=stdf1.condition, alpha=0.6, linewidth=0)
#dotplot!(stdf1.condition .+ 8, stdf1[:,:random], label="", group=stdf1.condition)

## add summary stats 
# ss = combine(groupby(stdf1, :condition), :random => mean, :random => sem, :global => mean, :global => sem, :overlapping => mean, :overlapping => sem)
# scatter!(ss.condition, ss.overlapping_mean, yerror=ss.overlapping_sem, color=:black, markershape=:circle, markerstrokewidth=3, markersize=6, label="")
# scatter!(ss.condition .+ 4, ss.global_mean, yerror=ss.global_sem, color=:black, markershape=:circle, markerstrokewidth=3, markersize=6, label="")
# scatter!(ss.condition .+ 8, ss.random_mean, yerror=ss.random_sem, color=:black, markershape=:circle, markerstrokewidth=3, markersize=6, label="")

xticks!([2, 6, 10, 14], ["2 simil.", "1 simil.", "0 simil.", "Random"])
#xticks!([2, 6, 10], ["Global", "Overlapping"])
xlabel!("New behavioral strategy")
#savefig("WMM/Figures/HMM/expe1_switchtype_condition.png")

## Only partial switches
std_p = stdf1[stdf1.condition .== 3, [:overlapping, :nonoverlapping]]
# Z = vec(sum(Matrix(std_p), dims=2))
# std_p.overlapping ./= Z
# std_p.nonoverlapping ./= Z

hline([0.5], color=:black, linewidth=3, linestyle= :dash, label="", size=(500, 500))
@df std_p violin!([1], :overlapping, label="", xticks=(), ylims=(0, 1), ylabel="", palette=StatsPlots.palette(:Dark2)[4:4], alpha=0.6, linewidth=0, size=(500, 500), tickfontsize=14, labelfontsize=20, dpi=300)
@df std_p dotplot!([1], :overlapping,label="", xticks=[], xlims = (0.2, 1.8), ylims=(0, 1), ylabel="", palette=StatsPlots.palette(:Dark2)[4:4], alpha=0.6, linewidth=0, size=(500, 500), dpi=300, tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)


## Env 2
df2[df2.condition .> 5, :condition] .-= 5  # Re-merge stable and changing stims in condition 3
stdf = df2[df2.newBlock .== 1,:]
gdf = groupby(stdf, [:subject, :condition])

stdf2 = combine(gdf, :switchType => (x -> mean(x .== -3)) => :random, :switchType => (x -> mean(x .== 3)) => :global, :switchType => (x -> mean(x .== 2)) => :paired, :switchType => (x -> mean(0 .< x .<= 2)) => :overlapping, :switchType => (x -> mean(abs.(x) .== 3)) => :nonoverlapping)
stdf2 = stdf2[stdf2.condition .> 0,:]
##
@df stdf2 violin(:condition, :overlapping, group=:condition, label="", ylims=(0, 1), ylabel="Proportion of strategic changes", palette=StatsPlots.palette(:Dark2)[1:4], alpha=0.6, linewidth=0, size=(500, 500), tickfontsize=14, labelfontsize=18, dpi=300)
@df stdf2 dotplot!(:condition, :overlapping, group=:condition, label="", ylims=(0, 1), ylabel="Proportion of strategic changes", palette=StatsPlots.palette(:Dark2)[1:4],background_color=:transparent, :foreground_color=:black)
violin!(stdf2.condition .+ 4, stdf2[:,:global], label="", group=stdf2.condition, alpha=0.6, linewidth=0)
dotplot!(stdf2.condition .+ 4, stdf2[:,:global], label="", group=stdf2.condition)
violin!(stdf2.condition .+ 8, stdf2[:,:random], label="", group=stdf2.condition, alpha=0.6, linewidth=0)
dotplot!(stdf2.condition .+ 8, stdf2[:,:random], label="", group=stdf2.condition)
xticks!([2, 6, 10], ["Overlapping", "Global", "Random"])
#xticks!([2, 6, 10], ["Global", "Overlapping"])
xlabel!("New strategy")
#savefig("WMM/Figures/HMM/expe1_env2_switchtype_condition.pdf")

## Env 1 VS Env 2 prop global

tmp = stdf1[stdf1.condition .== 3,:]
tmp2 = stdf2[stdf2.condition .== 1,:]
tmp.condition .= 1
tmp2.condition .= 2
stdf = vcat(tmp, tmp2)

@df stdf violin(:condition, :overlapping, group=:condition, label= "", ylims=(0, 1), yticks=0:0.25:1,ylabel="Proportion of switches", palette=StatsPlots.palette(:Dark2)[[3,1]], size=(500, 500), legend_position=:topleft, alpha=0.6, linewidth=0, tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)
@df stdf dotplot!(:condition, :overlapping, group=:condition, label= "", palette=StatsPlots.palette(:Dark2)[[3,1]])

violin!(stdf.condition .+ 3, stdf[:,:global], group=stdf.condition, label="", alpha=0.6, linewidth=0)
dotplot!(stdf.condition .+ 3, stdf[:,:global], group=stdf.condition, label="")

violin!(stdf.condition .+ 6, stdf[:,:random], group=stdf.condition, label="", alpha=0.6, linewidth=0)
dotplot!(stdf.condition .+ 6, stdf[:,:random], group=stdf.condition, label="")

xticks!([1.5, 4.5, 7.5], ["Overlapping", "Global", "Random"])
xlabel!("Switch type")
#savefig("WMM/Figures/HMM/expe1_env1V2_prop_switches.pdf")

## Env1 VS Env2 latencies
latdf = vcat(latdf1, latdf2)

@df latdf[0 .< latdf.switchType .<= 2,:] violin(:condition, :latency_mean, group=:condition,label="", ylabel="Rule change to strategic change\n latency (# presentations)", palette=StatsPlots.palette(:Dark2)[[3,5]], size=(500, 500), alpha=0.6, linewidth=0, tickfontsize=14, labelfontsize=18, background_color=:transparent, foreground_color=:black)
@df latdf[0 .< latdf.switchType .<= 2,:] dotplot!(:condition, :latency_mean, group=:condition,label="", palette=StatsPlots.palette(:Dark2)[[3,5]], size=(500, 500))

violin!(latdf[latdf.switchType .== 3,:condition] .+ 3, latdf[latdf.switchType .== 3,:latency_mean], group = latdf[latdf.switchType .== 3,:condition], label="", alpha=0.6, linewidth=0)
dotplot!(latdf[latdf.switchType .== 3,:condition] .+ 3, latdf[latdf.switchType .== 3,:latency_mean], group = latdf[latdf.switchType .== 3,:condition], label="")

violin!(latdf[latdf.switchType .== -3,:condition] .+ 6, latdf[latdf.switchType .== -3,:latency_mean], group = latdf[latdf.switchType .== -3,:condition], label="", alpha=0.6, linewidth=0)
dotplot!(latdf[latdf.switchType .== -3,:condition] .+ 6, latdf[latdf.switchType .== -3,:latency_mean], group = latdf[latdf.switchType .== -3,:condition], label="", dpi=300)
xticks!([1.5, 4.5, 7.5], ["Overlapping", "Global", "Random"])
xlabel!("New strategy")
#savefig("WMM/Figures/HMM/expe1_env1V2_latency.png")

## Env1 VS Env2 fb pre switch
cc2.cond .= 2
cc1.cond .= 1

cc = vcat(cc1, cc2)
@df cc[cc.st .== 0,:] violin(:cond, :fb_mean, group=:cond,label="", ylabel="Prop. of positive feedbacks\nsince rule change", palette=StatsPlots.palette(:Dark2)[[3,1]], size=(500, 500), alpha=0.6, linewidth=0,tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black) 
@df cc[cc.st .== 0,:] dotplot!(:cond, :fb_mean, group=:cond,label="", palette=StatsPlots.palette(:Dark2)[[3,1]]) 

violin!(cc[cc.st .== 1,:cond] .+ 3, cc[cc.st .== 1,:fb_mean], group = cc[cc.st .== 1,:cond],label="", alpha=0.6, linewidth=0)
dotplot!(cc[cc.st .== 1,:cond] .+ 3, cc[cc.st .== 1,:fb_mean], group = cc[cc.st .== 1,:cond],label="")

violin!(cc[cc.st .== 2,:cond] .+ 6, cc[cc.st .== 2,:fb_mean], group = cc[cc.st .== 2,:cond],label="", alpha=0.6, linewidth=0)
dotplot!(cc[cc.st .== 2,:cond] .+ 6, cc[cc.st .== 2,:fb_mean], group = cc[cc.st .== 2,:cond],label="")
xticks!([1.5, 4.5, 7.5], ["Random", "Global", "Overlapping"])
xlabel!("Switch type")
#savefig("WMM/Figures/HMM/expe1_env1V2_fb_preswitch.pdf")