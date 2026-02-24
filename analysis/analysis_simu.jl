using CSV, DataFrames, ProgressMeter, StatsPlots

include("preprocessing.jl")
include("stats_plots_funs.jl")

## Load simus
model = "SI_MultVol_SampleAction"
task = 1
tasknames = ["WMM1", "WMM2", "ALL"]

pathname = string("../../Model_Tasks_Data/Data/WMM/Simus/all_simus/Simus/", model, "_", tasknames[task])
fileslist = filter(x -> occursin(".csv", x), readdir(pathname))
big_df = DataFrame()

wb = Progress(length(fileslist), 1, "loading...")
for fi in eachindex(fileslist)
    fl = fileslist[fi]
    df = CSV.read(string(pathname, "/", fl), DataFrame)
    df.subject .+= fi * 1000
    df[!,:zrt] = (df.rt .- mean(df.rt)) ./ std(df.rt)
    pres_in_block!(df)
    trials_in_block!(df)
    df[!,:nextCondition] .= 0
    gdf = groupby(df, [:subject, :sessNum, :blockNum])
    for gi = 1:length(gdf)
        if gdf[gi].blockNum[1] < 39 
            gdf[gi].nextCondition .= gdf[gi+1].condition[1]
        end
    end
    if in("ExploOut", names(df))
         append!(big_df, df)
    end
    next!(wb)
end
# Split by Environment
envdf = groupby(big_df, :task);

## Environment 1
df_env1 = DataFrame(envdf[[en.task[1] == "WMM1" for en in envdf]])

df_env1[!, :isStable] = falses(nrow(df_env1)) # Re-code stable associations from rule change
for t = 1:nrow(df_env1)
    if df_env1[t, Symbol("isStable_$(df_env1.stimulus[t])")] && df_env1.condition[t] > 0
        df_env1.isStable[t] = true
    end
end


# Check no effect of rare vs frequent associations 
grpstats0 = grp_stats(df_env1)

## Learning (correct)
grp_plot!(grpstats0[(grpstats0.condition .== 1) .| (grpstats0.condition .== 4),:], "correct", "false", [1,6], "bar"; xlims=(-3, 9), ylims=(0, 1), xticks=-2:2:9, yticks=0:0.2:1,  xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), label="", labelfontsize=20, tickfontsize=14,  background_color = :transparent, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black, dpi=300)

##
#df_env1.condition[df_env1.condition .== 4] .= 1 # Merge conditions 1 and 4
grpstats1 = grp_stats(df_env1);


## Recurrence effect (correct)
grp_plot!(grpstats1[0 .< grpstats1.condition .<= 2,:], "correct", "false", [1,2], "bar"; xlims=(-3, 9), xticks=(-2:2:9), ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300)

## Recurrence effect (explo)
grp_plot!(grpstats1[0 .< grpstats1.condition .<= 2,:], "explo", "false", [1, 2], "bar"; xlims=(-3, 9), xticks=(-2:2:9), ylims=(0,0.3), yticks=0:0.1:1, label="", xlabel="Stimulus presentations", ylabel="Non perseverative incorrect choice", size=(500, 500))

##
grp_plot(grpstats1[grpstats1.condition .< 3,:], "persev", [1,2]; xlims=(-3, 12), ylims=(0, 1), yticks=0:0.2:1, label=["New rule" "Recurent rule"], xlabel="Stimulus presentations", ylabel="Perseveration", size=(500, 500))

## Partial effect (correct)
grp_plot!(grpstats1[(grpstats1.condition .== 1) .|(grpstats1.condition .== 3),:], "correct", "false", [1, 3], "bar"; linestyle=[:solid :solid], xlims=(-3, 9), ylims=(0, 1),xticks=-2:2:9, yticks=0:0.2:1, xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), label="", labelfontsize=20, tickfontsize=14,  background_color = :transparent, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black, dpi=300)
grp_plot!(grpstats1[grpstats1.condition .== 3,:], "correct", "true", [4], "bar"; linestyle=:dot, label="")

## Partial vs complete (correct)
grp_plot(grpstats1[(grpstats1.condition .== 1) .| (grpstats1.condition .== 3),:], "correct", [1,3]; xlims=(-3, 9), ylims=(0, 1), yticks=0:0.2:1, label=["Complete rule change" "Partial rule change"], xlabel="Stimulus presentations", ylabel="Correct choice", size=(500, 500))
savefig("WMM/Figures/model_based/expe1_$(model)_partial_correct.pdf")
## Partial effect (explo)
grp_plot(grpstats1[(grpstats1.condition .== 1) .| (grpstats1.condition .== 3),:], "explo", [1,3]; xlims=(-3, 9), ylims=(0, 0.3), yticks=0:0.2:1, label=["New rule" "Partial changing"], xlabel="Stimulus presentations", ylabel="Perseveration", size=(500, 500))
savefig("WMM/Figures/model_based/expe1_$(model)_partial_explo.pdf")

## Interference effect
df_inter = copy(df_env1)
group_level_interference!(df_inter)
df_inter = df_inter[(df_inter.condition .== 3) .* (df_inter.stableAS) .* (df_inter.presInBlock .== 1),:]
gdf = groupby(df_inter, [:subject, :cumNegFB])
summary_inter = combine(gdf, :diffFromPlateau => mean => :correct)

gdf = groupby(summary_inter, [:cumNegFB])
summary_inter = combine(gdf, :correct => mean, :correct => sem)

summary_inter = summary_inter[.!isnan.(summary_inter.correct_sem),:]

@df summary_inter bar(:cumNegFB, :correct_mean, yerror = :correct_sem, label="",  ylims=(-0.6, 0.2), xlabel="# negative feedback since rule change", ylabel="Prop. correct choice", dpi=300, size=(500,500), background_color=:white)


## Environment 2
df_env2 = DataFrame(envdf[[en.task[1] == "WMM2" for en in envdf]])
df_env2[!, :isStable] = falses(nrow(df_env2)) # Re-code stable associations from rule change
for t = 1:nrow(df_env2)
    if df_env2[t, Symbol("isStable_$(df_env2.stimulus[t])")] && df_env2.condition[t] > 0
        df_env2.isStable[t] = true
    end
end
grpstats2 = grp_stats(df_env2)


## Relearning (correct)
grp_plot(grpstats2[0 .< grpstats2.condition .<= 2,:], "correct","false", [5,6]; xlims=(-3, 9),xticks=([-2, 3, 5,7, 9], [-3, 3, 5,7, 9]), ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)
#savefig("WMM/Figures/model_based/expe1_$(model)_relearn_correct.pdf")

## Relearning (explo)
grp_plot(grpstats2[0 .< grpstats2.condition .<= 2,:], "explo", [1,2]; xlims=(-3, 9), xticks=([-2, 3, 5,7, 9], [-3, 3, 5,7, 9]), ylims=(0, 0.3), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. non perseverative\nincorrect choice", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)
#savefig("WMM/Figures/model_based/expe1_$(model)_relearn_explo.pdf")

## Partial effect (stable)
grp_plot(grpstats2[0 .< grpstats2.condition .<= 2,:], "correct", "true", [1,2]; xlims=(-3, 9), xticks=([-2, 3, 5,7, 9], [-3, 3, 5,7, 9]), ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)
#savefig("WMM/Figures/model_based/expe1_$(model)_partial2_stable.pdf")
## Noise effect
grp_plot(grpstats2[1 .< grpstats2.condition .<= 3,:], "correct", "false", [2,3]; xlims=(-3, 9), xticks=([-2, 3, 5,7, 9], [-3,  3, 5,7, 9]), ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)
#savefig("WMM/Figures/model_based/expe1_$(model)_noise_correct.pdf")

## Noise effect (persev)
grp_plot(grpstats2[1 .< grpstats2.condition .<= 3,:], "explo", "false", [2,3]; xlims=(-3, 9), xticks=([-2,  3, 5,7, 9], [-3, 3, 5,7, 9]), ylims=(0, 0.3), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. non perseverative\nincorrect choice", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)
#savefig("WMM/Figures/model_based/expe1_$(model)_noise_explo.pdf")
## Noise effect (stable)
df_noise = copy(df_env2) # Start with a new df to relock the blocks on noise change (10 trials before rule change)
gp = groupby(df_noise, [:subject, :sessNum])
for gi in gp
    gi.blockNum .= vcat( gi.blockNum[10:end], fill(gi.blockNum[end], 9))
    gi.condition .= vcat( gi.condition[10:end], fill(gi.condition[end], 9))
    gi.nextCondition .= vcat( gi.nextCondition[10:end], fill(gi.nextCondition[end], 9))
end

trials_in_block!(df_noise)
pres_in_block!(df_noise)

for t = 1:nrow(df_noise)
    if df_noise[t, :stableAS] && (df_noise[t, :trialsInBlock] <= 10) && (df_noise.condition[t] > 0)
        df_noise.isStable[t] = true
    end
end
grpstats3 = grp_stats(df_noise)

##
grp_plot(grpstats3[(2 .< grpstats3.condition .<= 4) ,:], "correct", "true", [3,4]; linestyle=[ :solid :dot],  xlims=(1, 14), xticks=(2:2:12), ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", legend=:bottomright, size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)
grp_plot!(grpstats3[(grpstats3.condition .== 4) ,:], "correct", "false", [4]; linestyle=:solid, xlims=(1, 14), xticks=(2:2:12), ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", legend=:bottomright, size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)
#savefig("WMM/Figures/model_based/expe1_$(model)_noise_stable.pdf")


## Partial effect env1 VS env2 (correct)
gg = vcat(grpstats1[grpstats1.condition .== 3,: ], grpstats2[grpstats2.condition .== 1,: ])
env1idx = gg.condition .== 3
gg.condition[gg.condition .== 1] .= 2
gg.condition[env1idx] .= 1
grp_plot(gg, "correct", "false", [3,5]; xlims=(-3, 9),xticks=-2:2:9, ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300)

## HMM analysis
# Recompute pre-swtich perseveration
big_df[!,:correct_or_persev] = copy(big_df.persev)
pres_in_block!(big_df, blockId = :blockNum) 
big_df[big_df.negpresInBlock .> -3,:correct_or_persev] = big_df[big_df.negpresInBlock .> -3,:correct]

# Recode for stable stims 
big_df[big_df.stableAS, :correct_or_persev] = big_df[big_df.stableAS, :correct]

# Block count locked on HMM switch
big_df[!,:switchNum] = zeros(Int, nrow(big_df)) 
gdf = groupby(big_df, [:subject, :task, :sessNum])
for g in gdf
    idx = diff(vcat(0, findall(g.HMMSwitch .== 1) .- 1, nrow(g)))
    
    t = 0
    for b in eachindex(idx)
        g[t+1:t+idx[b], :switchNum] .= b
        t += idx[b]
    end 
end
pres_in_block!(big_df, blockId = :switchNum)
trials_in_block!(big_df, blockId = :switchNum)


## Perseveration after switches
## For task 1
df1 = big_df[big_df.task .== "WMM1",:]
df1.condition[df1.condition .== 4] .= 1 # Merge conditions 1 and 4
for t = 1:nrow(df1)
    if df1[t, Symbol("isStable_$(df1.stimulus[t])")] && df1.condition[t] > 0 
        df1.condition[t] = 4 # Make stable stims a special condition
    end
end

## Show perseveration locked on strategic changes
alldf = copy(df1)
grpstats = grp_stats_hmm(alldf)
grpstats = grpstats[0 .< grpstats.condition,: ] 
grp_plot_hmm(grpstats, "persev", [1,2,3,4]; xlims=(-2, 9), xticks=-2:2:9,ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. perseverative choice", legend_position = :topright, linestyle=[:solid :solid :solid :solid], background_color=:transparent, foreground_color=:black, labelfontsize=20, tickfontsize=14, dpi=300)
hline!([1/3], color=:black, linewidth=3, linestyle=:dash, label="", size=(500, 500), legendfontsize=14)


## Partial paired switches in condition 3
compdf = df1[df1.switchType .== 2,:]
grpstats = grp_stats_hmm(compdf)
grpstats = grpstats[grpstats.condition .== 3,: ] 
# Remove other condition for the pre-switch trials

grp_plot_hmm(grpstats, "persev", [3, 3]; xlims=(-2.1, 8.1), xticks=([-2.1, 0.1, 2.1, 4.1, 6.1, 8.1], [-3, 1, 3, 5, 7, 9]),ylims=(0, 1), yticks=0:0.2:1, label="Partial rule change\n(changing associations)", xlabel="Stimulus presentations", ylabel="Perseveration")
hline!([1/3], color=:black, linewidth=3, linestyle=:dash, label="Chance level", size=(500, 500))


## Global switches in condition 3 (partial) for stable stims : split between switch types
partdf = df1[df1.condition .== 4,:]
partdf.persev .= partdf.correct#partdf.correct_or_persev
partdf.condition[abs.(partdf.switchType) .< 3] .= 1 # Recode conditions as switch types
partdf.condition[abs.(partdf.switchType) .== 3] .= 2
grpstats2 = grp_stats_hmm(partdf)
grpstats2 = grpstats2[0 .< grpstats2.condition .<= 2,: ]
grp_plot_hmm(grpstats2, "persev", [4,4]; linestyle=[:solid :dot],xlims=(-2, 9), xticks=-2:2:9, ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300)
#hline!([1/3], color=:black, linewidth=3, linestyle=:dash, label="", tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)


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
        continue
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
                if !isempty(postIdx)
                    gd.afterExplo[postIdx[min(z, length(postIdx))]] = z
                end
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
            if nrow(gd) < tin+z || tin + z - 1 > tout# gd.ExploOut[tin+z-1] == 1# 
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
@df before_df[1 .> before_df.beforeExplo .> -4,:] plot(:beforeExplo, :persev_mean, ribbon=:persev_sem, label="", linewidth = 3, color=:black, xticks=(-2:2:exploWin), size=(500, 500), background_color=:transparent, foreground_color=:black, labelfontsize=20, tickfontsize = 14, dpi=300, ylabel="Prop. perseverative choice", xlabel="Stimulus presentations")
@df in_df[0 .< in_df.inExplo .<= exploWin,:] plot!(:inExplo,  :persev_mean, ribbon=:persev_sem, label="", linewidth = 3, color=:black, xticks=(-2:2:exploWin), ylims=(0, 1.0))
hline!([1/3], color=:black, linewidth=3, linestyle=:dash, label="", size=(500, 500), background_color=:transparent, foreground_color=:black, labelfontsize=20, tickfontsize = 14, dpi=300)

##

@df dur_df[(dur_df.presInBlock .<= exploWin) .* (0 .< dur_df.condition .<= 2),:] plot(:presInBlock, :ExploOut, group=:condition, linewidth = 5, palette=StatsPlots.palette(:Dark2)[1:2], alpha=1.0, label="", xticks=(2:2:exploWin), xlims=(0, 9), yticks=(0:0.1:0.3), ylims=(0, 0.45),labelfontsize=20, tickfontsize = 14, xlabel="Stimulus presentations", ylabel="Probability of switching out\nof the random strategy", background_color=:transparent, foreground_color=:black, size=(500, 500), dpi=300)
#hline!([1/3], color=:black, linewidth=3, linestyle=:dash, label="")



##
tmp = out_df[out_df.outExplo .>= -4,:]
sort!(tmp, :outExplo)
@df tmp plot(:outExplo, :cor_mean, ribbon=:cor_sem, label="", linewidth = 3, color=[:blue :black], ylims=(0, 1.0), xticks=(-6:2:6))
@df after_df[after_df.afterExplo .< 6,:] plot!(:afterExplo,:cor_mean, ribbon = :cor_sem, label="", linewidth = 3, color=[:blue :black], ylabel="Prop. correct choice", xlabel="Stimulus presentations", legendfontsize=12)
hline!([1/3], color=:black, linewidth=3, linestyle=:dash, label="", size=(500, 500), background_color=:transparent, foreground_color=:black, labelfontsize=20, tickfontsize = 14, dpi=300)


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
@df latsum violin(:condition, :latency_mean, group=:condition, xticks = [],xlabel="", ylims = (0,20), title="Random strategy duration\n(# presentations)",ylabel = "", label = "", palette=StatsPlots.palette(:Dark2), alpha=0.6, linewidth=0)
@df latsum dotplot!(:condition, :latency_mean, group=:condition, label = "", palette=StatsPlots.palette(:Dark2)[1:2], size=(700, 700), labelfontsize=30, tickfontsize=18, titlefontsize=30,background_color=:transparent, foreground_color=:black, xaxis=:off, dpi=300, markersize=5)

# add summary stats 
ss = combine(groupby(latsum, :condition), :latency_mean => mean => :latency_mean, :latency_mean => sem => :latency_sem)
scatter!(ss.condition, ss.latency_mean, yerror=ss.latency_sem, color=:black, markershape=:circle, markerstrokewidth=3, markersize=12, label="")


## Latency
## Env 1 
gdf = groupby(df1, [:subject, :blockNum, :sessNum])
latdf = DataFrame(subject = zeros(Int, length(gdf)), blockNum = zeros(Int, length(gdf)), condition = zeros(Int, length(gdf)), switchType = zeros(Int, length(gdf)), latency_trials = zeros(Int, length(gdf)), latency_pres = fill(NaN, length(gdf)))

for gi in 1:length(gdf)
    g = gdf[gi]
    latdf[gi, :subject] = g.subject[1]
    latdf[gi, :blockNum] = g.blockNum[1]
    latdf[gi, :condition] = g.condition[1] < 4 ? g.condition[1] : 3 # Reclassify stable associations in condition 3 g.condition[1]#
    latdf[gi, :switchType] = g.switchType[1]
    if g.blockNum[1] > 1
        idx = findfirst(g.HMMSwitch .== 1)
        if !isnothing(idx)
            latdf[gi, :latency_trials] = idx
            mean([g.stimulus[1:idx] .== stim for stim = 1:3])
            latdf[gi, :latency_pres] = mean([count(g.stimulus[1:idx] .== stim) for stim = 1:3]) # Averagre presentations from the begining of the block to the switch
        end
    end
end

##Recode conditions as condtiion X switchType
latdf = latdf[(latdf.condition .> 0) .* .!isnan.(latdf.latency_pres), :]


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

