using CSV, DataFrames, ProgressMeter, StatsPlots
include("main.jl")
include("preprocessing.jl")
include("stats_plots_funs.jl")

## Load simus
model = "cfql"
task = 1
tasknames = ["WMM1", "WMM2", "ALL"]

pathname = string("WMM/Simus/all_simus/Simus/", model, "_", tasknames[task])
fileslist = filter(x -> occursin(".csv", x), readdir(pathname))
big_df = DataFrame()

wb = Progress(length(fileslist), 1, "loading...")
for fi in eachindex(fileslist)
    fl = fileslist[fi]
    df = CSV.read(string(pathname, "/", fl), DataFrame)
    df.subject .+= fi * 1000
    df[!,:zrt] = (df.rt .- mean(df.rt)) ./ std(df.rt)
    pres_in_block!(df)
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
df_env1.condition[df_env1.condition .== 4] .= 1 # Merge conditions 1 and 4
df_env1[!, :isStable] = falses(nrow(df_env1)) # Re-code stable associations from rule change
for t = 1:nrow(df_env1)
    if df_env1[t, Symbol("isStable_$(df_env1.stimulus[t])")] && df_env1.condition[t] > 0
        df_env1.isStable[t] = true
    end
end

grpstats1 = grp_stats(df_env1);
# recidx = findall(grpstats1.condition .== 2)
# grpstats1.condition[grpstats1.condition .== 3] .= 2
# grpstats1.condition[recidx] .= 3

## Recurrence effect (correct)
grp_plot!(grpstats1[0 .< grpstats1.condition .<= 2,:], "correct", "false", [1,2], "bar"; xlims=(-3, 9), xticks=(-2:2:9), ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), tickfontsize=14, labelfontsize=20, background_color=:white, foreground_color=:black, dpi=300)
#savefig("WMM/Figures/model_based/expe1_$(model)_recurrence_correct.pdf")
## Recurrence effect (explo)
grp_plot!(grpstats1[0 .< grpstats1.condition .<= 2,:], "explo", "false", [1,2], "bar"; xlims=(-3, 9), xticks=(-2:2:9), ylims=(0,0.3), yticks=0:0.1:1, label="", xlabel="Stimulus presentations", ylabel="Non perseverative incorrect choice", size=(500, 500))
#savefig("WMM/Figures/model_based/expe1_$(model)_recurrence_explo.pdf")
##
grp_plot(grpstats1[grpstats1.condition .< 3,:], "persev", [1,2]; xlims=(-3, 12), ylims=(0, 1), yticks=0:0.2:1, label=["New rule" "Recurent rule"], xlabel="Stimulus presentations", ylabel="Perseveration", size=(500, 500))

## Partial effect (correct)
grp_plot!(grpstats1[(grpstats1.condition .== 1) .|(grpstats1.condition .== 3),:], "correct", "false", [1, 3], "bar"; linestyle=[:solid :solid], xlims=(-3, 9), ylims=(0, 1),xticks=-2:2:9, yticks=0:0.2:1, xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), label="", labelfontsize=20, tickfontsize=14,  background_color = :transparent, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black, dpi=300)
grp_plot!(grpstats1[grpstats1.condition .== 3,:], "correct", "true", [4], "bar"; linestyle=:dot, xlims=(-3, 9), ylims=(0, 1),xticks=-2:2:9, yticks=0:0.2:1, xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), label="", labelfontsize=20, tickfontsize=14,  background_color = :white, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black, dpi=300)

#savefig("WMM/Figures/model_based/expe1_$(model)_partial_stable.pdf")

## Partial vs complete (correct)
grp_plot(grpstats1[(grpstats1.condition .== 1) .| (grpstats1.condition .== 3),:], "correct", [1,3]; xlims=(-3, 9), ylims=(0, 1), yticks=0:0.2:1, label=["Complete rule change" "Partial rule change"], xlabel="Stimulus presentations", ylabel="Correct choice", size=(500, 500))
savefig("WMM/Figures/model_based/expe1_$(model)_partial_correct.pdf")
## Partial effect (explo)
grp_plot(grpstats1[(grpstats1.condition .== 1) .| (grpstats1.condition .== 3),:], "explo", [1,3]; xlims=(-3, 9), ylims=(0, 0.3), yticks=0:0.2:1, label=["New rule" "Partial changing"], xlabel="Stimulus presentations", ylabel="Perseveration", size=(500, 500))
savefig("WMM/Figures/model_based/expe1_$(model)_partial_explo.pdf")


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
    # eidx = g.ExploOut .== 1
    # for i in eachindex(eidx)
    #     if eidx[i] && g.switchType[i] ≠ -3
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
pres_in_block!(big_df, blockId = :switchNum)
trials_in_block!(big_df, blockId = :switchNum)

## Switch type per condition
## Env1
df1 = big_df[big_df.task .== "WMM1",:]
df1.condition[df1.condition .== 4] .= 1
stdf = df1[df1.newBlock .== 1,:]
gdf = groupby(stdf, [:subject, :condition])

stdf1 = combine(gdf, :switchType => (x -> mean(x .== -3)) => :random, :switchType => (x -> mean(x .== 3)) => :global, :switchType => (x -> mean(x .== 2)) => :paired, :switchType => (x -> mean(0 .< x .<= 2)) => :overlapping, :switchType => (x -> mean(abs.(x) .== 3)) => :nonoverlapping)
stdf1 = stdf1[stdf1.condition .> 0,:]
##
@df stdf1 violin(:condition, :overlapping, group=:condition, label="", ylims=(0, 1), ylabel="Proportion of switches", palette=StatsPlots.palette(:Dark2)[1:3], alpha=0.6, linewidth=0, size=(500, 500), tickfontsize=14, labelfontsize=18, background_color=:transparent, foreground_color=:black)
@df stdf1 dotplot!(:condition, :overlapping, group=:condition, label="", ylims=(0, 1), ylabel="Proportion of strategic changes", palette=StatsPlots.palette(:Dark2)[1:3])
violin!(stdf1.condition .+ 4, stdf1[:,:global], label="", group=stdf1.condition, alpha=0.6, linewidth=0)
dotplot!(stdf1.condition .+ 4, stdf1[:,:global], label="", group=stdf1.condition)
violin!(stdf1.condition .+ 8, stdf1[:,:random], label="", group=stdf1.condition, alpha=0.6, linewidth=0)
dotplot!(stdf1.condition .+ 8, stdf1[:,:random], label="", group=stdf1.condition)
xticks!([2, 6, 10], ["Overlapping", "Global", "Random"])
xlabel!("New strategy")


## add summary stats 
ss = combine(groupby(stdf1, :condition), :random => mean, :random => sem, :global => mean, :global => sem, :overlapping => mean, :overlapping => sem)
scatter!(ss.condition, ss.overlapping_mean, yerror=ss.overlapping_sem, color=:black, markershape=:circle, markerstrokewidth=3, markersize=6, label="")
scatter!(ss.condition .+ 4, ss.global_mean, yerror=ss.global_sem, color=:black, markershape=:circle, markerstrokewidth=3, markersize=6, label="")
scatter!(ss.condition .+ 8, ss.random_mean, yerror=ss.random_sem, color=:black, markershape=:circle, markerstrokewidth=3, markersize=6, label="")

xticks!([2, 6, 10], ["Overlapping", "Global", "Random"])
#xticks!([2, 6, 10], ["Global", "Overlapping"])
xlabel!("New strategy")

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


## Latency
## Env 1 
gdf = groupby(df1, [:subject, :blockNum, :sessNum])
latdf = DataFrame(subject = zeros(Int, length(gdf)), blockNum = zeros(Int, length(gdf)), condition = zeros(Int, length(gdf)), switchType = zeros(Int, length(gdf)), latency_trials = zeros(Int, length(gdf)), latency_pres = zeros(length(gdf)))

for gi in 1:length(gdf)
    g = gdf[gi]
    latdf[gi, :subject] = g.subject[1]
    latdf[gi, :blockNum] = g.blockNum[1]
    latdf[gi, :condition] = g.condition[1]#g.condition[1] < 4 ? g.condition[1] : 3 # Reclassify stable associations in condition 3
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
latdf = latdf[latdf.condition .> 0, :]
#latdf = latdf[latdf.switchType .== -3,:] # keep only global switches

sdf = groupby(latdf, [:subject, :condition, :switchType])
latsum = combine(sdf, :latency_pres => mean => :latency_mean)
latdf1 = latsum[latsum.condition .== 3,:]
latdf1.condition .= 1
@df latsum[latsum.switchType .== 3,:] violin(:condition, :latency_mean, group=:condition, label = "", ylabel = "Rule change - Strategic switch\nlatency (# presentations)", palette=StatsPlots.palette(:Dark2)[1:3], alpha=0.6, linewidth=0, size=(500, 600),  background_color = :transparent, foreground_color=:black, labelfontsize=20, tickfontsize=14)
@df latsum[latsum.switchType .== 3,:] dotplot!(:condition, :latency_mean, group=:condition, palette=StatsPlots.palette(:Dark2)[1:3], label="")
violin!(latsum[latsum.switchType .== 2,:condition] .+ 4, latsum[latsum.switchType .== 2,:latency_mean], group=latsum[latsum.switchType .== 2,:condition], alpha=0.6, linewidth=0, label = "")
dotplot!(latsum[latsum.switchType .== 2,:condition] .+ 4, latsum[latsum.switchType .== 2,:latency_mean], group=latsum[latsum.switchType .== 2,:condition], label = "")
violin!(latsum[latsum.switchType .== -3,:condition] .+ 8, latsum[latsum.switchType .== -3,:latency_mean], group=latsum[latsum.switchType .== -3,:condition], alpha=0.6, linewidth=0, label = "")
dotplot!(latsum[latsum.switchType .== -3,:condition] .+ 8, latsum[latsum.switchType .== -3,:latency_mean], group=latsum[latsum.switchType .== -3,:condition], label="")
xticks!([2, 6, 10], ["Global", "Overlapping", "Random"])
#xticks!([2, 6], ["Global", "Overlapping"])
xlabel!("Switch type")



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


## Effect of joint probability on probability of repeat 
df_env1[!,:MutualInfo] = zeros(nrow(df_env1))
gdf = groupby(df_env1, [:sessNum, :subject])
probaRep = zeros(11,2) # probability of repetition of a given action, given the repeition of this action in the 5 previous presentations of the same stim (6 possible values) and the repetition of the joint event {action for the current stim and the previous choice} in the 5 previous presentations for each stim (11 possible values)

normCount = zeros(11,2) # Normalization count
# PJ = zeros(2,2,length(gdf))
# PF = zeros(2,length(gdf))
# PC = zeros(2,length(gdf))
# MI = zeros(length(gdf))

win = 10
wb = Progress(length(gdf), 1, "processing...")
for gi = 1:length(gdf)
    subdf = DataFrame(gdf[gi])
    t1 = findfirst(subdf.blockNum .> 1)
    # PF[1,gi] = mean(subdf[t1,:fb])
    # PF[2,gi] = 1 - mean(subdf[t1,:fb])
    # PC[1,gi] = mean(subdf[t1,:correct])
    # PC[2,gi] = 1 - mean(subdf[t1,:correct])
    for t = 45:nrow(subdf)
        # PF = [mean(subdf[t-win-1:t-1,:fb]), 1-mean(subdf[t-win-1:t-1,:fb])]
        # PC = [mean(subdf[t-win:t,:correct]) 1-mean(subdf[t-win:t,:correct])]
        # PJ = [
        #     mean(subdf[t-win-1:t-1,:fb] .* subdf[t-win:t,:correct]) mean(subdf[t-win-1:t-1,:fb] .* (1 .- subdf[t-win:t,:correct]));
        #     mean((1 .- subdf[t-win-1:t-1,:fb]) .* subdf[t-win:t,:correct]) mean((1 .- subdf[t-win-1:t-1,:fb]) .* (1 .- subdf[t-win:t,:correct]))
        # ]
        # subdf[t,:MutualInfo] = sum(PJ .* log.(PJ ./ PF ./ PC))

        # Find all relevant indices (no stimulus repeat)
        s = subdf[t,:stimulus]
        sp = subdf[t-1,:stimulus]
        ap = subdf[t-1, :choice]
        fb = Int(subdf[t-1,:fb])
        if s ≠ sp
            # PF[1,gi] += subdf[t-1,:fb]
            # PF[2,gi] += 1 - subdf[t-1,:fb]
            # PC[1,gi] += subdf[t,:correct]
            # PC[2,gi] += 1 - subdf[t,:correct]
            # PJ[1,1,gi] += subdf[t-1,:fb] * subdf[t,:correct]
            # PJ[1,2,gi] += subdf[t-1,:fb] * (1 - subdf[t,:correct])
            # PJ[2,1,gi] += (1 - subdf[t-1,:fb]) * subdf[t,:correct]
            # PJ[2,2,gi] += (1 - subdf[t-1,:fb]) * (1 - subdf[t,:correct])


            for a = 1:3

                # Find the number of previous repetitions in the 5 last presentations
                prevIdx = findall(subdf[1:t-1, :stimulus] .== s)[end-4:end]
                priorRep = sum(subdf[prevIdx, :choice] .== a)

                # Find the number of joint repetitions 
                prevIdxP = findall(subdf[1:t-1, :stimulus] .== sp)[end-4:end]
                priorJoint = 0
                for i in prevIdx
                    idx = findlast(subdf[1:i, :stimulus] .== sp)
                    priorJoint += subdf[idx, :choice] == ap && subdf[i,:choice] == a
                end
                for i in prevIdxP
                    idx = findlast(subdf[1:i, :stimulus] .== s)
                    priorJoint += subdf[idx,:choice] == a && subdf[i,:choice] == ap
                end

                probaRep[priorJoint+1, fb+1] += (subdf[t,:choice] == a ) - priorRep/5
                normCount[priorJoint+1, fb+1]  += 1
            end
        end
    end
    # PF[:,gi] ./= length(t1:nrow(subdf))
    # PC[:,gi] ./= length(t1:nrow(subdf))
    # PJ[:,:,gi] ./= length(t1:nrow(subdf))
    # MI[gi] = sum(PJ[:,:,gi] .* log.(PJ[:,:,gi] ./ PF[:,gi] ./ PC[:,gi]'))
    next!(wb)
end
probaRep ./= normCount

##
heatmap(0:0.1:1.0, 0:0.2:1.0, probaRep[:,:,2])

## Effect of trap in counterfactual inference 
gdf = groupby(df_env1, [:sessNum, :subject])
N = length(gdf)
probaRep = zeros(8,N)
normCount = zeros(8,N)
wb = Progress(length(gdf), 1, "processing...")
for gi = 1:length(gdf)
    subdf = DataFrame(gdf[gi])
    for t in 45:nrow(subdf)
        # Find all relevant indices (no stimulus repeat)
        s = subdf[t,:stimulus]
        sp = subdf[t-1,:stimulus]
        a = subdf[t,:choice]
        trap = subdf[t-1,:fb] ≠ subdf[t-1,:correct]
        fb = subdf[t-1,:fb]
        prevIdx = findlast(subdf[1:t-1,:stimulus] .== s)
        prevFb = subdf[prevIdx, :fb]
        prevC = subdf[prevIdx, :correct]
        sfb = sum(subdf[prevIdx:t-1,:fb])
        nfb = length(prevIdx:t-1)
        # if prevFb == 0 && nfb < 6
        #     probaRep[sfb+1, nfb, gi] += subdf[prevIdx,:choice] ≠ a
        #     normCount[sfb+1, nfb, gi] += 1
        # elseif prevFb == 1 && nfb < 6
        #     probaRep[10+sfb+1, nfb, gi] += subdf[prevIdx,:choice] == a
        #     normCount[10+sfb+1, nfb, gi] += 1
        # end
        if s ≠ sp
            idx = prevC * 2 + fb + 1
            if prevFb == 0 
                probaRep[idx, gi] += subdf[prevIdx,:choice] ≠ a
                normCount[idx, gi] += 1 
            elseif prevFb == 1
                probaRep[4 + idx, gi] += subdf[prevIdx,:choice] == a
                normCount[4 + idx, gi] += 1
            end 
        end
    end
    next!(wb)
end

##
V = Matrix((probaRep ./ normCount)')
dotplot([V[:,2] .- V[:,1], V[:,4] .- V[:,3]], label="", xticks=([1,2], ["Lose-shift", "Win-stay"]),ylims=(-0.14, 0.22), ylabel="Proba. action repeat difference\npositive vs negative fb", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)