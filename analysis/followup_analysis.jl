using CSV, DataFrames, StatsPlots, StatsFuns
include("preprocessing.jl")
include("stats_plots_funs.jl")
include("/Users/sami/PhD/Model_Tasks_Data/Data/WMM/pipeline/utils.jl")

## Load data
folder = "WMM/RAW_FU"
flist = filter(x -> occursin(".csv", x), readdir(folder))

df = DataFrame()

for f in flist
    tmp = CSV.read(string(folder, "/", f), DataFrame)
    tmp = preprocess(tmp)
    mrt = movstat(mean, tmp.rt, 12)
    srt = movstat(std, tmp.rt, 12)
    tmp[!,:zrt] = (tmp.rt .- mrt) ./ srt
    append!(df, tmp)
end

# Split by Environment
envdf = groupby(df, :task)

## Keep only subjects that completed all 3 Environment
keep_subs = intersect([en.subject for en in envdf]...)

## Environment 1 
df_env1 = DataFrame(envdf[[en.task[1] == "task1" for en in envdf]])
df_env1[!, :isStable] = falses(nrow(df_env1))
grpstats1 = grp_stats(df_env1)
## Recurence effect (correct)
grp_plot(grpstats1[grpstats1.condition .> 0,:], "correct", "false", [1,2]; xlims=(-3, 12), xticks=(-2:2:12), ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500),tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300)

## Statistical Significance (cluster based permutation test)
tmp = df_env1[(0 .< df_env1.condition) .* (df_env1.presInBlock .<= 10) ,:] # after switch
#tmp2 = df_env1[(0 .< df_env1.condition) .* (df_env1.negpresInBlock .>= -3) ,:] # before switch, trials that are both <10 trials post switch and <3trials pre switch will count twice
#tmp2.presInBlock .= tmp2.negpresInBlock
#tmp = vcat(tmp, tmp2)
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
plot!(res.clusters[1], 0.95 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[2], alpha=0.5)

#savefig("WMM/Figures/model_free/expe2_recurrence_correct.png")
## Recurence effect (explo)
grp_plot(grpstats1[grpstats1.condition .> 0,:], "explo", "false", [1,2]; xlims=(-3, 9), xticks=([-2,0,1,  3, 5,7, 9], [-3,-1,1, 3, 5,7, 9]), ylims=(0, 0.3), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. non perseverative incorrect choice", size=(500, 500),tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)
#savefig("WMM/Figures/model_free/expe2_recurrence_explo.pdf")
## Environment 2 
df_env2 = DataFrame(envdf[[en.task[1] == "task2" for en in envdf]])
df_env2[!, :isStable] = falses(nrow(df_env2)) # Re-code stable associations from rule change
for t = 1:nrow(df_env2)
    if df_env2[t, Symbol("isStable_$(df_env2.stimulus[t])")] && df_env2.condition[t] > 0
        df_env2.isStable[t] = true
    end
end

grpstats2 = grp_stats(df_env2)
## Complete/Partial effect (correct)
grp_plot(grpstats2[1 .< grpstats2.condition .<= 2,:], "correct", "false", [3]; xlims=(-3, 9), xticks=(-2:2:10), ylims=(0, 1), yticks=0:0.2:1, linestyle=[:solid :solid :dot], label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", legend_position = :bottomright, size=(500, 500),tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300)

grp_plot!(grpstats2[grpstats2.condition .== 2,:], "correct", "true", [4]; xlims=(-3, 9), xticks=(-2:2:10), ylims=(0, 1), yticks=0:0.2:1, linestyle=:dot, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", legend_position = :bottomright, size=(500, 500),tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)

## Statistical Significance (cluster based permutation test)
# 1/ Stable effect (compared to pre switch performance)
M = combine(groupby(df_env2[(df_env2.negpresInBlock .>= -3),:], :subject), :correct => mean)
tmp = df_env2[(df_env2.condition .== 2) .* (df_env2.isStable) .* (df_env2.presInBlock .<= 10),:]
#tmp2 = df_env2[(df_env2.condition .== 2) .* (df_env2.stableAS) .* (df_env2.negpresInBlock .> -3) ,:] # before switch, trials that are both <10 trials post switch and <3trials pre switch will count twice
#tmp2.presInBlock .= tmp2.negpresInBlock .+ 1
#tmp = vcat(tmp, tmp2)
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, :correct => mean)
X = zeros(length(unique(cc.subject)), 10)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:10
        jidx = gp[i].presInBlock[j]
        X[i,jidx] = gp[i].correct_mean[j] - M.correct_mean[i]
    end
end
res = cluster_perm_test(X; niter=1e5)
plot!(0:9, fill(mean(M.correct_mean), 10), color = :grey, linewidth=3, linestyle=:dash, label="", alpha=0.5)
plot!(res.clusters[1] .- 2, 0.95 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[4], alpha=0.5)

## t-test on first presentation only :
OneSampleTTest(X[:,1])
hline([0.0], color=:grey, linewidth=6, linestyle=:dash, alpha=0.5,label="", xlims=(0.4, 1.6), xticks=[], ylabel="Prop. correct choice", background_color=:transparent, size=(500, 500), dpi=300, foreground_color=:black, labelfontsize=32, tickfontsize=14)
violin!([1], X[:,1], alpha=0.6, linewidth=0, color=StatsPlots.palette(:Dark2)[4], label="")
dotplot!([1], X[:,1], color=StatsPlots.palette(:Dark2)[4], label="")
scatter!([1], [mean(X[:,1])], yerror = [sem(X[:,1])],color=:black, markershape=:circle, markerstrokewidth=5, markersize=12, label="")

## 2/ Partial learning effect (condition 1 VS 2)
tmp = df_env2[(0 .< df_env2.condition .<= 2) .* (df_env2.presInBlock .<= 10) .* .!df_env2.isStable,:] # after switch
#tmp2 = df_env2[(0 .< df_env2.condition .<= 2) .* (df_env2.negpresInBlock .>= -3) ,:] # before switch, trials that are both <10 trials post switch and <3trials pre switch will count twice
#tmp2.presInBlock .= tmp2.negpresInBlock
#tmp = vcat(tmp, tmp2)
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
plot!(res.clusters[1], 0.9 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[3], alpha=0.5)
#savefig("WMM/Figures/model_free/expe2_partial_correct.png")
## Complete/Partial effect (persev)
grp_plot(grpstats2[0 .< grpstats2.condition .< 3,:], "explo", "false", [1, 3]; xlims=(-3, 9), xticks=([-2, 0, 1, 3, 5,7, 9], [-3, -1, 1, 3, 5,7, 9]), ylims=(0, 0.3), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Non perseverative incorrect choice", size=(500, 500))
#savefig("WMM/Figures/model_free/expe2_partial_explo.pdf")
## Complete/Partial effect (stable)
grp_plot(grpstats2[grpstats2.condition .>= 2,:], "correct", [3, 4];linestyle=[:solid :dot], xlims=(-3, 9), ylims=(0, 1), yticks=0:0.2:1, label=["Changing associations" "Stable associations"], xlabel="Stimulus presentations", ylabel="Correct choice", size=(500, 500))
#savefig("WMM/Figures/model_free/expe2_partial_stable.pdf")

## Environment 3 
df_env3 = DataFrame(envdf[[en.task[1] == "task3" for en in envdf]])

df_env3[!, :isStable] = falses(nrow(df_env3)) # Re-code stable associations from rule change
# for t = 1:nrow(df_env3)
#     if df_env3[t, Symbol("isStable_$(df_env3.stimulus[t])")] && df_env3.condition[t] > 0
#         df_env3.isStable[t] = true
#     end
# end
# grpstats3 = grp_stats(df_env3)
##
df_chg = deepcopy(df_env3)
df_chg[!,:blockStim] = ones(Int, nrow(df_chg))

# Split per stim
perstim = groupby(df_chg, [:subject, :stimulus])
blidx = 1
for ps in perstim
    stim = ps.stimulus[1]
    # For each stim, find the first changing stim of each block
    V = falses(nrow(ps))
    tmp = DataFrame(ps)
    gdf = groupby(tmp, :blockNum)
    for g in gdf
        idx = DataFrames.rows(g)[1] # first trial index in the whole dataframe
        V[idx] = !ps[idx,Symbol("isStable_$stim")]
    end

    blen = findall(V) .- 1 # block lengths
    for bl in 1:length(blen)-1
        ps[blen[bl] + 1:blen[bl+1],:blockStim] .= blidx # Block number for the corresponding stim
        blidx += 1
    end
    ps[(blen[end] + 1):nrow(ps),:blockStim] .= blidx
    blidx += 1
end
nextCondition!(df_chg)
df_chg.blockNum = df_chg.blockStim
grpstats3 = grp_stats(df_chg)
## Relearning (correct)
grp_plot(grpstats3[0 .< grpstats3.condition ,:], "correct", "false", [1,2]; xlims=(-3, 9), ylims=(0, 1), yticks=0:0.2:1, label=["2/3 rule change" "1/3 rule change"], xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500),tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)
# savefig("WMM/Figures/model_free/expe2_relearn_correct.pdf")
## Relearning (persev)
grp_plot(grpstats3[0 .< grpstats3.condition ,:], "explo","false", [1,2]; xlims=(-3, 9), ylims=(0, 0.3), yticks=0:0.2:1, label=["2/3 rule change" "1/3 rule change"], xlabel="Stimulus presentations", ylabel="Non perseverative incorrect choice", size=(500, 500))
#savefig("WMM/Figures/model_free/expe2_relearn_explo.pdf")


## Each condition 
grp_plot(grpstats3[0 .< grpstats3.condition .<= 1 ,:], "correct", "false", [1]; xlims=(0.5, 9), ylims=(0, 1), yticks=0:0.2:1, label="2/3 rule change", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500),tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)

## Stable stims
df_stable = df_env3[df_env3.stableAS .* (2 .>= df_env3.condition .> 0),:]#df_env3[df_env3.stableAS,:]#
df_stable.isStable .= true
grpstats4 = grp_stats(df_stable)
grp_plot(grpstats4[0 .< grpstats4.condition .<= 2,:], "correct", "true", [4,6]; xlims=(-3, 9), xticks=(-2:2:9), ylims=(0, 1), yticks=0:0.2:1, linestyle=:dot, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300)

## Statistical Significance (cluster based permutation test)
# Stable effect (compared to pre switch performance)
M = combine(groupby(df_stable[(df_stable.negpresInBlock .>= -3),:], :subject), :correct => mean)
tmp = df_stable[(df_stable.presInBlock .<= 10) .* (df_stable.condition .== 1),:]
#tmp2 = df_env2[(0 .< df_env2.condition .== 3) .* (df_env2.negpresInBlock .> -3) ,:] # before switch, trials that are both <10 trials post switch and <3trials pre switch will count twice
#tmp2.presInBlock .= tmp2.negpresInBlock .+ 1
#tmp = vcat(tmp, tmp2)
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, :correct => mean)
X = zeros(length(unique(cc.subject)), 10)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:nrow(gp[i])
        jidx = gp[i].presInBlock[j]
        X[i,jidx] = gp[i].correct_mean[j] - M.correct_mean[i]
    end
end
res = cluster_perm_test(X; niter=1e5)
#plot!(0:9, fill(mean(M.correct_mean), 10), color = :grey, linewidth=3, linestyle=:dash, label="", alpha=0.5)
plot!(res.clusters[3] .+ [-0.2, 0.2], 0.95 .* ones(length(res.clusters[3]) + 1), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[4], alpha=0.5)

## t-test on first presentation only :
SignedRankTest(X[:,1]) # p = 0.002
hline([0.0], color=:grey, linewidth=6, linestyle=:dash, alpha=0.5,label="", xlims=(0.4, 1.6), xticks=[], ylabel="Prop. correct choice", background_color=:transparent, size=(500, 500), dpi=300, foreground_color=:black, labelfontsize=32, tickfontsize=14)
violin!([1], X[:,1], alpha=0.6, linewidth=0, color=StatsPlots.palette(:Dark2)[4], label="")
dotplot!([1], X[:,1], color=StatsPlots.palette(:Dark2)[4], label="")
scatter!([1], [mean(X[:,1])], yerror = [sem(X[:,1])],color=:black, markershape=:circle, markerstrokewidth=5, markersize=12, label="")


#savefig("WMM/Figures/model_free/expe2_partial2_stable.pdf")
## Env1 VS Env2 on complete non recurent
tmp = grpstats2[grpstats2.condition .== 1,:]
tmp.condition .= 2
gg = vcat(grpstats1[grpstats1.condition .== 1,:], tmp)
## Correct
grp_plot(gg, "correct", [1,2]; xlims=(-3, 12), ylims=(0, 1), yticks=0:0.2:1, label=["New rule in env1" "New rule in env2"], xlabel="Stimulus presentations", ylabel="Correct choice", size=(500, 500))
# savefig("WMM/Figures/model_free/expe2_relearn1v2_correct.pdf")
## Persev
grp_plot(gg, "persev", [1,2]; xlims=(-3, 12), ylims=(0, 1), yticks=0:0.2:1, label=["New rule in env1" "New rule in env2"], xlabel="Stimulus presentations", ylabel="Perseveration", size=(500, 500))
savefig("WMM/Figures/model_free/expe2_relearn1v2_persev.pdf")

## Env2 VS Env3 on partial 1
grp_plot(grpstats2[1 .< grpstats2.condition .<= 2,:], "correct", "false", [ 3]; xlims=(-3, 9), xticks=(-2:2:10), ylims=(0, 1), yticks=0:0.2:1, linestyle=[:solid :solid :dot], label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", legend_position = :bottomright, size=(500, 500),tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)
grp_plot!(grpstats3[0 .< grpstats3.condition .<= 1 ,:], "correct", "false", [1]; xlims=(-3, 9), ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500),tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)

## Stable
grp_plot(grpstats2[grpstats2.condition .== 2,:], "correct", "true", [4]; xlims=(-3, 9), xticks=(-2:2:10), ylims=(0, 1), yticks=0:0.2:1, linestyle=:dot, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", legend_position = :bottomright, size=(500, 500),tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)

grp_plot!(grpstats4[0 .< grpstats4.condition .<= 1,:], "correct", "true", [5]; xlims=(-3, 9), ylims=(0, 1), yticks=0:0.2:1,linestyle=:dot, label="", size=(500, 500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)

##
tmp = df_chg[df_chg.condition .== 1,:]
select!(tmp, Not(:blockStim))
tmp.isStable .= false 
tmp = vcat(df_env2[df_env2.condition .== 2,:], tmp)

gg2 = grp_stats(tmp)
#gg2 = vcat(grpstats2[grpstats2.condition .== 1,:], tmp)
## Correct
grp_plot(gg2[gg2.condition .> 0,:], "correct", "false", [1,3]; xlims=(-3, 9), xticks=([-2,0,1, 3, 5,7, 9], [-3,-1,1, 3, 5,7, 9]), ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500),tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)
# savefig("WMM/Figures/model_free/expe2_partial1v3_correct.pdf")
## Persev
grp_plot(gg2, "explo", [3,1]; xlims=(-3, 9), ylims=(0, 0.3), yticks=0:0.2:1, label=["2/3 rule change in env. A3" "2/3 rule change in env. B2"], xlabel="Stimulus presentations", ylabel="Perseveration", size=(500, 500))
#savefig("WMM/Figures/model_free/expe2_partial1v3_persev.pdf")
## Stable
gg3 = vcat(grpstats2[grpstats2.condition .== 3,:], grpstats4[grpstats4.condition .== 1,:])
grp_plot(gg3, "correct", [1,3]; xlims=(-3, 9), xticks=([-2,0,1, 3, 5,7, 9], [-3,-1,1,  3, 5,7, 9]), ylims=(0, 1), yticks=0:0.2:1, linestyle=[:dot :dot], label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500),tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)
# savefig("WMM/Figures/model_free/expe2_partial1v3_stable.pdf")

## Comparision interference effect between task 2 and 3 :
hline([0.0], color=:grey, linewidth=6, linestyle=:dash, alpha=0.5,label="", xlims=(0.4, 3.6), xticks=[], xaxis=false, ylabel="Prop. correct choice", background_color=:transparent, size=(650, 500), dpi=300, foreground_color=:black, labelfontsize=32, tickfontsize=14)
violin!([1], Y1[:,1], alpha=0.6, linewidth=0, color=StatsPlots.palette(:Dark2)[4], label="")
dotplot!([1], Y1[:,1], color=StatsPlots.palette(:Dark2)[4], label="")
scatter!([1], [mean(Y1[:,1])], yerror = [sem(Y1[:,1])],color=:black, markershape=:circle, markerstrokewidth=5, markersize=12, label="")

violin!([2], Y2[:,1], alpha=0.6, linewidth=0, color=StatsPlots.palette(:Dark2)[4], label="")
dotplot!([2], Y2[:,1], color=StatsPlots.palette(:Dark2)[4], label="")
scatter!([2], [mean(Y2[:,1])], yerror = [sem(Y2[:,1])],color=:black, markershape=:circle, markerstrokewidth=5, markersize=12, label="")

violin!([3], Y3[:,1], alpha=0.6, linewidth=0, color=StatsPlots.palette(:Dark2)[6], label="")
dotplot!([3], Y3[:,1], color=StatsPlots.palette(:Dark2)[6], label="")
scatter!([3], [mean(Y3[:,1])], yerror = [sem(Y3[:,1])],color=:black, markershape=:circle, markerstrokewidth=5, markersize=12, label="")

## Comparision interference effect between conditions in task 3 :
hline([0.0], color=:grey, linewidth=6, linestyle=:dash, alpha=0.5,label="", xlims=(0.4, 2.6), xticks=[], xaxis=false, ylims=(-0.28, 0.18), yticks=(-0.2:0.1:0.1), ylabel="Prop. correct choice", background_color=:transparent, size=(500, 500), dpi=300, foreground_color=:black, labelfontsize=32, tickfontsize=14)
violin!([1], Y2[:,1], alpha=0.6, linewidth=0, color=StatsPlots.palette(:Dark2)[4], label="")
dotplot!([1], Y2[:,1], color=StatsPlots.palette(:Dark2)[4], label="")
scatter!([1], [mean(Y1[:,1])], yerror = [sem(Y1[:,1])],color=:black, markershape=:circle, markerstrokewidth=5, markersize=12, label="")

violin!([2], Y3[:,1], alpha=0.6, linewidth=0, color=StatsPlots.palette(:Dark2)[6], label="")
dotplot!([2], Y3[:,1], color=StatsPlots.palette(:Dark2)[6], label="")
scatter!([2], [mean(Y2[:,1])], yerror = [sem(Y2[:,1])],color=:black, markershape=:circle, markerstrokewidth=5, markersize=12, label="")


## Mutual information task 1 VS task 2 VS task 3
function mut_inf(X)
    pj = zeros(2,2)
    pj[1,1] = mean((X[1:end-1] .== 1) .* (X[2:end] .== 1))
    pj[1,2] = mean((X[1:end-1] .== 1) .* (X[2:end] .== 0))
    pj[2,1] = mean((X[1:end-1] .== 0) .* (X[2:end] .== 1))
    pj[2,2] = mean((X[1:end-1] .== 0) .* (X[2:end] .== 0))
    v = [pj[i,j] * log(pj[i,j] / sum(pj[i,:]) / sum(pj[:,j])) for i = 1:2, j = 1:2]
    mi = sum([pj[i,j] * log(pj[i,j] / sum(pj[i,:]) / sum(pj[:,j])) for i = 1:2, j = 1:2])
    return mi, v[:]
end

sdf = groupby(df_env1, [:subject, :sessNum])

MI_env1 = zeros(length(sdf))
PJ_env1 = zeros(4, length(sdf))
for i = 1:length(sdf)
    MI_env1[i], PJ_env1[:,i] = mut_inf(sdf[i].fb)
end

sdf = groupby(df_env2, [:subject, :sessNum])

MI_env2 = zeros(length(sdf))
PJ_env2 = zeros(4, length(sdf))
for i = 1:length(sdf)
    MI_env2[i], PJ_env2[:,i] = mut_inf(sdf[i].fb)
end

sdf = groupby(df_env3, [:subject, :sessNum])

MI_env3 = zeros(length(sdf))
PJ_env3 = zeros(4, length(sdf))
for i = 1:length(sdf)
    MI_env3[i], PJ_env3[:,i] = mut_inf(sdf[i].fb)
end

boxplot([MI_env1, MI_env2, MI_env3])

## Mutual information per stimulus 


sdf = groupby(df_env1, [:subject, :sessNum, :stimulus])

MI_env1_ps = zeros(length(sdf))
PJ_env1_ps = zeros(4, length(sdf))
for i = 1:length(sdf)
    MI_env1_ps[i], PJ_env1_ps[:,i] = mut_inf(sdf[i].fb)
end

sdf = groupby(df_env2, [:subject, :sessNum, :stimulus])

MI_env2_ps = zeros(length(sdf))
PJ_env2_ps = zeros(4, length(sdf))
for i = 1:length(sdf)
    MI_env2_ps[i], PJ_env2_ps[:,i] = mut_inf(sdf[i].fb)
end

sdf = groupby(df_env3, [:subject, :sessNum, :stimulus])

MI_env3_ps = zeros(length(sdf))
PJ_env3_ps = zeros(4, length(sdf))
for i = 1:length(sdf)
    MI_env3_ps[i], PJ_env3_ps[:,i] = mut_inf(sdf[i].fb)
end

boxplot([MI_env1_ps, MI_env2_ps, MI_env3_ps])