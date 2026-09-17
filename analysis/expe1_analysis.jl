using CSV, DataFrames
include("utils.jl")
include("preprocessing.jl")
include("stats_plots_funs.jl")
include("memory_analysis.jl")

## Load data
folder = "data/main_study"
flist = filter(x -> occursin(".csv", x), readdir(folder))

df = DataFrame()

for f in flist
    tmp = CSV.read(string(folder, "/", f), DataFrame)
    tmp = preprocess(tmp)
    tmp.sessNum .= parse(Int, f[7]) + 1
    mrt = mean(tmp.rt)
    srt = std(tmp.rt)
    tmp[!,:zrt] = (tmp.rt .- mrt) ./ srt
    append!(df, tmp)
end

# Session order 
gps = groupby(df, [:subject, :task])
for i = 1:length(gps)
    firstSess = gps[i].sessNum[1]
    secondSess = gps[i].sessNum[end]
    gps[i][!, :sessOrder] = 2 .- (gps[i].sessNum .== firstSess)
end

# Split by Environment
envdf = groupby(df, :task)

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
grp_plot(grpstats0[(grpstats0.condition .== 1) .| (grpstats0.condition .== 4),:], "correct", "false", [1,6]; xlims=(-3, 9), ylims=(0, 1), xticks=-2:2:9, yticks=0:0.2:1,  xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), label="", labelfontsize=20, tickfontsize=14,  background_color = :transparent, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black, dpi=300)

## Statistical Significance (cluster based permutation test)
tmp = df_env1[((df_env1.condition .== 1) .| (df_env1.condition .== 4)) .* (df_env1.presInBlock .<= 10) ,:] # after switch
tmp2 = df_env1[((df_env1.condition .== 1) .| (df_env1.condition .== 4)) .* (df_env1.negpresInBlock .>= -3) ,:] # before switch, trials that are both <10 trials post switch and <3trials pre switch will count twice
tmp2.presInBlock .= tmp2.negpresInBlock
tmp = vcat(tmp, tmp2)
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, [:correct, :condition] => ((x, y) -> mean(x[y .== 1]) - mean( x[y .== 4])) => :diff1_2)
X = zeros(length(unique(cc.subject)), 13)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:10
        jidx = gp[i].presInBlock[j]
        X[i,jidx] = gp[i].diff1_2[j]
    end
end
res = cluster_perm_test(X; niter=1e5)

## Compare new vs recurrent
df_env1.condition[df_env1.condition .== 4] .= 1 # Merge conditions 1 and 4
grpstats1 = grp_stats(df_env1)

## Learning (correct)
grp_plot(grpstats1[0 .< grpstats1.condition .<= 2,:], "correct", "false", [1,2]; xlims=(-3, 9), ylims=(0, 1), xticks=-2:2:9, yticks=0:0.2:1,  xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), label="", labelfontsize=20, tickfontsize=14,  background_color = :transparent, dpi=300)

## Statistical Significance (cluster based permutation test)
tmp = df_env1[(0 .< df_env1.condition .<= 2) .* (df_env1.presInBlock .<= 10) ,:] # after switch
tmp2 = df_env1[(0 .< df_env1.condition .<= 2) .* (df_env1.negpresInBlock .>= -3) ,:] # before switch, trials that are both <10 trials post switch and <3trials pre switch will count twice
tmp2.presInBlock .= tmp2.negpresInBlock
tmp = vcat(tmp, tmp2)
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, [:correct, :condition] => ((x, y) -> mean(x[y .== 1]) - mean( x[y .== 2])) => :diff1_2)
X = zeros(length(unique(cc.subject)), 13)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:10
        jidx = gp[i].presInBlock[j]
        X[i,jidx] = gp[i].diff1_2[j]
    end
end
res = cluster_perm_test(X; niter=1e5)
plot!(res.clusters[1], 0.95 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[2], alpha=0.5)


## Learning (explo)
grp_plot(grpstats1[0 .< grpstats1.condition .<= 2,:], "explo","false", [1,2]; xlims=(-3, 9), ylims=(0,0.3), xticks=-2:2:9, yticks=0:0.1:1, label="", xlabel="Stimulus presentations", ylabel="Prop. exploratory choice", size=(500, 500), labelfontsize=20, tickfontsize=14, background_color = :transparent, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black, dpi=300)

## Statistical Significance (cluster based permutation test)
tmp = df_env1[(0 .< df_env1.condition .<= 2) .* (df_env1.presInBlock .<= 10) ,:] # after switch
tmp2 = df_env1[(0 .< df_env1.condition .<= 2) .* (df_env1.negpresInBlock .>= -3) ,:] # before switch, trials that are both <10 trials post switch and <3trials pre switch will count twice
tmp2.presInBlock .= tmp2.negpresInBlock
tmp = vcat(tmp, tmp2)
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, [:explo, :condition] => ((x, y) -> mean(x[y .== 1]) - mean( x[y .== 2])) => :diff1_2)
X = zeros(length(unique(cc.subject)), 13)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    X[i,:] .= gp[i].diff1_2
end
res = cluster_perm_test(X; niter=1e5)
plot!(res.clusters[1], 0.285 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[2], alpha=0.5)


## Partial effect (correct)
grp_plot(grpstats1[(grpstats1.condition .== 1) .|(grpstats1.condition .== 3),:], "correct", "false", [1, 3]; linestyle=[:solid :solid], xlims=(-3, 9), ylims=(0, 1),xticks=-2:2:9, yticks=0:0.2:1, xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), label="", labelfontsize=20, tickfontsize=14,  background_color = :transparent, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black, dpi=300)
grp_plot!(grpstats1[grpstats1.condition .== 3,:], "correct", "true", [4]; linestyle=:dot, xlims=(-3, 9), ylims=(0, 1),xticks=-2:2:9, yticks=0:0.2:1, xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), label="", labelfontsize=20, tickfontsize=14,  background_color = :transparent, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black, dpi=300)

## Statistical Significance (cluster based permutation test)
# 1/ Stable effect (compared to pre switch performance)
M = combine(groupby(df_env1[df_env1.negtrialsInBlock .>= -10,:], :subject), :correct => mean)
tmp = df_env1[(df_env1.condition .== 3) .* (df_env1.isStable) .* (df_env1.presInBlock .<= 10),:]

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
plot!(res.clusters[1], 0.95 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[4], alpha=0.5)

## t-test on first presentation only :
# OneSampleTTest(X[:,3]) # p < 1e-3, t-stat = -3.81, df = 50
hline([0.0], color=:grey, linewidth=6, linestyle=:dash, alpha=0.5,label="", xlims=(0.4, 1.6), xticks=[], ylabel="Prop. correct choice", background_color=:transparent, size=(500, 500), dpi=300, foreground_color=:black, labelfontsize=32, tickfontsize=14)
bar!([1], [mean(X[:,1])], yerror = [sem(X[:,1])], linewidth=0, msw=5, color=StatsPlots.palette(:Dark2)[4], alpha=0.6, label="")
dotplot!([1], X[:,1], bar_width=0.3,color=:grey, msw=0, label="")
scatter!([1], [mean(X[:,1])], yerror = [sem(X[:,1])],color=StatsPlots.palette(:Dark2)[4], msc=:black, markershape=:circle, markerstrokewidth=5, markersize=0, label="")


## Direct comparison on successive trials for stable associations around rule change

firstStableIdx = findall((df_env1.stableAS) .* (df_env1.presInBlock .== 1) .* (df_env1.condition .== 3)) # find first presentation of stable stim after partial rule change

D = DataFrame(subject = df_env1.subject[firstStableIdx], tvt_1 = zeros(length(firstStableIdx)), t_1vt_2 = zeros(length(firstStableIdx)), fb_inbetween = zeros(length(firstStableIdx)))

for i in eachindex(firstStableIdx)
    idx = firstStableIdx[i]
    t_1idx = findlast(df_env1.stimulus[1:idx-1] .== df_env1.stimulus[idx]) # find index of last stimulus presentation before rule change
    t_2idx = findlast(df_env1.stimulus[1:t_1idx-1] .== df_env1.stimulus[idx]) # find index of penultimate stimulus presentation before rule change
    tp1idx = findfirst(df_env1.stimulus[idx+1:end] .== df_env1.stimulus[idx]) # find index of next stimulus presentation after rule change
    if !isnothing(t_1idx)
        D[i,:tvt_1] = df_env1.correct[idx] - df_env1.correct[t_1idx]
        if t_1idx == idx - 1 # if consecutive
            D[i,:fb_inbetween] = -1
        else
            D[i,:fb_inbetween] = sum(((df_env1.persev[t_1idx + 1 : idx - 1]) .* (.!df_env1.fb[t_1idx + 1 : idx - 1])))
        end
    else
        D[i,:tvt_1] = NaN
        D[i,:fb_inbetween] = NaN
    end

    if !isnothing(t_2idx)
        D[i,:t_1vt_2] = df_env1.correct[t_1idx] - df_env1.correct[t_2idx]
    else
        D[i,:t_1vt_2] = NaN
    end

end
# Per subject
gdf = combine(groupby(D, :subject), :tvt_1 => mean, :t_1vt_2 => mean)
bar([1 2], [mean(gdf.t_1vt_2_mean) mean(gdf.tvt_1_mean)], ylims=(-0.28, 0.22), color=[:white :orange],label=["Pre rule change" "Post vs Pre rule change"])
@df gdf dotplot!([1 2], [:t_1vt_2_mean, :tvt_1_mean], bar_width=0.3,color=:grey, msw=0,label="")
scatter!([1, 2], [mean(gdf.t_1vt_2_mean), mean(gdf.tvt_1_mean)], yerror = [sem(gdf.t_1vt_2_mean), sem(gdf.tvt_1_mean)], msc=:black, markershape=:circle, markerstrokewidth=5, markersize=0, label="", xticks=[], ylabel="Performance difference", size=(500, 500), dpi=300, labelfontsize=18, legendfontsize=12, tickfontsize=12, legend_position=:topright)
##
gdf = combine(groupby(D, [:subject, :fb_inbetween]), :tvt_1 => mean, :t_1vt_2 => mean)
sdf = combine(groupby(gdf, :fb_inbetween), :tvt_1_mean => mean => :tvt_1_mean, :tvt_1_mean => sem => :tvt_1_sem, :t_1vt_2_mean => mean => :t_1vt_2_mean, :t_1vt_2_mean => sem => :t_1vt_2_sem)

ggdf = groupby(gdf, :fb_inbetween)
bar(1:5, [mean(ggdf[i].tvt_1_mean) for i = 1:5], color=:orange, ylims=(-0.25, 0.05), xticks=(1:5, ["Consecutive", 0, 1, 2, 3]), ylabel = "Performance difference", xlabel="# Disconfirmatory feedback\nbetween presentations",label="",size=(500, 500), dpi=300, labelfontsize=18, legendfontsize=12, tickfontsize=12)

scatter!(1:5, [mean(ggdf[i].tvt_1_mean) for i = 1:5], yerror = [sem(ggdf[i].tvt_1_mean) for i = 1:5], msc=:black, markershape=:circle, markerstrokewidth=5, markersize=0, label="")



# @df sdf[1:5,:] plot(:fb_inbetween, :t_1vt_2_mean, ribbon=:t_1vt_2_sem, label="Pre rule change", ylabel = "Performance difference", xlabel="# Disconfirmatory feedback in between", color=:black, linewidth=3)
# @df sdf[1:5,:] plot!(:fb_inbetween, :tvt_1_mean, ribbon=:tvt_1_sem, label="Post vs Pre rule change", ylabel = "Performance difference", xlabel="# Disconfirmatory feedback\nbetween presentations", xticks=(-1:3, ["Consecutive", 0, 1, 2, 3]), color=:orange, linewidth=3, size=(500, 500), dpi=300, labelfontsize=18, legendfontsize=12, tickfontsize=12, legend_position=:bottomleft)




## 2/ Partial learning effect (condition 1 VS 2)
tmp = df_env1[(0 .< df_env1.condition .<= 3) .* (df_env1.presInBlock .<= 10) .* .!df_env1.isStable ,:] 
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, [:correct, :condition] => ((x, y) -> mean(x[y .== 1]) - mean( x[y .== 3])) => :diff1_2)
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


## Check inter-session effects 
df_sess1 = df_env1[df_env1.sessOrder .== 1,:]
df_sess2 = df_env1[df_env1.sessOrder .== 2,:]
grpstats_sess1 = grp_stats(df_sess1)
grpstats_sess2 = grp_stats(df_sess2)

# First session
# grp_plot(grpstats_sess1[0 .< grpstats_sess1.condition .<= 2,:], "correct", "false", [1,2];xlims=(-3, 9), ylims=(0, 1), xticks=-2:2:9, yticks=0:0.2:1,  xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), label="", title = "First session", titlefontsize=20, labelfontsize=20, tickfontsize=14,  background_color = :transparent, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black, dpi=300)
grp_plot(grpstats_sess1[(grpstats_sess1.condition .== 1) .|(grpstats_sess1.condition .== 3),:], "correct", "false", [1, 3]; linestyle=[:solid :solid], xlims=(-3, 9), ylims=(0, 1),xticks=-2:2:9, yticks=0:0.2:1, xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), label="", labelfontsize=20, tickfontsize=14,  background_color = :transparent, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black, dpi=300)
grp_plot!(grpstats_sess1[grpstats_sess1.condition .== 3,:], "correct", "true", [4]; linestyle=:dot, xlims=(-3, 9), ylims=(0, 1),xticks=-2:2:9, yticks=0:0.2:1, xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), label="", labelfontsize=20, tickfontsize=14,  background_color = :transparent, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black, dpi=300)

## Statistical significance : recurrence
tmp = df_sess1[(0 .< df_sess1.condition .<= 2) .* (df_sess1.presInBlock .<= 10) ,:] # 

gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, [:correct, :condition] => ((x, y) -> mean(x[y .== 1]) - mean( x[y .== 2])) => :diff1_2)
X = zeros(length(unique(cc.subject)), 13)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:10
        jidx = gp[i].presInBlock[j]
        X[i,jidx] = gp[i].diff1_2[j]
    end
end
res = cluster_perm_test(X; niter=1e5)
plot!(res.clusters[1], 0.95 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[2], alpha=0.5)

## Stats : partial
# 1/ Stable effect (compared to pre switch performance)
M = combine(groupby(df_sess1[df_sess1.negtrialsInBlock .>= -10,:], :subject), :correct => mean)
tmp = df_sess1[(df_sess1.condition .== 3) .* (df_sess1.isStable) .* (df_sess1.presInBlock .<= 10),:]

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
plot!(res.clusters[1], 0.95 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[4], alpha=0.5)

## t-test on first presentation only :
hline([0.0], color=:grey, linewidth=6, linestyle=:dash, alpha=0.5,label="", xlims=(0.4, 1.6), xticks=[], ylabel="Prop. correct choice", background_color=:transparent, size=(500, 500), dpi=300, foreground_color=:black, labelfontsize=32, tickfontsize=14)
bar!([1], [mean(X[:,1])], yerror = [sem(X[:,1])], linewidth=0, msw=5, color=StatsPlots.palette(:Dark2)[4], alpha=0.6, label="")
dotplot!([1], X[:,1], bar_width=0.3,color=:grey, msw=0, label="")
scatter!([1], [mean(X[:,1])], yerror = [sem(X[:,1])],color=StatsPlots.palette(:Dark2)[4], msc=:black, markershape=:circle, markerstrokewidth=5, markersize=0, label="")


## 2/ Partial learning effect (condition 1 VS 2)
tmp = df_sess1[(0 .< df_sess1.condition .<= 3) .* (df_sess1.presInBlock .<= 10) .* .!df_sess1.isStable ,:] 
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, [:correct, :condition] => ((x, y) -> mean(x[y .== 1]) - mean( x[y .== 3])) => :diff1_2)
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

## Second session
# grp_plot(grpstats_sess2[0 .< grpstats_sess2.condition .<= 2,:], "correct", "false", [1,2]; xlims=(-3, 9), ylims=(0, 1), xticks=-2:2:9, yticks=0:0.2:1,  xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), label="", title = "Second session", titlefontsize=20, labelfontsize=20, tickfontsize=14,  background_color = :transparent, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black, dpi=300)

grp_plot(grpstats_sess2[(grpstats_sess2.condition .== 1) .|(grpstats_sess2.condition .== 3),:], "correct", "false", [1, 3]; linestyle=[:solid :solid], xlims=(-3, 9), ylims=(0, 1),xticks=-2:2:9, yticks=0:0.2:1, xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), label="", labelfontsize=20, tickfontsize=14,  background_color = :transparent, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black, dpi=300)
grp_plot!(grpstats_sess2[grpstats_sess2.condition .== 3,:], "correct", "true", [4]; linestyle=:dot, xlims=(-3, 9), ylims=(0, 1),xticks=-2:2:9, yticks=0:0.2:1, xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), label="", labelfontsize=20, tickfontsize=14,  background_color = :transparent, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black, dpi=300)

## Stats : recurrence
tmp = df_sess2[(0 .< df_sess2.condition .<= 2) .* (df_sess2.presInBlock .<= 10) ,:] # 

gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, [:correct, :condition] => ((x, y) -> mean(x[y .== 1]) - mean( x[y .== 2])) => :diff1_2)
X = zeros(length(unique(cc.subject)), 13)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:10
        jidx = gp[i].presInBlock[j]
        X[i,jidx] = gp[i].diff1_2[j]
    end
end
res = cluster_perm_test(X; niter=1e5)
plot!(res.clusters[1], 0.95 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[2], alpha=0.5)


## Stats : partial
# 1/ Stable effect (compared to pre switch performance)
M = combine(groupby(df_sess2[df_sess2.negtrialsInBlock .>= -10,:], :subject), :correct => mean)
tmp = df_sess2[(df_sess2.condition .== 3) .* (df_sess2.isStable) .* (df_sess2.presInBlock .<= 10),:]

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
plot!(res.clusters[1], 0.95 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[4], alpha=0.5)

## t-test on first presentation only :
hline([0.0], color=:grey, linewidth=6, linestyle=:dash, alpha=0.5,label="", xlims=(0.4, 1.6), xticks=[], ylabel="Prop. correct choice", background_color=:transparent, size=(500, 500), dpi=300, foreground_color=:black, labelfontsize=32, tickfontsize=14)
bar!([1], [mean(X[:,1])], yerror = [sem(X[:,1])], linewidth=0, msw=5, color=StatsPlots.palette(:Dark2)[4], alpha=0.6, label="")
dotplot!([1], X[:,1], bar_width=0.3,color=:grey, msw=0, label="")
scatter!([1], [mean(X[:,1])], yerror = [sem(X[:,1])],color=StatsPlots.palette(:Dark2)[4], msc=:black, markershape=:circle, markerstrokewidth=5, markersize=0, label="")


## 2/ Partial learning effect (condition 1 VS 2)
tmp = df_sess2[(0 .< df_sess2.condition .<= 3) .* (df_sess2.presInBlock .<= 10) .* .!df_sess2.isStable ,:] 
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, [:correct, :condition] => ((x, y) -> mean(x[y .== 1]) - mean( x[y .== 3])) => :diff1_2)
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

## Memory effect at the strategy level 
df_mem = copy(df_env1)
group_level_mem!(df_mem, 1)

validIdx = df_mem.blockNum .> 1
# Keep only transitions from and to non recurrent rules 
for i in findall(validIdx)
    prevBlock = df_mem.blockNum[i] - 1
    stim = df_mem.stimulus[i]
    prevIdx = findfirst((df_mem.blockNum .== prevBlock) .* (df_mem.stimulus .== stim))
    validIdx[i] = (df_mem.condition[prevIdx] ≠ 2) * (df_mem.condition[i] == 1) 
end

gdf = groupby(df_mem[validIdx ,:], [:subject, :evForRec, :presInBlock, :consistent])
summary_sub = combine(gdf, :congruentDiff => mean => :congruentDiff)

gdf = groupby(summary_sub, [:evForRec, :presInBlock, :consistent])
summary_mem = combine(gdf, :congruentDiff => mean, :congruentDiff => sem)

summary_mem = summary_mem[.!isnan.(summary_mem.congruentDiff_sem),:]

@df summary_mem[summary_mem.consistent .== 1,:] plot(:presInBlock, :congruentDiff_mean, ribbon=:congruentDiff_sem, group=:evForRec, linewidth=3, xlim=(0,12),legend_title="Feedback", labels=["Negative" "Positive"], xlabel="Position of feedback in episode\n(presentation # from rule change)", ylabel="Δ choice as recurrent rule for other stims", size=(500, 500), dpi=300, background_color=:transparent)

hline!([0.0], label="", color=:black, linewidth=3, linestyle=:dash)

## Statistical Significance (cluster based permutation test) 

uniqueSubs = unique(summary_sub.subject)
X = fill(NaN, length(uniqueSubs), 12)
for i in eachindex(uniqueSubs)
    for j = 1:12
        sPos = findfirst((summary_sub.subject .== uniqueSubs[i]) .* (summary_sub.presInBlock .== j) .* (summary_sub.evForRec .== 1.0) .* (summary_sub.consistent .== 0.0))
        sNeg = findfirst((summary_sub.subject .== uniqueSubs[i]) .* (summary_sub.presInBlock .== j) .* (summary_sub.evForRec .== 0.0) .* (summary_sub.consistent .== 0.0))
        if !isnothing(sPos) && !isnothing(sNeg)
            X[i, j] = summary_sub.congruentDiff[sPos] - summary_sub.congruentDiff[sNeg]
        end
    end
end

res = cluster_perm_test(X; niter=1e5)
# cluster = [2], pval = 0.0325

plot!([-0.5, 0.5] .+ res.clusters[1], 0.1 .* ones(length(res.clusters[1]) + 1), linewidth=5, label="", color = StatsPlots.palette(:auto)[1], alpha=0.5, ylims=(-0.13, 0.105)) 

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
grp_plot(grpstats2[0 .< grpstats2.condition .<= 2,:], "correct","false", [5,2]; xlims=(-3, 9),xticks=(-2:2:9), ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300)

## Statistical Significance (cluster based permutation test)
tmp = df_env2[(0 .< df_env2.condition .<= 2) .* (df_env2.presInBlock .<= 10) .* .!df_env2.isStable,:] # after switch
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, [:correct, :condition] => ((x, y) -> mean(x[y .== 1]) - mean( x[y .== 2])) => :diff1_2)
X = zeros(length(unique(cc.subject)), 10)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:10
        jidx = findfirst(gp[i].presInBlock .== j)
        X[i,jidx] = gp[i].diff1_2[j]
    end
end
res = cluster_perm_test(X; niter=1e5)


## Relearning (Explo)
grp_plot(grpstats2[0 .< grpstats2.condition .<= 2,:], "explo","false", [5,2]; xlims=(-3, 9),xticks=(-2:2:9), ylims=(0, 0.3), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. non perseverative\nincorrect choice", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300)

## Statistical Significance (cluster based permutation test)
tmp = df_env2[(0 .< df_env2.condition .<= 2) .* (df_env2.presInBlock .<= 10) .* .!df_env2.isStable,:] # after switch
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, [:explo, :condition] => ((x, y) -> mean(x[y .== 1]) - mean( x[y .== 2])) => :diff1_2)
X = zeros(length(unique(cc.subject)), 10)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:10
        jidx = gp[i].presInBlock[j]
        X[i,jidx] = gp[i].diff1_2[j]
    end
end
res = cluster_perm_test(X; niter=1e5)
plot!(res.clusters[1], 0.25 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[2], alpha=0.5)


## Relearning (Persev)
grp_plot(grpstats2[0 .< grpstats2.condition .<= 2,:], "persev","false", [5,2]; xlims=(-3, 9),xticks=(-2:2:9), ylims=(0, 1.0), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. perseverative choice", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300)

## Statistical Significance (cluster based permutation test)
tmp = df_env2[(0 .< df_env2.condition .<= 2) .* (df_env2.presInBlock .<= 10) .* .!df_env2.isStable,:] # after switch
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, [:persev, :condition] => ((x, y) -> mean(x[y .== 1]) - mean( x[y .== 2])) => :diff1_2)
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

## Partial effect (stable)
grp_plot(grpstats2[0 .< grpstats2.condition .<= 2,:], "correct", "true", [5,2]; xlims=(-3, 9), xticks=(-2:2:9), ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300)

## Statistical Significance (cluster based permutation test)
## 1/ Difference between conditions
tmp = df_env2[(0 .< df_env2.condition .<= 2) .* (df_env2.isStable) .* (df_env2.presInBlock .<= 10) ,:] # after switch
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
plot!(res.clusters[1], 0.97 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[3], alpha=0.5)

## 2/ Difference from plateau condition 1
M = combine(groupby(df_env2[(0 .< df_env2.nextCondition .<= 2) .* (df_env2.negtrialsInBlock .> -10),:], :subject), :correct => mean)
tmp = df_env2[(df_env2.condition .== 1) .* (df_env2.isStable) .* (df_env2.presInBlock .<= 15),:]
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, :correct => mean)
X = fill(NaN, length(unique(cc.subject)), 15)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:nrow(gp[i])
        jidx = findfirst(gp[i].presInBlock .== j)
        X[i,jidx] = gp[i].correct_mean[j] - M.correct_mean[i]
    end
end
res = cluster_perm_test(X; niter=1e5)
plot!(0:9, fill(mean(M.correct_mean), 10), color = :grey, linewidth=3, linestyle=:dash, label="", alpha=0.5)
plot!(res.clusters[1], 0.95 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[5], alpha=0.5)

## t-test on first presentation only :
OneSampleTTest(X[:,1]) # p = 0.1224, t-stat = 1.57, df = 50
hline([0.0], color=:grey, linewidth=6, linestyle=:dash, alpha=0.5,label="", xlims=(0.4, 1.6), xticks=[], ylabel="Prop. correct choice", background_color=:transparent, size=(500, 500), dpi=300, foreground_color=:black, labelfontsize=32, tickfontsize=14)
violin!([1], X[:,1], alpha=0.6, linewidth=0, color=StatsPlots.palette(:Dark2)[4], label="")
dotplot!([1], X[:,1], color=StatsPlots.palette(:Dark2)[4], label="")
scatter!([1], [mean(X[:,1])], yerror = [sem(X[:,1])],color=:black, markershape=:circle, markerstrokewidth=5, markersize=12, label="")

## 3/ Difference from plateau condition 2
M = combine(groupby(df_env2[(0 .< df_env2.nextCondition .<= 2) .* (df_env2.negtrialsInBlock .> -10),:], :subject), :correct => mean)
tmp = df_env2[(df_env2.condition .== 2) .* (df_env2.isStable) .* (df_env2.presInBlock .<= 15),:]
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, :correct => mean)
X = fill(NaN, length(unique(cc.subject)), 15)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:nrow(gp[i])
        jidx = gp[i].presInBlock[j]
        X[i,jidx] = gp[i].correct_mean[j] - M.correct_mean[i]
    end
end
res = cluster_perm_test(X; niter=1e5)
plot!(res.clusters[1], 0.9 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[2], alpha=0.5)

## t-test on first presentation only :
OneSampleTTest(X[:,1]) # p = 0.16, t-stat = 1.43, df = 50
hline([0.0], color=:grey, linewidth=6, linestyle=:dash, alpha=0.5,label="", xlims=(0.4, 1.6), xticks=[], ylabel="Prop. correct choice", background_color=:transparent, size=(500, 500), dpi=300, foreground_color=:black, labelfontsize=32, tickfontsize=14)
violin!([1], X[:,1], alpha=0.6, linewidth=0, color=StatsPlots.palette(:Dark2)[4], label="")
dotplot!([1], X[:,1], color=StatsPlots.palette(:Dark2)[4], label="")
scatter!([1], [mean(X[:,1])], yerror = [sem(X[:,1])],color=:black, markershape=:circle, markerstrokewidth=5, markersize=12, label="")

## Noise effect
grp_plot(grpstats2[1 .< grpstats2.condition .<= 3,:], "correct", "false", [2,3]; xlims=(-3, 9), xticks=(-2:2:9), ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300)

## Statistical Significance (cluster based permutation test)
tmp = df_env2[(1 .< df_env2.condition .<= 3) .* (df_env2.presInBlock .<= 10) .* .!df_env2.isStable,:] # after switch

gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, [:correct, :condition] => ((x, y) -> mean(x[y .== 3]) - mean( x[y .== 2])) => :diff1_2)
X = fill(NaN, length(unique(cc.subject)), 10)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:10
        jidx = findfirst(gp[i].presInBlock .== j)
        X[i,jidx] = gp[i].diff1_2[jidx]
    end
end
res = cluster_perm_test(X; niter=1e5)


## Noise effect (explo)
grp_plot(grpstats2[1 .< grpstats2.condition .<= 3,:], "explo", "false", [2,3]; xlims=(-3, 9), xticks=(-2:2:12),ylims=(0, 0.3), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. non perseverative\nincorrect choice", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)

## Statistical Significance (cluster based permutation test)
tmp = df_env2[(1 .< df_env2.condition .<= 3) .* (df_env2.presInBlock .<= 10) .* .!df_env2.isStable ,:] # after switch
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, [:explo, :condition] => ((x, y) -> mean(x[y .== 3]) - mean( x[y .== 2])) => :diff1_2)
X = fill(NaN, length(unique(cc.subject)), 10)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:10
        jidx = gp[i].presInBlock[j]
        X[i,jidx] = gp[i].diff1_2[j]
    end
end
res = cluster_perm_test(X; niter=1e5)
plot!(res.clusters[1], 0.95 * 0.3 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[2], alpha=0.5)


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
grp_plot!(grpstats3[(grpstats3.condition .== 4) ,:], "correct", "false", [4]; linestyle=:solid, xlims=(1, 14), xticks=(2:2:12), ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", legend=:bottomright, size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300)


## Statistical Significance (cluster based permutation test)
## 1/ Difference between conditions
tmp = df_noise[(2 .< df_noise.condition .<= 4) .* (df_noise.isStable) .* (df_noise.presInBlock .<= 15) ,:] # after switch
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, [:correct, :condition] => ((x, y) -> mean(x[y .== 3]) - mean( x[y .== 4])) => :diff1_2)
X = zeros(length(unique(cc.subject)), 15)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:15
        jidx = gp[i].presInBlock[j]
        X[i,jidx] = gp[i].diff1_2[j]
    end
end
res = cluster_perm_test(X; niter=1e5)
plot!(res.clusters[2], 0.97 .* ones(length(res.clusters[2])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[3], alpha=0.5)

## 2/ Difference from plateau condition 3 (switch + noise)
M = combine(groupby(df_env2[((0 .< df_env2.condition .<= 2) .| (5 .< df_env2.condition .<= 7)) .* (df_env2.negtrialsInBlock .> -10),:], :subject), :correct => mean)
tmp = df_env2[(df_env2.condition .== 8) .* (df_env2.presInBlock .<= 15),:]
tmp2 = df_env2[(df_env2.nextCondition .== 3) .* (df_env2.stableAS) .* (df_env2.negpresInBlock .> -3) ,:] # before switch, trials that are both <10 trials post switch and <3trials pre switch will count twice
tmp2.presInBlock .= tmp2.negpresInBlock .+ 1
tmp = vcat(tmp, tmp2)
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, :correct => mean)
X = fill(NaN, length(unique(cc.subject)), 17)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:nrow(gp[i])
        jidx = gp[i].presInBlock[j]
        X[i,jidx+2] = M.correct_mean[i] - gp[i].correct_mean[j]
    end
end
res = cluster_perm_test(X; niter=1e5)
plot!(0:9, fill(mean(M.correct_mean), 10), color = :grey, linewidth=3, linestyle=:dash, label="", alpha=0.5)
plot!(res.clusters[1] .- 2, 0.95 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[3], alpha=0.5)


## 3/ Difference from plateau condition 2
M = combine(groupby(df_env2[((0 .< df_env2.condition .<= 2) .| (5 .< df_env2.condition .<= 7)) .* (df_env2.negtrialsInBlock .> -10),:], :subject), :correct => mean)
tmp = df_env2[(df_env2.condition .== 9) .* (df_env2.presInBlock .<= 15),:]
tmp2 = df_env2[(df_env2.nextCondition .== 4) .* (df_env2.stableAS) .* (df_env2.negpresInBlock .> -3) ,:] # before switch, trials that are both <10 trials post switch and <3trials pre switch will count twice
tmp2.presInBlock .= tmp2.negpresInBlock .+ 1
tmp = vcat(tmp, tmp2)
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, :correct => mean)
X = fill(NaN, length(unique(cc.subject)), 17)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:nrow(gp[i])
        jidx = gp[i].presInBlock[j]
        X[i,jidx+2] = M.correct_mean[i] - gp[i].correct_mean[j]
    end
end
res = cluster_perm_test(X; niter=1e5)
plot!(0:9, fill(mean(M.correct_mean), 10), color = :grey, linewidth=3, linestyle=:dash, label="", alpha=0.5)
plot!(res.clusters[1] .- 2, 0.93 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[4], alpha=0.5)

## 2/ Difference from plateau condition 3 (switch + noise)
M = combine(groupby(df_env2[((0 .< df_env2.condition .<= 2) .| (5 .< df_env2.condition .<= 7)) .* (df_env2.negtrialsInBlock .> -10),:], :subject), :correct => mean)
tmp = df_env2[(df_env2.condition .== 4) .* (df_env2.presInBlock .<= 15),:]
tmp2 = df_env2[(df_env2.nextCondition .== 4) .* (df_env2.stableAS) .* (df_env2.negpresInBlock .> -3) ,:] # before switch, trials that are both <10 trials post switch and <3trials pre switch will count twice
tmp2.presInBlock .= tmp2.negpresInBlock .+ 1
tmp = vcat(tmp, tmp2)
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, :correct => mean)
X = fill(NaN, length(unique(cc.subject)), 17)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:nrow(gp[i])
        jidx = gp[i].presInBlock[j]
        X[i,jidx+2] = M.correct_mean[i] - gp[i].correct_mean[j]
    end
end
res = cluster_perm_test(X; niter=1e5)
plot!(0:9, fill(mean(M.correct_mean), 10), color = :grey, linewidth=3, linestyle=:dash, label="", alpha=0.5)
plot!(res.clusters[1] .- 2, 0.93 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[6], alpha=0.5)


## Partial effect env1 VS env2 (correct)
gg = vcat(grpstats1[grpstats1.condition .== 3,: ], grpstats2[grpstats2.condition .== 1,: ])
env1idx = gg.condition .== 3
gg.condition[gg.condition .== 1] .= 2
gg.condition[env1idx] .= 1
grp_plot(gg, "correct", "false", [3,5]; xlims=(-3, 9),xticks=-2:2:9, ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300)

## Statistical Significance (cluster based permutation test)
tmp = vcat(df_env1[(df_env1.condition .== 3) .* (.!df_env1.isStable),:], df_env2[(df_env2.condition .== 1) .* (.!df_env2.isStable),:])
tmp = tmp[tmp.presInBlock .<= 10,:]
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, [:correct, :condition] => ((x, y) -> mean(x[y .== 1]) - mean( x[y .== 3])) => :diff1_2)
cc = cc[cc.presInBlock .<= 10,:]
X = zeros(length(unique(cc.subject)), 10)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:10
        jidx = gp[i].presInBlock[j]
        X[i,jidx] = gp[i].diff1_2[j]
    end
end
res = cluster_perm_test(X; niter=1e5)
pl = Plots.current()
for i in eachindex(res.clusters)
    plot!(vcat(res.clusters[i][1]-0.25, res.clusters[i], res.clusters[i][1]+0.25), 0.95 .* ones(length(res.clusters[i])+2), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[5], alpha=0.5)
end
plot(pl)

## Check inter-session effects : check if the effect exists in the first session of the second environment (3rd session) - equivalent to removing the 4th session
df1_no4 = df_env1[df_env1.sessNum .< 4,:]
df2_no4 = df_env2[df_env2.sessNum .< 4,:]
grpstats_env1 = grp_stats(df1_no4)
grpstats_env2 = grp_stats(df2_no4)

gg = vcat(grpstats_env1[grpstats_env1.condition .== 3,: ], grpstats_env2[grpstats_env2.condition .== 1,: ])
env1idx = gg.condition .== 3
gg.condition[gg.condition .== 1] .= 2
gg.condition[env1idx] .= 1
grp_plot(gg, "correct", "false", [3,5]; xlims=(-3, 9),xticks=-2:2:9, ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300)

## Statistical Significance (cluster based permutation test)
tmp = vcat(df1_no4[(df1_no4.condition .== 3) .* (.!df1_no4.isStable),:], df2_no4[(df2_no4.condition .== 1) .* (.!df2_no4.isStable),:])
tmp = tmp[tmp.presInBlock .<= 10,:]
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, [:correct, :condition] => ((x, y) -> mean(x[y .== 1]) - mean( x[y .== 3])) => :diff1_2)
cc = cc[cc.presInBlock .<= 10,:]
X = zeros(length(unique(cc.subject)), 10)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:10
        jidx = gp[i].presInBlock[j]
        X[i,jidx] = gp[i].diff1_2[j]
    end
end
res = cluster_perm_test(X; niter=1e5)
pl = Plots.current()
for i in eachindex(res.clusters)
    if res.pvalues[i] <= 0.05
        plot!(vcat(res.clusters[i][1]-0.25, res.clusters[i], res.clusters[i][1]+0.25), 0.95 .* ones(length(res.clusters[i])+2), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[5], alpha=0.5)
    end
end
plot(pl)

## Partial effect env1 VS env2 (explo)
grp_plot(gg, "explo", "false", [3,1]; xlims=(-3, 9), xticks=([-2,0,1, 3, 5,7, 9], [-3,-1,1,  3, 5,7, 9]), ylims=(0, 0.3), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. non perseverative\nincorrect choice", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)

## Partial effect env1 VS env2 (stable)
gg = vcat(grpstats1[grpstats1.condition .== 3,: ], grpstats2[grpstats2.condition .== 1,: ])
grp_plot(gg, "correct", "true", [5,3]; xlims=(-3, 9), xticks=-2:2:9, ylims=(0, 1), yticks=0:0.2:1, linestyle=:dot, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300)

## Difference between env 1 and env 2
## Statistical Significance (cluster based permutation test)
tmp = vcat(df_env1[(df_env1.condition .== 3) .* (df_env1.isStable),:], df_env2[(df_env2.condition .== 1) .* (df_env2.isStable),:])
tmp = tmp[tmp.presInBlock .<= 10,:]
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, [:correct, :condition] => ((x, y) -> mean(x[y .== 3]) - mean( x[y .== 1])) => :diff1_2)
cc = cc[cc.presInBlock .<= 10,:]
X = zeros(length(unique(cc.subject)), 10)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:10
        jidx = gp[i].presInBlock[j]
        X[i,jidx] = gp[i].diff1_2[j]
    end
end
res = cluster_perm_test(X; niter=1e5)

plot!(res.clusters[1], 0.95 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[5], alpha=0.5)

## 2/ Difference from plateau env 2
M = combine(groupby(df_env2[(0 .< df_env2.nextCondition .<= 2) .* (df_env2.negtrialsInBlock .> -10),:], :subject), :correct => mean)
tmp = df_env2[(df_env2.condition .== 1) .* (df_env2.isStable) .* (df_env2.presInBlock .<= 15),:]

gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, :correct => mean)
X = fill(NaN, length(unique(cc.subject)), 15)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:nrow(gp[i])
        jidx = gp[i].presInBlock[j]
        X[i,jidx] = gp[i].correct_mean[j] - M.correct_mean[i]
    end
end
res = cluster_perm_test(X; niter=1e5)
plot!(0:9, fill(mean(M.correct_mean), 10), color = StatsPlots.palette(:Dark2)[5], linewidth=3, linestyle=:dash, label="", alpha=0.5)
plot!(res.clusters[1], 0.95 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[5], alpha=0.5)
X_env2 = X

## 3/ Difference from plateau env 1
M = combine(groupby(df_env1[df_env1.negtrialsInBlock .>= -10,:], :subject), :correct => mean)
tmp = df_env1[(df_env1.condition .== 3) .* (df_env1.isStable) .* (df_env1.presInBlock .<= 10),:]

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
plot!(0:9, fill(mean(M.correct_mean), 10), color = StatsPlots.palette(:Dark2)[3], linewidth=3, linestyle=:dash, label="", alpha=0.5)
plot!(res.clusters[1], 0.9 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[3], alpha=0.5)
X_env1 = X

## t-test on first presentation only :
# Env 1 : 
bar([1], [mean(X_env1[:,1])],color=StatsPlots.palette(:Dark2)[3], alpha=0.6, linewidth=0, label="", xlims=(0.4, 2.6), xticks=[], ylabel="Prop. correct choice", background_color=:transparent, size=(500, 500), dpi=300, foreground_color=:black, labelfontsize=32, tickfontsize=14)
dotplot!([1], X_env1[:,1], color=:grey, msw=0.0, bar_width=0.3,label="")
scatter!([1], [mean(X_env1[:,1])], yerror = [sem(X_env1[:,1])],color=:black, markershape=:circle, markerstrokewidth=5, markersize=0, label="")

## Env 2 : 

bar!([2], [mean(X_env2[:,1])],color=StatsPlots.palette(:Dark2)[5], alpha=0.6, linewidth=0, label="", xlims=(0.4, 2.6), xticks=[], ylabel="Prop. correct choice", background_color=:transparent, size=(500, 500), dpi=300, foreground_color=:black, labelfontsize=32, tickfontsize=14)
dotplot!([2], X_env2[:,1], color=:grey, msw=0.0, bar_width=0.3,label="")
scatter!([2], [mean(X_env2[:,1])], yerror = [sem(X_env1[:,1])],color=:black, markershape=:circle, markerstrokewidth=5, markersize=0, label="")

## Same but removing session 4 (checking if the effect exists already in the first session of the second environment)
gg = vcat(grpstats_env1[grpstats_env1.condition .== 3,: ], grpstats_env2[grpstats_env2.condition .== 1,: ])
grp_plot(gg, "correct", "true", [5,3]; xlims=(-3, 9), xticks=-2:2:9, ylims=(0, 1), yticks=0:0.2:1, linestyle=:dot, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300)

## Difference between env 1 and env 2
## Statistical Significance (cluster based permutation test)
tmp = vcat(df1_no4[(df1_no4.condition .== 3) .* (df1_no4.isStable),:], df2_no4[(df2_no4.condition .== 1) .* (df2_no4.isStable),:])
tmp = tmp[tmp.presInBlock .<= 10,:]
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, [:correct, :condition] => ((x, y) -> mean(x[y .== 3]) - mean( x[y .== 1])) => :diff1_2)
cc = cc[cc.presInBlock .<= 10,:]
X = zeros(length(unique(cc.subject)), 10)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:10
        jidx = gp[i].presInBlock[j]
        X[i,jidx] = gp[i].diff1_2[j]
    end
end
res = cluster_perm_test(X; niter=1e5)

plot!(res.clusters[1], 0.95 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[5], alpha=0.5)

## 2/ Difference from plateau env 2
M = combine(groupby(df2_no4[(0 .< df2_no4.nextCondition .<= 2) .* (df2_no4.negtrialsInBlock .> -10),:], :subject), :correct => mean)
tmp = df2_no4[(df2_no4.condition .== 1) .* (df2_no4.isStable) .* (df2_no4.presInBlock .<= 15),:]

gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, :correct => mean)
X = fill(NaN, length(unique(cc.subject)), 15)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:nrow(gp[i])
        jidx = gp[i].presInBlock[j]
        X[i,jidx] = gp[i].correct_mean[j] - M.correct_mean[i]
    end
end
res = cluster_perm_test(X; niter=1e5)
plot!(0:9, fill(mean(M.correct_mean), 10), color = StatsPlots.palette(:Dark2)[5], linewidth=3, linestyle=:dash, label="", alpha=0.5)
plot!(res.clusters[1], 0.95 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[5], alpha=0.5)
X_env2 = X

## 3/ Difference from plateau env 1
M = combine(groupby(df1_no4[df1_no4.negtrialsInBlock .>= -10,:], :subject), :correct => mean)
tmp = df1_no4[(df1_no4.condition .== 3) .* (df1_no4.isStable) .* (df1_no4.presInBlock .<= 10),:]

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
plot!(0:9, fill(mean(M.correct_mean), 10), color = StatsPlots.palette(:Dark2)[3], linewidth=3, linestyle=:dash, label="", alpha=0.5)
plot!(res.clusters[1], 0.9 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[3], alpha=0.5)
X_env1 = X

## t-test on first presentation only :
# Env 1 : 
bar([1], [mean(X_env1[:,1])],color=StatsPlots.palette(:Dark2)[3], alpha=0.6, linewidth=0, label="", xlims=(0.4, 2.6), xticks=[], ylabel="Prop. correct choice", background_color=:transparent, size=(500, 500), dpi=300, foreground_color=:black, labelfontsize=32, tickfontsize=14)
dotplot!([1], X_env1[:,1], color=:grey, msw=0.0, bar_width=0.3,label="")
scatter!([1], [mean(X_env1[:,1])], yerror = [sem(X_env1[:,1])],color=:black, markershape=:circle, markerstrokewidth=5, markersize=0, label="")

## Env 2 : 

bar!([2], [mean(X_env2[:,1])],color=StatsPlots.palette(:Dark2)[5], alpha=0.6, linewidth=0, label="", xlims=(0.4, 2.6), xticks=[], ylabel="Prop. correct choice", background_color=:transparent, size=(500, 500), dpi=300, foreground_color=:black, labelfontsize=32, tickfontsize=14)
dotplot!([2], X_env2[:,1], color=:grey, msw=0.0, bar_width=0.3,label="")
scatter!([2], [mean(X_env2[:,1])], yerror = [sem(X_env1[:,1])],color=:black, markershape=:circle, markerstrokewidth=5, markersize=0, label="")

