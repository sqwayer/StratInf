include("preprocessing.jl")
include("stats_plots_funs.jl")


## Load data
folder = "WMM/RAW2/complete"
flist = filter(x -> occursin(".csv", x), readdir(folder))

df = DataFrame()

for f in flist
    tmp = CSV.read(string(folder, "/", f), DataFrame)
    tmp = preprocess(tmp)
    tmp.sessNum .= parse(Int, f[7]) + 1
    mrt = movstat(mean, tmp.rt, 12)
    srt = movstat(std, tmp.rt, 12)
    tmp[!,:zrt] = (tmp.rt .- mrt) ./ srt
    append!(df, tmp)
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

#savefig("WMM/Figures/model_free/expe1_env1_learning_correct.png")
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

#savefig("WMM/Figures/model_free/expe1_env1_learning_explo.png")
## Learning (perseveration)
grp_plot(grpstats1[0 .< grpstats1.condition .<= 2,:], "persev", [1,2]; xlims=(-3, 9), ylims=(0, 1), xticks=([-2,0,1,  3, 5,7, 9], [-3,-1,1, 3, 5,7, 9]), yticks=0:0.2:1,  xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), label="", labelfontsize=20, tickfontsize=14,  background_color = :transparent, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black)


## Partial effect (correct)
# grp_plot(grpstats1[(grpstats1.condition .== 1) .|(grpstats1.condition .== 2) .| (grpstats1.condition .== 4),:], "correct", [1, 3,4]; linestyle=[:solid :solid :dot], xlims=(-3, 9), ylims=(0, 1),xticks=-2:2:9, yticks=0:0.2:1, xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), label="", labelfontsize=20, tickfontsize=14,  background_color = :transparent, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black, dpi=300)
grp_plot(grpstats1[(grpstats1.condition .== 1) .|(grpstats1.condition .== 3),:], "correct", "false", [1, 3]; linestyle=[:solid :solid], xlims=(-3, 9), ylims=(0, 1),xticks=-2:2:9, yticks=0:0.2:1, xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), label="", labelfontsize=20, tickfontsize=14,  background_color = :transparent, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black, dpi=300)
grp_plot!(grpstats1[grpstats1.condition .== 3,:], "correct", "true", [4]; linestyle=:dot, xlims=(-3, 9), ylims=(0, 1),xticks=-2:2:9, yticks=0:0.2:1, xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), label="", labelfontsize=20, tickfontsize=14,  background_color = :transparent, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black, dpi=300)

## Statistical Significance (cluster based permutation test)
# 1/ Stable effect (compared to pre switch performance)
M = combine(groupby(df_env1[df_env1.negtrialsInBlock .>= -10,:], :subject), :correct => mean)
tmp = df_env1[(df_env1.condition .== 3) .* (df_env1.isStable) .* (df_env1.presInBlock .<= 10),:]
#tmp2 = df_env1[(0 .< df_env1.nextCondition .== 3) .* (df_env1.negpresInBlock .> -3) ,:] # before switch, trials that are both <10 trials post switch and <3trials pre switch will count twice
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
plot!(res.clusters[1], 0.95 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[4], alpha=0.5)

## t-test on first presentation only :
# OneSampleTTest(X[:,3]) # p < 1e-3, t-stat = -3.81, df = 50
hline([0.0], color=:grey, linewidth=6, linestyle=:dash, alpha=0.5,label="", xlims=(0.4, 1.6), xticks=[], ylabel="Prop. correct choice", background_color=:transparent, size=(500, 500), dpi=300, foreground_color=:black, labelfontsize=32, tickfontsize=14)
#violin!([1], X[:,1], alpha=0.6, linewidth=0, color=StatsPlots.palette(:Dark2)[4], label="")
bar!([1], [mean(X[:,1])], yerror = [sem(X[:,1])], linewidth=0, msw=5, color=StatsPlots.palette(:Dark2)[4], alpha=0.6, label="")
dotplot!([1], X[:,1], bar_width=0.3,color=:grey, msw=0, label="")
scatter!([1], [mean(X[:,1])], yerror = [sem(X[:,1])],color=StatsPlots.palette(:Dark2)[4], msc=:black, markershape=:circle, markerstrokewidth=5, markersize=0, label="")



## 2/ Partial learning effect (condition 1 VS 2)
tmp = df_env1[(0 .< df_env1.condition .<= 3) .* (df_env1.presInBlock .<= 10) .* .!df_env1.isStable ,:] # after switch
#tmp2 = df_env1[(0 .< df_env1.condition .<= 2) .* (df_env1.negpresInBlock .>= -3) ,:] # before switch, trials that are both <10 trials post switch and <3trials pre switch will count twice
#tmp2.presInBlock .= tmp2.negpresInBlock
#tmp = vcat(tmp, tmp2)
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

#label=["Changing associations" "Stable associations"]
#savefig("WMM/Figures/model_free/expe1_partial_stable.pdf")

## Partial effect (stable)
grp_plot(grpstats1[(0 .< grpstats1.condition .<= 2) .| (grpstats1.condition .== 4),:], "correct", [1,3,4]; xlims=(-3, 9), ylims=(0, 1), xticks=-2:2:9, yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), linestyle=[:solid :solid :dot], labelfontsize=20, tickfontsize=14,  background_color = :white, tickfontcolor=:black, guidefontcolor=:black)
#savefig("WMM/Figures/model_free/expe1_partial_stable.pdf")

## Oddbal effect (need to keep condition 4)
grp_plot(grpstats1[(grpstats1.condition .== 1) .|(grpstats1.condition .== 4) ,:], "correct", "false", [1,6]; xlims=(-3, 9), ylims=(0, 1), xticks=-2:2:9, yticks=0:0.2:1,  xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), label="", labelfontsize=20, tickfontsize=14,  background_color = :transparent, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black, dpi=300)

## Explo
grp_plot(grpstats1[(grpstats1.condition .== 1) .|(grpstats1.condition .== 4) ,:], "explo", "false", [1,6]; xlims=(-3, 9), ylims=(0, 0.3), xticks=-2:2:9, yticks=0:0.2:1,  xlabel="Stimulus presentations", ylabel="Prop. exploratory choice", size=(500, 500), label="", labelfontsize=20, tickfontsize=14,  background_color = :transparent, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black, dpi=300)


## RT
#df_env1[df_env1.condition .== 4,:condition] .= 3 # Merge stable and changing for rt plot
grpstats1 = grp_stats(df_env1)
grp_plot(grpstats1[grpstats1.condition .> 0,:], "rt", [1,2, 3]; xlims=(-3, 9), ylims=(-0.25, 0.25), label=["New rule" "Recurent rule" "Partial new"], xlabel="Stimulus presentations", ylabel="RT (z-score)", size=(500, 500))
#savefig("WMM/Figures/model_free/expe1_rt.pdf")

## Learning complete normal VS Oddball
grp_plot(grpstats1[in.(grpstats1.condition, [[1, 4]]),:], "correct", [1,6]; xlims=(-3, 9), ylims=(0, 1), xticks=([-2, 0, 1, 3, 5,7, 9], [-3, -1, 1, 3, 5,7, 9]), yticks=0:0.2:1,  xlabel="Stimulus presentations", ylabel="Correct choice", size=(500, 500), label="")
savefig("WMM/Figures/model_free/expe1_oddball_correct.pdf")
##
grp_plot(grpstats1[in.(grpstats1.condition, [[1, 4]]),:], "explo", [1,6]; xlims=(-3, 9), ylims=(0,0.3), xticks=([-2, 0, 1, 3, 5,7, 9], [-3, -1, 1, 3, 5,7, 9]), yticks=0:0.1:1, label="", xlabel="Stimulus presentations", ylabel="Non perseverative incorrect choice", size=(500, 500))
savefig("WMM/Figures/model_free/expe1_oddball_explo.pdf")

## Per epoch 
df_env1[!,:epoch] = (df_env1.blockNum .< 14) .+ 2 .* (14 .<= df_env1.blockNum .< 27) .+ 3 .* (27 .<= df_env1.blockNum)
df_epoch = copy(df_env1)
df_epoch.condition .= df_epoch.epoch .+ 3 .* (df_epoch.condition .- 1)

grpstats_epoch1 = grp_stats(df_epoch)

grp_plot(grpstats_epoch1[(grpstats_epoch1.condition .== 10) .| (grpstats_epoch1.condition .== 11) .| (grpstats_epoch1.condition .== 12),:], "correct", [4]; linestyle = [:solid :dash :dot], xlims=(-3, 9), ylims=(0, 1), xticks=-2:2:9, yticks=0:0.2:1,  xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), label="", labelfontsize=20, tickfontsize=14,  background_color = :transparent, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black, dpi=300)


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
dotplot([V[:,2] .- V[:,1], V[:,4] .- V[:,3]], label="", xticks=([1,2],["Lose-shift", "Win-stay"]), ylims = (-0.14, 0.22), ylabel="Proba. action repeat difference\npositive vs negative fb", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)
# X = repeat(0:9,2) ./ (1:5)'
# X[X .> 1] .= NaN
# uX = sort(unique(X))
# Mloss = zeros(length(uX), length(gdf))
# for i in axes(Mloss, 1), j in axes(Mloss, 2)
#     idx = findall(X[1:10,:] .== uX[i])
#     Mloss[i,j] = mean(filter(!isnan, getindex.([V[1:10,:,:]], idx, j)))
# end
# Mwin = zeros(length(uX), length(gdf))
# for i in axes(Mwin, 1), j in axes(Mwin, 2)
#     idx = findall(X[11:20,:] .== uX[i])
#     Mwin[i,j] = mean(filter(!isnan, getindex.([V[11:20,:,:]], idx, j)))
# end

## Environment 2
df_env2 = DataFrame(envdf[[en.task[1] == "WMM2" for en in envdf]])
#sub_split = [15121,17172,18969,19807,22685,32111,32186,44194,57493,78689,83476,91572,96048]
#df_env2 = df_env2[in.(df_env2.subject, [sub_split]),:]
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
plot!(res.clusters[1], 0.95 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[6], alpha=0.5)

#savefig("WMM/Figures/model_free/expe1_relearn_correct.pdf")
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
plot!(res.clusters[1], 0.95 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[2], alpha=0.5)


#savefig("WMM/Figures/model_free/expe1_relearn_explo.pdf")

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
#tmp2 = df_env2[(0 .< df_env2.condition .== 4) .* (df_env2.negpresInBlock .> -3) ,:] # before switch, trials that are both <10 trials post switch and <3trials pre switch will count twice
#tmp2.presInBlock .= tmp2.negpresInBlock .+ 1
#tmp = vcat(tmp, tmp2)
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
#tmp2 = df_env2[(0 .< df_env2.condition .== 4) .* (df_env2.negpresInBlock .> -3) ,:] # before switch, trials that are both <10 trials post switch and <3trials pre switch will count twice
#tmp2.presInBlock .= tmp2.negpresInBlock .+ 1
#tmp = vcat(tmp, tmp2)
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
#plot!(0:9, fill(mean(M.correct_mean), 10), color = :grey, linewidth=3, linestyle=:dash, label="", alpha=0.5)
plot!(res.clusters[1], 0.9 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[2], alpha=0.5)

## t-test on first presentation only :
OneSampleTTest(X[:,1]) # p = 0.16, t-stat = 1.43, df = 50
hline([0.0], color=:grey, linewidth=6, linestyle=:dash, alpha=0.5,label="", xlims=(0.4, 1.6), xticks=[], ylabel="Prop. correct choice", background_color=:transparent, size=(500, 500), dpi=300, foreground_color=:black, labelfontsize=32, tickfontsize=14)
violin!([1], X[:,1], alpha=0.6, linewidth=0, color=StatsPlots.palette(:Dark2)[4], label="")
dotplot!([1], X[:,1], color=StatsPlots.palette(:Dark2)[4], label="")
scatter!([1], [mean(X[:,1])], yerror = [sem(X[:,1])],color=:black, markershape=:circle, markerstrokewidth=5, markersize=12, label="")

#savefig("WMM/Figures/model_free/expe1_partial2_stable.pdf")
## Noise effect
grp_plot(grpstats2[1 .< grpstats2.condition .<= 3,:], "correct", "false", [2,3]; xlims=(-3, 9), xticks=(-2:2:9), ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300)

## Statistical Significance (cluster based permutation test)
tmp = df_env2[(1 .< df_env2.condition .<= 3) .* (df_env2.presInBlock .<= 10) .* .!df_env2.isStable,:] # after switch
# tmp2 = df_env2[(1 .< df_env2.condition .<= 3) .* (df_env2.negpresInBlock .> -3) .* .!df_env2.isStable,:] # before switch, trials that are both <10 trials post switch and <3trials pre switch will count twice
#tmp2.presInBlock .= tmp2.negpresInBlock .+ 1
#tmp = vcat(tmp, tmp2)
gp = groupby(tmp, [:presInBlock, :subject])
cc = combine(gp, [:explo, :condition] => ((x, y) -> mean(x[y .== 3]) - mean( x[y .== 2])) => :diff1_2)
X = fill(NaN, length(unique(cc.subject)), 10)
gp = groupby(cc, :subject)
for i in 1:length(gp)
    for j = 1:10
        jidx = findfirst(gp[i].presInBlock .== j)
        X[i,jidx] = gp[i].diff1_2[jidx]
    end
end
res = cluster_perm_test(X; niter=1e5)
plot!(res.clusters[1], 0.28 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[2], alpha=0.5)


#savefig("WMM/Figures/model_free/expe1_noise_correct.pdf")

## Noise effect (explo)
grp_plot(grpstats2[1 .< grpstats2.condition .<= 3,:], "explo", "false", [2,3]; xlims=(-3, 9), xticks=(-2:2:12),ylims=(0, 0.3), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. non perseverative\nincorrect choice", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)

## Statistical Significance (cluster based permutation test)
tmp = df_env2[(1 .< df_env2.condition .<= 3) .* (df_env2.presInBlock .<= 10) .* .!df_env2.isStable ,:] # after switch
#tmp2 = df_env2[(1 .< df_env2.condition .<= 3) .* (df_env2.negpresInBlock .> -3) ,:] # before switch, trials that are both <10 trials post switch and <3trials pre switch will count twice
#tmp2.presInBlock .= tmp2.negpresInBlock .+ 1
#tmp = vcat(tmp, tmp2)
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


#savefig("WMM/Figures/model_free/expe1_noise_explo.pdf")
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

# grp_plot(grpstats2[(grpstats2.condition .== 4) .| (6 .< grpstats2.condition .<= 9),:], "correct", [6,2,3,4]; linestyle=[:dot :solid :solid :dot], xlims=(-3, 9), xticks=([-2,0,1, 3, 5,7, 9], [-3,-1,1,  3, 5,7, 9]), ylims=(0, 1), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. correct choice", legend=:bottomright, size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)

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


#savefig("WMM/Figures/model_free/expe1_noise_stable.pdf")

## Per epoch
df_env2[!,:epoch] = (df_env2.blockNum .< 14) .+ 2 .* (14 .<= df_env2.blockNum .< 27) .+ 3 .* (27 .<= df_env2.blockNum)
df_epoch = copy(df_env2)
df_epoch.condition .= df_epoch.epoch .+ 3 .* (df_epoch.condition .- 1)

grpstats_epoch2 = grp_stats(df_epoch)

grp_plot(grpstats_epoch2[(27 .< grpstats_epoch2.condition .<= 30),:], "correct", [1]; linestyle = [:solid :dash :dot], xlims=(-3, 9), ylims=(0, 1), xticks=-2:2:9, yticks=0:0.2:1,  xlabel="Stimulus presentations", ylabel="Prop. correct choice", size=(500, 500), label="", labelfontsize=20, tickfontsize=14,  background_color = :transparent, tickfontcolor=:black, guidefontcolor=:black, foreground_color=:black, dpi=300)

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


#savefig("WMM/Figures/model_free/expe1_partial1v2_correct.png")
## Partial effect env1 VS env2 (explo)
grp_plot(gg, "explo", "false", [3,1]; xlims=(-3, 9), xticks=([-2,0,1, 3, 5,7, 9], [-3,-1,1,  3, 5,7, 9]), ylims=(0, 0.3), yticks=0:0.2:1, label="", xlabel="Stimulus presentations", ylabel="Prop. non perseverative\nincorrect choice", size=(500,500), tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)
#savefig("WMM/Figures/model_free/expe1_partial1v2_explo.pdf")
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
#tmp2 = df_env2[(0 .< df_env2.condition .== 4) .* (df_env2.negpresInBlock .> -3) ,:] # before switch, trials that are both <10 trials post switch and <3trials pre switch will count twice
#tmp2.presInBlock .= tmp2.negpresInBlock .+ 1
#tmp = vcat(tmp, tmp2)
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
#tmp2 = df_env1[(0 .< df_env1.nextCondition .== 3) .* (df_env1.negpresInBlock .> -3) ,:] # before switch, trials that are both <10 trials post switch and <3trials pre switch will count twice
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
plot!(0:9, fill(mean(M.correct_mean), 10), color = StatsPlots.palette(:Dark2)[3], linewidth=3, linestyle=:dash, label="", alpha=0.5)
plot!(res.clusters[1], 0.9 .* ones(length(res.clusters[1])), linewidth=5, label="", color = StatsPlots.palette(:Dark2)[3], alpha=0.5)
X_env1 = X

## t-test on first presentation only :
# Env 1 : 
# hline([0.0], color=:grey, linewidth=6, linestyle=:dash, alpha=0.5,label="", xlims=(0.4, 2.6), xticks=[], ylabel="Prop. correct choice", background_color=:transparent, size=(500, 500), dpi=300, foreground_color=:black, labelfontsize=32, tickfontsize=14)
# violin!([1], X_env1[:,1], alpha=0.6, linewidth=0, color=StatsPlots.palette(:Dark2)[3], label="")
# dotplot!([1], X_env1[:,1], color=StatsPlots.palette(:Dark2)[3], label="")
# scatter!([1], [mean(X_env1[:,1])], yerror = [sem(X_env1[:,1])],color=:black, markershape=:circle, markerstrokewidth=5, markersize=12, label="")
bar([1], [mean(X_env1[:,1])],color=StatsPlots.palette(:Dark2)[3], alpha=0.6, linewidth=0, label="", xlims=(0.4, 2.6), xticks=[], ylabel="Prop. correct choice", background_color=:transparent, size=(500, 500), dpi=300, foreground_color=:black, labelfontsize=32, tickfontsize=14)
dotplot!([1], X_env1[:,1], color=:grey, msw=0.0, bar_width=0.3,label="")
scatter!([1], [mean(X_env1[:,1])], yerror = [sem(X_env1[:,1])],color=:black, markershape=:circle, markerstrokewidth=5, markersize=0, label="")

## Env 2 : 
# violin!([2], X_env2[:,1], alpha=0.6, linewidth=0, color=StatsPlots.palette(:Dark2)[5], label="")
# dotplot!([2], X_env2[:,1], color=StatsPlots.palette(:Dark2)[5], label="")
# scatter!([2], [mean(X_env2[:,1])], yerror = [sem(X_env2[:,1])],color=:black, markershape=:circle, markerstrokewidth=5, markersize=12, label="")

bar!([2], [mean(X_env2[:,1])],color=StatsPlots.palette(:Dark2)[5], alpha=0.6, linewidth=0, label="", xlims=(0.4, 2.6), xticks=[], ylabel="Prop. correct choice", background_color=:transparent, size=(500, 500), dpi=300, foreground_color=:black, labelfontsize=32, tickfontsize=14)
dotplot!([2], X_env2[:,1], color=:grey, msw=0.0, bar_width=0.3,label="")
scatter!([2], [mean(X_env2[:,1])], yerror = [sem(X_env1[:,1])],color=:black, markershape=:circle, markerstrokewidth=5, markersize=0, label="")

#OneSampleTTest(X_env1[:,1] .- X_env2[:,1]) # p = 0.01, tstat = -2.47, df= 50



#savefig("WMM/Figures/model_free/expe1_partial1v2_firstpres.png")

## Mutual information task 1 VS task 2
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
PJ_env1 = zeros(4,length(sdf))
for i = 1:length(sdf)
    MI_env1[i], PJ_env1[:,i] = mut_inf(sdf[i].fb)
end

sdf = groupby(df_env2, [:subject, :sessNum])

MI_env2 = zeros(length(sdf))
PJ_env2 = zeros(4,length(sdf))
for i = 1:length(sdf)
    MI_env2[i], PJ_env2[:,i] = mut_inf(sdf[i].fb)
end

## Per state 

sdf = groupby(df_env1, [:subject, :sessNum, :stimulus])

MI_env1_ps = zeros(length(sdf))
PJ_env1_ps = zeros(4,length(sdf))
for i = 1:length(sdf)
    MI_env1_ps[i], PJ_env1_ps[:,i] = mut_inf(sdf[i].fb)
end

sdf = groupby(df_env2, [:subject, :sessNum, :stimulus])

MI_env2_ps = zeros(length(sdf))
PJ_env2_ps = zeros(4,length(sdf))
for i = 1:length(sdf)
    MI_env2_ps[i], PJ_env2_ps[:,i] = mut_inf(sdf[i].fb)
end