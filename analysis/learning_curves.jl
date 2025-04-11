using Turing, Optim, StatsBase, ProgressMeter, CSV, DataFrames
include("preprocessing.jl")

""" Fitting the learning curves of each condition with a saturating exponential function """

## Saturating exponential model
@model sat_expo_WMM1(X, C, Y) = begin
    ρ₀ ~ Beta()
    ρ₁ ~ Beta()    
    pmax ~ Beta()  
    βi ~ MvLogNormal(4, 1.0)
    
    pmin₀ = ρ₀
    pmin₁ = ρ₁
    for i in eachindex(X)
        β = min(50.0, βi[C[i]])#50 * logistic(logβ[C[i]])
        pmin = C[i] < 4 ? pmin₀ : pmin₁
        P = pmin + (pmax - pmin) * (1 - exp(- β * (X[i]-1)))
        if !(zero(P) <= P <= one(P))
            @show P
            @show pmin
            @show pmax
            @show β
        end
        Y[i] ~ Bernoulli(P)
    end
    return
end

@model sat_expo_WMM2(X, C, Y) = begin
    ρ₀ ~ Beta() # pmin changing stims
    ρ₁ ~ Beta() # pmin stable stims 1 (2 changing)
    ρ₂ ~ Beta() # pmin stable stims 2 (1 changing)
    ρ₃ ~ Beta() # pmin stable stims 3 (1 changing + noise)
    pmax ~ Beta()  
    βi ~ MvLogNormal(6, 1.0) # only 3*2 slopes (remove cond 4)
    
    pmin₀ = ρ₀
    pmin₁ = ρ₁
    pmin₂ = ρ₂
    pmin₃ = ρ₃
    for i in eachindex(X)
        β = min(50.0, βi[C[i]])#50 * logistic(logβ[C[i]])
        if C[i] < 4
            pmin = pmin₀
        elseif C[i] == 4
            pmin = pmin₁
        elseif C[i] == 5
            pmin = pmin₂
        elseif C[i] == 6
            pmin = pmin₃
        end
        P = pmin + (pmax - pmin) * (1 - exp(- β * (X[i]-1)))
        if !(zero(P) <= P <= one(P))
            @show P
            @show pmin
            @show pmax
            @show β
        end
        Y[i] ~ Bernoulli(P)
    end
    return
end

@model sat_expo_FU(X, C, Y) = begin
    ρ₀ ~ Beta()
    ρ₁ ~ Beta()
    ρ₂ ~ Beta()    
    pmax ~ Beta()  
    βi ~ MvLogNormal(4, 1.0)
    
    pmin₀ = ρ₀
    pmin₁ = ρ₁
    pmin₂ = ρ₂
    for i in eachindex(X)
        β = min(50.0, βi[C[i]])
        if C[i] < 3
            pmin = pmin₀
        elseif C[i] == 3
            pmin = pmin₁
        elseif C[i] == 4
            pmin = pmin₂
        end
        P = pmin + (pmax - pmin) * (1 - exp(- β * (X[i]-1)))
        if !(zero(P) <= P <= one(P))
            @show P
            @show pmin
            @show pmax
            @show β
        end
        Y[i] ~ Bernoulli(P)
    end
    return
end

## Load data
# Data : 
# folder = "WMM/RAW2/complete"
# flist = filter(x -> occursin(".csv", x), readdir(folder))

# df = DataFrame()

# for f in flist
#     tmp = CSV.read(string(folder, "/", f), DataFrame)
#     tmp = preprocess(tmp)
#     append!(df, tmp)
# end

# Simu : 
model = "SI_MultVol_SampleAction"
task = 1
tasknames = ["WMM1", "WMM2"]

pathname = string("WMM/Simus/", model, "_", tasknames[task])
fileslist = filter(x -> occursin(".csv", x), readdir(pathname))
df = DataFrame()

wb = Progress(length(fileslist), 1, "loading...")
for fi in eachindex(fileslist)
    fl = fileslist[fi]
    tmp = CSV.read(string(pathname, "/", fl), DataFrame)
    tmp.subject .+= fi * 1000
    pres_in_block!(tmp)
    if in("ExploOut", names(tmp))
        append!(df, tmp)
    end
    next!(wb)
end


# Split by Environment
envdf = groupby(df, :task)

## Environment 1 
df_env1 = DataFrame(envdf[[en.task[1] == "WMM1" for en in envdf]])
df_env1 = df_env1[df_env1.condition .> 0,:] # Remove condition 0 (starting blocks)
df_env1.condition[df_env1.condition .== 4] .= 1
for t = 1:nrow(df_env1)
    if df_env1[t, Symbol("isStable_$(df_env1.stimulus[t])")] && df_env1.condition[t] > 0
        df_env1.condition[t] = 4 # Make stable stims a special condition
    end
end


## Split per subject and condition
gdf = groupby(df_env1, :subject)

# Fit per subject and condition
function fit_learning_curve_WMM1(gdf)
    N = length(gdf)

    fit_results = DataFrame(subject = zeros(Int, N), pmin₀ = zeros(N), pmin₁ = zeros(N), pmax = zeros(N), beta_1 = zeros(N), beta_2 = zeros(N), beta_3 = zeros(N), beta_4 = zeros(N), rec_effect = zeros(N), part_effect = zeros(N), stable_effect = zeros(N), converged = falses(N))

    wb = Progress(N, 1, "fitting....")
    for gi in 1:N
        g = gdf[gi]
        x = g.presInBlock
        y = g.correct
        c = g.condition

        mdl = sat_expo_WMM1(x,c,y)
        chn = sample(mdl, NUTS(), MCMCThreads(), 1000, 4)
        # nt = mean(chn).nt
        # res = NamedTuple{(nt.parameters...,)}(nt.mean)
        nt = get(chn, [:ρ₀, :ρ₁, :pmax, :βi])
        
        res = (ρ₀ = mean(logit.(nt.ρ₀)), 
        ρ₁ = mean(logit.(nt.ρ₁)),
        pmax = mean(logit.(nt.pmax)),
        βi = [mean(log.(nt.βi[i])) for i = 1:4])
        
        fit_results[gi, :subject] = g.subject[1]
        fit_results[gi, :pmin₀] = res[:ρ₀]
        fit_results[gi, :pmin₁] = res[:ρ₁]
        fit_results[gi, :pmax] = res[:pmax]
        for j = 1:4
            fit_results[gi, Symbol("beta_$j")] = res[:βi][j]#res[Symbol("βi[$j]")]
        end
        fit_results[gi, :converged] = all(ess_rhat(chn).nt.rhat .< 1.05)
        fit_results[gi, :rec_effect] = mean(chn[Symbol("βi[1]")] .< chn[Symbol("βi[2]")])
        fit_results[gi, :part_effect] = mean(chn[Symbol("βi[3]")] .< chn[Symbol("βi[1]")])
        fit_results[gi, :stable_effect] = mean(chn[:ρ₁] .< chn[:pmax])
        next!(wb)
    end
    return fit_results
end
fit_results = fit_learning_curve_WMM1(gdf)

## Environment 2
df_env2 = DataFrame(envdf[[en.task[1] == "WMM2" for en in envdf]])
df_env2 = df_env2[0 .< df_env2.condition .< 4,:] # Remove conditions 0 (starting blocks) and 4 (only noise)
for t = 1:nrow(df_env2)
    if df_env2[t, Symbol("isStable_$(df_env2.stimulus[t])")]
        df_env2.condition[t] += 3 # Make stable stims a special condition
    end
end

## Split per subject and condition
gdf = groupby(df_env2, :subject)

# Fit per subject and condition
function fit_learning_curve_WMM2(gdf)
    N = length(gdf)

    fit_results = DataFrame(subject = zeros(Int, N), pmin₀ = zeros(N), pmin₁ = zeros(N), pmin₂ = zeros(N), pmin₃ = zeros(N), pmax = zeros(N), beta_1 = zeros(N), beta_2 = zeros(N), beta_3 = zeros(N), beta_4 = zeros(N), beta_5 = zeros(N), beta_6 = zeros(N), part_effect_1 = zeros(N), part_effect_2 = zeros(N), stable_effect_1 = zeros(N), stable_effect_2 = zeros(N), stable_effect_3 = zeros(N), converged = falses(N))

    wb = Progress(N, 1, "fitting....")
    for gi in 1:N
        g = gdf[gi]
        x = g.presInBlock
        y = g.correct
        c = g.condition

        mdl = sat_expo_WMM2(x,c,y)
        chn = sample(mdl, NUTS(), MCMCThreads(), 1000, 4)
        nt = mean(chn).nt
        res = NamedTuple{(nt.parameters...,)}(nt.mean)
        
        fit_results[gi, :subject] = g.subject[1]
        fit_results[gi, :pmin₀] = res[:ρ₀]
        fit_results[gi, :pmin₁] = res[:ρ₁]
        fit_results[gi, :pmin₂] = res[:ρ₂]
        fit_results[gi, :pmin₃] = res[:ρ₃]
        fit_results[gi, :pmax] = res[:pmax]
        for j = 1:6
            fit_results[gi, Symbol("beta_$j")] = res[Symbol("βi[$j]")]
        end
        fit_results[gi, :converged] = all(ess_rhat(chn).nt.rhat .< 1.05)
        fit_results[gi, :part_effect_1] = mean(chn[Symbol("βi[2]")] .< chn[Symbol("βi[1]")])
        fit_results[gi, :part_effect_2] = mean(chn[Symbol("βi[3]")] .< chn[Symbol("βi[1]")])
        fit_results[gi, :stable_effect_1] = mean(chn[:ρ₁] .< chn[:pmax])
        fit_results[gi, :stable_effect_2] = mean(chn[:ρ₂] .< chn[:pmax])
        fit_results[gi, :stable_effect_3] = mean(chn[:ρ₃] .< chn[:pmax])
        next!(wb)
    end
    return fit_results
end
fit_results = fit_learning_curve_WMM2(gdf)

## Follow-up

# Data : 
folder = "WMM/RAW_FU"
flist = filter(x -> occursin(".csv", x), readdir(folder))

df = DataFrame()

for f in flist
    tmp = CSV.read(string(folder, "/", f), DataFrame)
    tmp = preprocess(tmp)
    append!(df, tmp)
end

# Simu : 
# df = CSV.read("/Users/sami/PhD/Model_Tasks_Data/Data/WMM/Simus/all_simus/cross_simus/SI_MultVol_SampleAction_FU3_CV.csv", DataFrame)

# Split by Environment
envdf = groupby(df, :task)

## Environment 1 
df_env1 = DataFrame(envdf[[en.task[1] == "task3" for en in envdf]])
# df_env1.condition[df_env1.condition .== 4] .= 1
df_env1 = df_env1[df_env1.condition .> 0,:] # Remove condition 0 (starting blocks)
# for t = 1:nrow(df_env1)
#     if df_env1[t, Symbol("isStable_$(df_env1.stimulus[t])")] && df_env1.condition[t] > 0
#         df_env1.condition[t] += 2 # Make stable stims a special condition
#     end
# end

## FOR FU3
df_chg = deepcopy(df_env1)
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
pres_in_block!(df_chg)
df_env1 = deepcopy(df_chg)

## Split per subject and condition
gdf = groupby(df_env1, :subject)

# Fit per subject and condition
function fit_learning_curve_FU(gdf)
    N = length(gdf)

    fit_results = DataFrame(subject = zeros(Int, N), pmin₀ = zeros(N), pmin₁ = zeros(N), pmin₂ = zeros(N), pmax = zeros(N), beta_1 = zeros(N), beta_2 = zeros(N), beta_3 = zeros(N), beta_4 = zeros(N), rec_effect = zeros(N), part_effect_1 = zeros(N), part_effect_2 = zeros(N), stable_effect_1 = zeros(N), stable_effect_2 = zeros(N), converged = falses(N))

    wb = Progress(N, 1, "fitting....")
    for gi in 1:N
        g = gdf[gi]
        x = g.presInBlock
        y = g.correct
        c = g.condition

        mdl = sat_expo_FU(x,c,y)
        chn = sample(mdl, NUTS(), MCMCThreads(), 1000, 4)
        nt = mean(chn).nt
        res = NamedTuple{(nt.parameters...,)}(nt.mean)
        
        fit_results[gi, :subject] = g.subject[1]
        fit_results[gi, :pmin₀] = res[:ρ₀]
        fit_results[gi, :pmin₁] = res[:ρ₁]
        fit_results[gi, :pmin₂] = res[:ρ₂]
        fit_results[gi, :pmax] = res[:pmax]
        for j = 1:4
            fit_results[gi, Symbol("beta_$j")] = res[Symbol("βi[$j]")]
        end
        fit_results[gi, :converged] = all(ess_rhat(chn).nt.rhat .< 1.05)
        fit_results[gi, :rec_effect] = mean(chn[Symbol("βi[1]")] .< chn[Symbol("βi[2]")])
        fit_results[gi, :part_effect_1] = mean(chn[Symbol("βi[3]")] .< chn[Symbol("βi[1]")])
        fit_results[gi, :part_effect_2] = mean(chn[Symbol("βi[4]")] .< chn[Symbol("βi[1]")])
        fit_results[gi, :stable_effect_1] = mean(chn[:ρ₁] .< chn[:pmax])
        fit_results[gi, :stable_effect_2] = mean(chn[:ρ₂] .< chn[:pmax])
        next!(wb)
    end
    return fit_results
end
fit_results = fit_learning_curve_FU(gdf)

## Figures
# df = CSV.read("WMM/learning_curves_data_WMM1.csv", DataFrame)
#df2 = CSV.read("WMM/learning_curves_data_WMM2.csv", DataFrame)
#df3 = CSV.read("WMM/learning_curves_data_FU2.csv", DataFrame)
#df4 = CSV.read("WMM/learning_curves_data_FU3.csv", DataFrame)
# df2 = CSV.read("WMM/Simus/learning_curves_adaptive_alpha.csv", DataFrame)
# df3 = CSV.read("WMM/Simus/learning_curves_adaptive_beta.csv", DataFrame)
df_data1 = CSV.read("WMM/learning_curves_data_FU1.csv", DataFrame)
df_data2 = CSV.read("WMM/learning_curves_data_FU2.csv", DataFrame)
df_data3 = CSV.read("WMM/learning_curves_data_FU3.csv", DataFrame)
#df_model = CSV.read("WMM/Simus/learning_curves_full_adaptive_cfql2.csv", DataFrame)
# df_model = CSV.read("WMM/Simus/learning_curves_probe_LRRel.csv", DataFrame)
#df_model = CSV.read("WMM/Simus/learning_curves_SI_MultVol_SampleAction.csv", DataFrame)

## Basal learning speed and recurrence
simu_m = [mean(exp.(df_model.beta_1)), mean(exp.(df_model.beta_2 )), mean(exp.(df_model.beta_3 ))]
simu_s = [std(exp.(df_model.beta_1))/sqrt(51), std(exp.(df_model.beta_2 ))/sqrt(51), std(exp.(df_model.beta_3 ))/sqrt(51)]

scatter([1], [simu_m[1]], yerror=[simu_s[1]], markersize=20, color=:white, msw=10, msc = StatsPlots.palette(:Dark2)[1], label="")
scatter!([1.7], [simu_m[2]], yerror=[simu_s[2]], markersize=20, color=:white, msw=10, msc = StatsPlots.palette(:Dark2)[2], label="")

data_m = [mean(exp.(df_data.beta_1)), mean(exp.(df_data.beta_2 )), mean(exp.(df_data.beta_3 ))]
data_s = [std(exp.(df_data.beta_1))/sqrt(51), std(exp.(df_data.beta_2 ))/sqrt(51), std(exp.(df_data.beta_3 ))/sqrt(51)]

scatter!([0.8], [data_m[1]], yerror=[data_s[1]], markersize=20, color=StatsPlots.palette(:Dark2)[1], msw=10, msc = StatsPlots.palette(:Dark2)[1], label="")
scatter!([1.5], [data_m[2]], yerror=[data_s[2]], markersize=20, color=StatsPlots.palette(:Dark2)[2], msw=10, msc = StatsPlots.palette(:Dark2)[2], label="", xaxis=:off, xticks=[], xlims=(0.6, 1.9), ylims=(0.39, 0.74),ylabel="Learning speed", labelfontsize=32, tickfontsize = 14, background_color=:transparent, foreground_color=:black, size=(500, 500), dpi=300)


## Interference (changing)
bar([mean( (exp.(df_model.beta_1))) mean( (exp.(df_model.beta_3)) )  ], alpha=1)

#scatter!([1], [mean( (exp.(df_data.beta_3) .- exp.(df_data.beta_1 ))  )])

##
V = log.(df.beta_2 ./ df.beta_1)
hline([0], color=:black, linewidth=3, linestyle= :dash, label="", size=(500, 500))
violin!([1], V, label="", xlim=(0.2, 1.8), color=RGBA(217/255,95/255,2/255,0.6), linewidth=0, yticks=[-1, 0, 1, 2], ytickfontsize=20, size=(500, 500), xticks=[],axis=false)
dotplot!([1], V, label="", color=RGBA(217/255,95/255,2/255,1), xticks=[], yticks=[-1, 0, 1, 2], ytickfontsize=20, size=(500, 500), axis=false,tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)

##
V = log.(df.beta_3 ./ df.beta_1)
hline([0], color=:black, linewidth=3, linestyle= :dash, label="", size=(500, 500))
violin!([1], V, label="", xlim=(0.2, 1.8), palette=StatsPlots.palette(:Dark2)[[3]], linewidth=0, yticks=[-1, 0, 1, 2], ytickfontsize=20, size=(500, 500), xticks=[],axis=false)
dotplot!([1], V, label="", palette=StatsPlots.palette(:Dark2)[[3]], xticks=[], yticks=[-1, 0, 1, 2], ytickfontsize=20, size=(500, 500), axis=false,tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)
##
V1 = log.(df.beta_2 ./ df.beta_1)
V2 = log.(df2.beta_2 ./ df2.beta_1)
hline([0], color=:black, linewidth=3, linestyle= :dash, label="", size=(500, 500))
violin!([1 2], [V1, V2], label="", xlim=(0.2, 2.8), color=RGBA(217/255,95/255,2/255,0.6), linewidth=0, yticks=[-1, 0, 1, 2], ytickfontsize=20, size=(500, 500), xticks=[],axis=false)
dotplot!([1 2], [V1, V2], label="", color=RGBA(217/255,95/255,2/255,1), xticks=[], yticks=[-1, 0, 1, 2], ytickfontsize=20, size=(500, 500), axis=false,tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)

##
X = log.(df.pmax ./ df.pmin₁)
hline([0], color=:black, linewidth=3, linestyle= :dash, label="", size=(500, 500))
violin!([1], X, label="", xlim=(0.2, 1.8), color=RGBA(231/255,41/255,138/255,0.6), linewidth=0)
dotplot!([1], X, label="", color=RGBA(231/255,41/255,138/255,1), xticks=[], yticks=-0.2:0.2:0.8,ytickfontsize=20, size=(500, 500), axis=false,tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)

##
X1 = log.(df.pmax ./ df.pmin₁)
X2 = log.(df3.pmax ./ df3.pmin₁)
X3 = log.(df4.pmax ./ df4.pmin₁)
X4 = log.(df2.pmax ./ df2.pmin₁)
hline([0], color=:black, linewidth=3, linestyle= :dash, label="", size=(500, 500))
violin!([1 2 3 4], [X1, X2, X3, X4], label="", xlim=(0.2, 4.8), palette=StatsPlots.palette(:Dark2)[[4, 4, 3, 3]], alpha=0.6, linewidth=0)
dotplot!([1 2 3 4], [X1, X2, X3, X4], label="", palette=StatsPlots.palette(:Dark2)[[4, 4, 3, 3]], xticks=[], yticks=-0.6:0.3:0.9,ytickfontsize=20, axis=false,tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, size=(1000, 500))

##
MP = [mean(df.pmax), mean(df2.pmax), mean(df3.pmax), mean(df4.pmax), mean(df5.pmax), mean(df6.pmax)]
SP = [std(df.pmax), std(df2.pmax), std(df3.pmax), std(df4.pmax), std(df5.pmax), std(df6.pmax)]
MB = [mean(df.beta_1), mean(df2.beta_1), mean(df3.beta_1), mean(df4.beta_1), mean(df5.beta_1), mean(df6.beta_1)]
SB = [std(df.beta_1), std(df2.beta_1), std(df3.beta_1), std(df4.beta_1), std(df5.beta_1), std(df6.beta_1)]

MR = [mean(log.(df.beta_2 ./ df.beta_1)), mean(log.(df2.beta_2 ./ df2.beta_1)), mean(log.(df3.beta_2 ./ df3.beta_1)), mean(log.(df4.beta_2 ./ df4.beta_1)), mean(log.(df5.beta_2 ./ df5.beta_1)), mean(log.(df6.beta_2 ./ df6.beta_1))]
SR = [std(log.(df.beta_2 ./ df.beta_1)), std(log.(df2.beta_2 ./ df2.beta_1)), std(log.(df3.beta_2 ./ df3.beta_1)), std(log.(df4.beta_2 ./ df4.beta_1)), std(log.(df5.beta_2 ./ df5.beta_1)), std(log.(df6.beta_2 ./ df6.beta_1))]

MS = [mean(log.(df.pmax ./ df.pmin₁)), mean(log.(df2.pmax ./ df2.pmin₁)), mean(log.(df3.pmax ./ df3.pmin₁)), mean(log.(df4.pmax ./ df4.pmin₁)), mean(log.(df5.pmax ./ df5.pmin₁)), mean(log.(df6.pmax ./ df6.pmin₁))]
SS = [std(log.(df.pmax ./ df.pmin₁)), std(log.(df2.pmax ./ df2.pmin₁)), std(log.(df3.pmax ./ df3.pmin₁)), std(log.(df4.pmax ./ df4.pmin₁)), std(log.(df5.pmax ./ df5.pmin₁)), std(log.(df6.pmax ./ df6.pmin₁))]

##
scatter(MB', MP',xerror = SB' ./ sqrt(51), yerror=SP' ./ sqrt(51), label=["Subjects" "Adaptive α" "Adaptive β" "Adaptive α + β" "Probe" "Strategic inference"], legend_position=:topleft, palette = vcat(colorant"black", ColorBrewer.palette("Set1", 5)),  msc=vcat(colorant"black", ColorBrewer.palette("Set1", 5))', msw=3, markershape = [:diamond :circle :star5 :hexagon :dtriangle :rect], markersize=8, xlabel="Learning speed", ylabel="Asymptotic performance", ylims=(0.72, 0.81))
#savefig("WMM/Figures/model_based/expe1_learning_curves1.pdf")
##
scatter(MB', MP',xerror = SB' ./ sqrt(51), yerror=SP' ./ sqrt(51), label=["Subjects" "Adaptive α + Reset" "Strategic inference"], legend_position=:topleft, palette = vcat(colorant"black", ColorBrewer.palette("Set1", 5)[[1,5]]),  msc=vcat(colorant"black", ColorBrewer.palette("Set1", 5)[[1,5]])', msw=3, markershape = [:diamond :circle :star5 :hexagon :dtriangle :rect], markersize=8, xlabel="Learning speed", ylabel="Asymptotic performance", ylims=(0.72, 0.81))
savefig("WMM/Figures/model_based/expe1_learning_curves_adaptive_alpha_resetVS_strategic_inference.pdf")
##

scatter(MR', MS',xerror = SR' ./sqrt(51) , yerror=SS' ./ sqrt(51), label=["Subjects" "Adaptive α" "Adaptive β" "Adaptive α + β" "Probe" "Strategic inference"], legend_position=:topleft, palette = vcat(colorant"black", ColorBrewer.palette("Set1", 5)),  msc=vcat(colorant"black", ColorBrewer.palette("Set1", 5))', msw=3, markershape = [:diamond :circle :star5 :hexagon :dtriangle :rect], markersize=8, xlabel="Recurence effect", ylabel="Global inteference effect")
savefig("WMM/Figures/model_based/expe1_learning_curves2.pdf")

## Non parametric statistics for simulations
function test_simu(df, varname1, varname2; niter=10_000)
    N = 51
    rep = 10
    simgrp = zeros(N)
    res = falses(niter)
    for ni = 1:niter
        for i = 1:N
            si = rand(1:N*rep)#(i-1) + rand(1:rep)#
            simgrp[i] = (df[si,varname1]) - (df[si,varname2])
        end
        tt = SignedRankTest(simgrp)#OneSampleTTest(simgrp)#
        res[ni] = pvalue(tt; tail=:both) < 0.05
    end
    return mean(res)
end

function simuVSdata(df_data, df_simu, varname; niter=10_000, nsub = 54, rep=20)
    simgrp = zeros(nsub)
    res = falses(niter)
    for ni = 1:niter
        for i = 1:nsub
            si = (i-1)*rep + rand(1:rep)
            simgrp[i] = df_data[i, varname] - df_simu[si, varname]
        end
        tt = SignedRankTest(simgrp)
        res[ni] = pvalue(tt) < 0.05
    end
    return mean(res)
end