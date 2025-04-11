parmdf = CSV.read("/Users/sami/PhD/Model_Tasks_Data/Data/WMM/Fits/Raw2/SI_MultVol_SampleAction_WMM1/Summary/summarystats.csv", DataFrame)
parmdf2 = CSV.read("/Users/sami/PhD/Model_Tasks_Data/Data/WMM/Fits/Raw2/SI_MultVol_SampleAction_WMM2/Summary/summarystats.csv", DataFrame)

parmdf[!,:logβ] = log.(parmdf[:,:β])
parmdf2[!,:logβ] = log.(parmdf2[:,:β])
parmdf[!,:logρ] = log.(parmdf[:,:ρ])
parmdf2[!,:logρ] = log.(parmdf2[:,:ρ])

##
#ω = log.(softmax(Matrix(parmdf[:,Symbol.(["ω[$i]" for i = 1:5])]), dims=2))
ω = Matrix(parmdf[:,Symbol.(["ω[$i]" for i = 1:5])])
#ω .-= ω[:,5]
##
# ω2 = log.(softmax(Matrix(parmdf2[:,Symbol.(["ω[$i]" for i = 1:5])]), dims=2))
ω2 = Matrix(parmdf2[:,Symbol.(["ω[$i]" for i = 1:5])])
#ω2 .-= ω2[:,5]


##
bar([1], [mean(parmdf.logρ)], color= StatsPlots.palette(:Dark2)[3], linewidth=0, alpha=0.6, bar_width= 0.7, label="", xlims=(0.2, 8.8), xticks=[], xaxis=false, ylabel="Evidence weight ρ", ylims=(-0.6, 0.0), yticks=-0.5:0.2:0, size=(500, 500), bottommargin=10mm, tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300)
@df parmdf dotplot!([1], :logρ, color = :gray, bar_width=0.2, msw=0.0, label="")
scatter!([1], [mean(parmdf.logρ)], yerror=[std(parmdf.logρ)/sqrt(nrow(parmdf))], color=:black, markersize=0.0, markerstrokewidth=5, label="")
bar!([2], [mean(parmdf2.logρ)], color= StatsPlots.palette(:Dark2)[5], linewidth=0, alpha=0.6, bar_width=0.7, label="")
@df parmdf2 dotplot!([2], :logρ, color = :gray, bar_width=0.2, msw=0.0, label="")
scatter!([2], [mean(parmdf2.logρ)], yerror=[std(parmdf2.logρ)/sqrt(nrow(parmdf2))], color=:black, markersize=0.0, markerstrokewidth=5, label="")
#@df parmdf violin([1], :logρ, side = :left, color = StatsPlots.palette(:Dark2)[3], xticks=[], xlims=(0.2, 1.8), xaxis=false, ylabel="log(ρ)", label="", alpha=0.6, linewidth=0)
#@df parmdf dotplot!([1], :logρ, side = :left, color = StatsPlots.palette(:Dark2)[3], label="", size=(500, 500))
#@df parmdf2 violin!([1], :logρ, side = :right, color = StatsPlots.palette(:Dark2)[5], xticks=[], xlims=(0.2, 1.8), label="", alpha=0.6, linewidth=0)
#@df parmdf2 dotplot!([1], :logρ, side = :right, color = StatsPlots.palette(:Dark2)[5], label="",  tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300)

#scatter!([0.85], [mean(parmdf.logρ)], yerror=[std(parmdf.logρ)/sqrt(nrow(parmdf))], color=:black, markersize=8, markerstrokewidth=3, label="")
#scatter!([1.15], [mean(parmdf2.logρ)], yerror=[std(parmdf2.logρ)/sqrt(nrow(parmdf2))], color=:black, markersize=8, msw=3,label="")
#savefig("WMM/Figures/model_based/expe1_SI_MV_SA_ρ.png")

##
bar([0, 3, 6], vec(mean(ω[:,1:3], dims=1)), color=StatsPlots.palette(:Dark2)[3], label="", xlims=(-0.8, 7.8), xticks=([0.5, 3.5, 6.5], ["2", "1", "0"]), xlabel="Strategies similarity", ylims=(-3.5, 1.5), yticks=-6:1:2, ylabel="Transition bias ω", alpha=0.6, linewidth=0, bar_width=0.7,size=(500, 500), tickfontsize=18, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300)
bar!([1, 4, 7], vec(mean(ω2[:,1:3], dims=1)), color=StatsPlots.palette(:Dark2)[5], label="", alpha=0.6, linewidth=0,  bar_width=0.7)

dotplot!([0 3 6], ω[:,1:3], color=:grey, markersize=3, msw=0.0, bar_width=0.2, label="")
dotplot!([1 4 7], ω2[:,1:3], color=:grey, markersize=3, msw=0.0, bar_width=0.2, label="")

scatter!([0, 3, 6], vec(mean(ω[:,1:3], dims=1)), yerror=vec(std(ω[:,1:3], dims=1))./sqrt(size(ω,1)), color=:black, markersize=0, markerstrokewidth=5, label="")
scatter!([1, 4, 7], vec(mean(ω2[:,1:3], dims=1)), yerror=vec(std(ω2[:,1:3], dims=1))./sqrt(size(ω2,1)), color=:black, markersize=0, msw=5,label="")


# violin((1:3)', ω[:,1:3], side=:left, color=StatsPlots.palette(:Dark2)[3], label="", xticks=1:3, xlabel="Distance between strategies", ylabel="log(ω)", alpha=0.6, linewidth=0)
# dotplot!((1:3)', ω[:,1:3], side=:left, color=StatsPlots.palette(:Dark2)[3], label="", tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)

# violin!((1:3)', ω2[:,1:3], side=:right, color=StatsPlots.palette(:Dark2)[5], label="", alpha=0.6, linewidth=0)
# dotplot!((1:3)', ω2[:,1:3], side=:right, color=StatsPlots.palette(:Dark2)[5], label="", tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300, size=(700, 600))


# scatter!([0.85, 1.85, 2.85], vec(mean(ω[:,1:3], dims=1)), yerror=vec(std(ω[:,1:3], dims=1))./sqrt(size(ω,1)), color=:black, markersize=8, markerstrokewidth=3, label="")
# scatter!([1.15, 2.15, 3.15], vec(mean(ω2[:,1:3], dims=1)), yerror=vec(std(ω2[:,1:3], dims=1))./sqrt(size(ω2,1)), color=:black, markersize=8, msw=3,label="")
#savefig("WMM/Figures/model_based/expe1_SI_MV_SA_ω.png")

##
@df parmdf violin([1], :ωᵣ, side = :left, color = StatsPlots.palette(:Dark2)[3], xticks=[], xlims=(0.2, 1.8), xaxis=false, ylims= (0.0, 1), label="", alpha=0.6, linewidth=0)
@df parmdf dotplot!([1], :ωᵣ, side = :left, color = StatsPlots.palette(:Dark2)[3], label="", size=(500, 500))
@df parmdf2 violin!([1], :ωᵣ, side = :right, color = StatsPlots.palette(:Dark2)[1], xticks=[], xlims=(0.2, 1.8), label="", alpha=0.6, linewidth=0)
@df parmdf2 dotplot!([1], :ωᵣ, side = :right, color = StatsPlots.palette(:Dark2)[1], label="",  tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)

##
@df parmdf violin([1], :ω₀, side = :left, color = StatsPlots.palette(:Dark2)[3], xticks=[], xlims=(0.2, 1.8), xaxis=false, ylims= (0.0, 1), label="", alpha=0.6, linewidth=0)
@df parmdf dotplot!([1], :ω₀, side = :left, color = StatsPlots.palette(:Dark2)[3], label="", size=(500, 500))
@df parmdf2 violin!([1], :ω₀, side = :right, color = StatsPlots.palette(:Dark2)[1], xticks=[], xlims=(0.2, 1.8), label="", alpha=0.6, linewidth=0)
@df parmdf2 dotplot!([1], :ω₀, side = :right, color = StatsPlots.palette(:Dark2)[1], label="",  tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)

##
@df parmdf violin([1], :ϵ, side = :left, color = StatsPlots.palette(:Dark2)[3], xticks=[], xlims=(0.2, 1.8), xaxis=false, ylims= (0.0, 0.8), label="", alpha=0.6, linewidth=0)
@df parmdf dotplot!([1], :ϵ, side = :left, color = StatsPlots.palette(:Dark2)[3], label="", size=(500, 500))
@df parmdf2 violin!([1], :ϵ, side = :right, color = StatsPlots.palette(:Dark2)[1], xticks=[], xlims=(0.2, 1.8), label="", alpha=0.6, linewidth=0)
@df parmdf2 dotplot!([1], :ϵ, side = :right, color = StatsPlots.palette(:Dark2)[1], label="",  tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)

##
@df parmdf violin([1], :logβ, side = :left, color = StatsPlots.palette(:Dark2)[3], xticks=[], xlims=(0.2, 1.8), xaxis=false, label="", alpha=0.6, linewidth=0)
@df parmdf dotplot!([1], :logβ, side = :left, color = StatsPlots.palette(:Dark2)[3], label="", size=(500, 500))
@df parmdf2 violin!([1], :logβ, side = :right, color = StatsPlots.palette(:Dark2)[1], xticks=[], xlims=(0.2, 1.8), label="", alpha=0.6, linewidth=0)
@df parmdf2 dotplot!([1], :logβ, side = :right, color = StatsPlots.palette(:Dark2)[1], label="",  tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)

##
bar([1], vec(mean(ω[:,5], dims=1)), color=StatsPlots.palette(:Dark2)[3], label="", xlims=(0.2, 8.8), xticks=[], xlabel="", ylims=(-0.5, 4.5), yticks=-1:1:5, ylabel="Perceived stability", alpha=0.6, linewidth=0, bar_width=0.7,size=(500, 500), tickfontsize=18, labelfontsize=20, background_color=:transparent, foreground_color=:black, dpi=300)
bar!([2], vec(mean(ω2[:,5], dims=1)), color=StatsPlots.palette(:Dark2)[5], label="", alpha=0.6, linewidth=0,  bar_width=0.7)

dotplot!([1], ω[:,5], color=:grey, markersize=3, msw=0.0, bar_width=0.2, label="")
dotplot!([2], ω2[:,5], color=:grey, markersize=3, msw=0.0, bar_width=0.2, label="")

scatter!([1], vec(mean(ω[:,5], dims=1)), yerror=vec(std(ω[:,5], dims=1))./sqrt(size(ω,1)), color=:black, markersize=0, markerstrokewidth=5, label="")
scatter!([2], vec(mean(ω2[:,5], dims=1)), yerror=vec(std(ω2[:,5], dims=1))./sqrt(size(ω2,1)), color=:black, markersize=0, msw=5,label="")


##
@df parmdf violin([1], :WAIC, side = :left, color = StatsPlots.palette(:Dark2)[3], xticks=[], xlims=(0.2, 1.8), xaxis=false, label="", alpha=0.6, linewidth=0)
@df parmdf dotplot!([1], :WAIC, side = :left, color = StatsPlots.palette(:Dark2)[3], label="", size=(500, 500))
@df parmdf2 violin!([1], :WAIC, side = :right, color = StatsPlots.palette(:Dark2)[1], xticks=[], xlims=(0.2, 1.8), label="", alpha=0.6, linewidth=0)
@df parmdf2 dotplot!([1], :WAIC, side = :right, color = StatsPlots.palette(:Dark2)[1], label="",  tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)

##
#@df parmdf violin([1], :logβ, side = :left, color = StatsPlots.palette(:Dark2)[3], xticks=[], xlims=(0.2, 1.8), xaxis=false, ylims= (0.0, 0.8), label="", alpha=0.6, linewidth=0)
@df parmdf dotplot([1], :logβ, side = :left, color = StatsPlots.palette(:Dark2)[3], label="", size=(500, 500))
#@df parmdf2 violin!([1], :logβ, side = :right, color = StatsPlots.palette(:Dark2)[1], xticks=[], xlims=(0.2, 1.8), label="", alpha=0.6, linewidth=0)
@df parmdf2 dotplot!([1], :logβ, side = :right, color = StatsPlots.palette(:Dark2)[1], label="",  tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black)


## Follow up task
parmdf = CSV.read("/Users/sami/PhD/Model_Tasks_Data/Data/WMM/Fits/Raw2/SI_MultVol_SampleAction_task1/Summary/summarystats.csv", DataFrame)
parmdf2 = CSV.read("/Users/sami/PhD/Model_Tasks_Data/Data/WMM/Fits/Raw2/SI_MultVol_SampleAction_task2/Summary/summarystats.csv", DataFrame)
parmdf3 = CSV.read("/Users/sami/PhD/Model_Tasks_Data/Data/WMM/Fits/Raw2/SI_MultVol_SampleAction_task3/Summary/summarystats.csv", DataFrame)
keep_subjects = intersect(unique(parmdf.subject),unique(parmdf2.subject),unique(parmdf3.subject))

parmdf = parmdf[in.(parmdf.subject, [keep_subjects]),:]
parmdf2 = parmdf2[in.(parmdf2.subject, [keep_subjects]),:]
parmdf3 = parmdf3[in.(parmdf3.subject, [keep_subjects]),:]

parmdf[!,:logρ] = log.(parmdf[:,:ρ])
parmdf2[!,:logρ] = log.(parmdf2[:,:ρ])
parmdf3[!,:logρ] = log.(parmdf3[:,:ρ])

# ω = log.(softmax(Matrix(parmdf[:,Symbol.(["ω[$i]" for i = 1:5])]), dims=2))
# ω2 = log.(softmax(Matrix(parmdf2[:,Symbol.(["ω[$i]" for i = 1:5])]), dims=2))
# ω3 = log.(softmax(Matrix(parmdf3[:,Symbol.(["ω[$i]" for i = 1:5])]), dims=2))


ω = Matrix(parmdf[:,Symbol.(["ω[$i]" for i = 1:5])])
ω2 = Matrix(parmdf2[:,Symbol.(["ω[$i]" for i = 1:5])])
ω3 = Matrix(parmdf3[:,Symbol.(["ω[$i]" for i = 1:5])])

## Rho
@df parmdf violin([1], :logρ, color = StatsPlots.palette(:Dark2)[1],  ylabel="log(ρ)", label="", alpha=0.6, linewidth=0)
@df parmdf dotplot!([1], :logρ, color = StatsPlots.palette(:Dark2)[1], label="")
@df parmdf2 violin!([2], :logρ, color = StatsPlots.palette(:Dark2)[4],  label="", alpha=0.6, linewidth=0)
@df parmdf2 dotplot!([2], :logρ, color = StatsPlots.palette(:Dark2)[4], label="",  tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, size=(500, 500), dpi=300)
@df parmdf3 violin!([3], :logρ, color = StatsPlots.palette(:Dark2)[6], xticks=1:3, xlabel="Context", label="", alpha=0.6, linewidth=0)
@df parmdf3 dotplot!([3], :logρ, color = StatsPlots.palette(:Dark2)[6], label="")
scatter!([mean(parmdf.logρ), mean(parmdf2.logρ), mean(parmdf3.logρ)], yerror=[std(parmdf.logρ)/sqrt(nrow(parmdf)), std(parmdf2.logρ)/sqrt(nrow(parmdf2)), std(parmdf3.logρ)/sqrt(nrow(parmdf3))], color=:black, markersize=8, markerstrokewidth=3, label="")


## Epsilon
@df parmdf violin([1], :ϵ, color = StatsPlots.palette(:Dark2)[4], label="", alpha=0.6, linewidth=0)
@df parmdf dotplot!([1], :ϵ, color = StatsPlots.palette(:Dark2)[4], label="")
@df parmdf2 violin!([2], :ϵ, color = StatsPlots.palette(:Dark2)[3],  label="", alpha=0.6, linewidth=0)
@df parmdf2 dotplot!([2], :ϵ, color = StatsPlots.palette(:Dark2)[3], label="",  tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, size=(500, 500), dpi=300)
@df parmdf3 violin!([3], :ϵ, color = StatsPlots.palette(:Dark2)[5], ylabel="ϵ",  xticks=1:3, xlabel="Context", label="", alpha=0.6, linewidth=0)
@df parmdf3 dotplot!([3], :ϵ, color = StatsPlots.palette(:Dark2)[5], label="")

## Omega
bar(mean.([ω[:,1], ω2[:,1], ω3[:,1]]), color=StatsPlots.palette(:Dark2)[[1, 4, 6]], alpha= 0.6, linewidth=0, label="", xticks=1:3, xaxis=:false, ylims=(-2.5, 1.8), yticks=-3:1:2, ylabel="ω₂")
dotplot!([1 2 3], [ω[:,1], ω2[:,1], ω3[:,1]], color=:grey, bar_width=0.3, msw=0, markersize=3,label="", tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, size=(500, 500), dpi=300)

scatter!(mean.([ω[:,1], ω2[:,1], ω3[:,1]]), yerror=std.([ω[:,1], ω2[:,1], ω3[:,1]])./sqrt(nrow(parmdf)), color=:black, markersize=0, markerstrokewidth=10, label="")
#savefig("WMM/Figures/model_based/expe2_omega1.png")
##
bar(mean.([ω[:,2], ω2[:,2], ω3[:,2]]), color=StatsPlots.palette(:Dark2)[[1, 4, 6]], label="", xticks=1:3, xaxis=:false, ylims=(-2.5, 0.8), yticks=-6:1:2,ylabel="ω₁", alpha=0.6, linewidth=0)
dotplot!([1 2 3], [ω[:,2], ω2[:,2], ω3[:,2]], color=:grey, bar_width=0.3, msw=0, markersize=3, label="", tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, size=(500, 500), dpi=300)

scatter!(mean.([ω[:,2], ω2[:,2], ω3[:,2]]), yerror=std.([ω[:,2], ω2[:,2], ω3[:,2]])./sqrt(nrow(parmdf)), color=:black, markersize=0, markerstrokewidth=10, label="")
#savefig("WMM/Figures/model_based/expe2_omega2.png")

##
bar(mean.([ω[:,3], ω2[:,3], ω3[:,3]]), color=StatsPlots.palette(:Dark2)[[1, 4, 6]], label="", xticks=1:3, xaxis=:false, ylims=(-3, 1.5), yticks=-6:1:2,ylabel="ω₀", alpha=0.6, linewidth=0)
dotplot!([1 2 3], [ω[:,3], ω2[:,3], ω3[:,3]], color=:grey, bar_width=0.3, msw=0, markersize=3, label="", tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, size=(500, 500), dpi=300)

scatter!(mean.([ω[:,3], ω2[:,3], ω3[:,3]]), yerror=std.([ω[:,3], ω2[:,3], ω3[:,3]])./sqrt(nrow(parmdf)), color=:black, markersize=0, markerstrokewidth=10, label="")

##
bar(mean.([ω[:,4], ω2[:,4], ω3[:,4]]), color=StatsPlots.palette(:Dark2)[[1, 4, 6]], label="", xticks=1:3, xaxis=:false, ylims=(-3, 2), yticks=-6:2:-2,ylabel="ω₁", alpha=0.6, linewidth=0)
dotplot!([1 2 3], [ω[:,4], ω2[:,4], ω3[:,4]], color=:grey, bar_width=0.3, msw=0, markersize=3, label="", tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, size=(500, 500), dpi=300)

scatter!(mean.([ω[:,4], ω2[:,4], ω3[:,4]]), yerror=std.([ω[:,4], ω2[:,4], ω3[:,4]])./sqrt(nrow(parmdf)), color=:black, markersize=0, markerstrokewidth=10, label="")

##
bar(mean.([ω[:,5], ω2[:,5], ω3[:,5]]), color=StatsPlots.palette(:Dark2)[[1, 4, 6]], label="", xticks=1:3, xaxis=:false, ylims=(-1, 5), yticks=-6:2:-2,ylabel="ω₁", alpha=0.6, linewidth=0)
dotplot!([1 2 3], [ω[:,5], ω2[:,5], ω3[:,5]], color=:grey, bar_width=0.3, msw=0, markersize=3, label="", tickfontsize=14, labelfontsize=20, background_color=:transparent, foreground_color=:black, size=(500, 500), dpi=300)

scatter!(mean.([ω[:,5], ω2[:,5], ω3[:,5]]), yerror=std.([ω[:,5], ω2[:,5], ω3[:,5]])./sqrt(nrow(parmdf)), color=:black, markersize=0, markerstrokewidth=10, label="")

