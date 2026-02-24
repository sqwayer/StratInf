## Data
folder = "data/HMM_data"
datadf = load_HMM_data(folder)

# Recode perseverations before switches for models 
pres_in_block!(datadf, blockId = :blockNum) 
datadf[datadf.negpresInBlock .< -3,:correctOrPersevProb_cfql] = datadf[datadf.negpresInBlock .< -3,:persevProb_cfql]

datadf[datadf.negpresInBlock .< -3,:correctOrPersevProb_full_adaptive_cfql2] = datadf[datadf.negpresInBlock .< -3,:persevProb_full_adaptive_cfql2]

datadf[datadf.negpresInBlock .< -3,:correctOrPersevProb_probe_LRRel] = datadf[datadf.negpresInBlock .< -3,:persevProb_probe_LRRel]

datadf[datadf.negpresInBlock .< -3,:correctOrPersevProb_SI_MultVol_SampleAction] = datadf[datadf.negpresInBlock .< -3,:persevProb_SI_MultVol_SampleAction]

df1 = process_HMM_data!(datadf)
df1.negpresInBlock .+= 1 # So the last trial is 0 and not -1

## Compare models predictions around HMM switch points 

summarySubAfter = combine(groupby(df1, [:presInBlock, :condition, :subject]), 
    :persevProb_full_adaptive_cfql2 => mean => :QL_pmodel_mean,     :persevProb_SI_MultVol_SampleAction => mean => :SI_pmodel_mean, 
    :persevProb_probe_LRRel => mean => :PL_pmodel_mean)
# summarySubAfter = combine(groupby(df1, [:presInBlock, :condition, :subject]), 
#     :L_full_adaptive_cfql2 => mean => :QL_pmodel_mean,     :H_SI_MultVol_SampleAction => (x -> mean(x ./ log(28))) => :SI_pmodel_mean, 
#     :H_probe_LRRel => (x -> mean(x ./ log(3))) => :PL_pmodel_mean)

summaryAfter = combine(groupby(summarySubAfter, [:presInBlock, :condition]), 
:QL_pmodel_mean => mean => :QL_pmodel_mean, 
:QL_pmodel_mean => sem => :QL_pmodel_sem, 
:SI_pmodel_mean => mean => :SI_pmodel_mean, 
:SI_pmodel_mean => sem => :SI_pmodel_sem, 
:PL_pmodel_mean => mean => :PL_pmodel_mean, 
:PL_pmodel_mean => sem => :PL_pmodel_sem)

summaryAfter = summaryAfter[.!isnan.(summaryAfter.SI_pmodel_sem),:]

summarySubBefore = combine(groupby(df1, [:negpresInBlock, :nextCondition, :subject]), 
    :correctOrPersevProb_full_adaptive_cfql2 => mean => :QL_pmodel_mean, :correctOrPersevProb_SI_MultVol_SampleAction => mean => :SI_pmodel_mean, :correctOrPersevProb_probe_LRRel => mean => :PL_pmodel_mean)
# summarySubBefore = combine(groupby(df1, [:negpresInBlock, :nextCondition, :subject]), 
#     :L_full_adaptive_cfql2 => mean => :QL_pmodel_mean, 
#     :H_SI_MultVol_SampleAction => (x -> mean(x ./ log(28))) => :SI_pmodel_mean, 
#     :H_probe_LRRel => (x -> mean(x ./ log(3))) => :PL_pmodel_mean)

summaryBefore = combine(groupby(summarySubBefore, [:negpresInBlock, :nextCondition]), :QL_pmodel_mean => mean => :QL_pmodel_mean, 
:QL_pmodel_mean => sem => :QL_pmodel_sem, 
:SI_pmodel_mean => mean => :SI_pmodel_mean, 
:SI_pmodel_mean => sem => :SI_pmodel_sem, 
:PL_pmodel_mean => mean => :PL_pmodel_mean, 
:PL_pmodel_mean => sem => :PL_pmodel_sem)

summaryBefore = summaryBefore[.!isnan.(summaryBefore.SI_pmodel_sem),:]

rename!(summaryBefore, :negpresInBlock => :presInBlock)
rename!(summaryBefore, :nextCondition => :condition)
grpstats = vcat(summaryBefore, summaryAfter)

grpstats = grpstats[0 .< grpstats.condition,: ]

grp_plot_hmm(grpstats, "SI_pmodel", xlims=(-3, 9), ylims=(0,1), size=(500, 500), dpi = 300, background_color = :transparent)





