using StatsPlots, ColorBrewer
function grp_stats(df)

    # Remove condition 0 (first block)
    df = df[df.condition .> 0,:]

    # Compute negative trials indices 
    pres_in_block!(df)
    #nextCondition!(df)
    df.negpresInBlock .+= 1 # So the last trial is 0 and not -1

    # Compute stats per subject
    gdf = groupby(df, [:subject, :condition, :isStable, :presInBlock])
    persub_forward = combine(gdf, :correct => mean => :correct, :persev => mean => :persev, :explo => mean => :explo, :zrt => mean => :rt)

    gdf = groupby(df, [:subject, :nextCondition, :stableAS, :negpresInBlock])
    persub_backward = combine(gdf, :correct => mean => :correct, :persev => mean => :persev, :explo => mean => :explo, :zrt => mean => :rt)

    # Compute stats for the whole group
    gdf = groupby(persub_forward, [:condition, :isStable, :presInBlock])
    grpstats_forward = combine(gdf, :correct => mean, :correct => sem, :persev => mean, :persev => sem, :explo => mean, :explo => sem, :rt => mean, :rt => sem)
    rename!(grpstats_forward, :isStable => :stableAS)

    gdf = groupby(persub_backward, [:nextCondition, :stableAS, :negpresInBlock])
    grpstats_backward = combine(gdf, :correct => mean, :correct => sem, :persev => mean, :persev => sem, :explo => mean, :explo => sem, :rt => mean, :rt => sem)

    rename!(grpstats_backward, :negpresInBlock => :presInBlock)
    rename!(grpstats_backward, :nextCondition => :condition)
    grpstats_backward.persev_mean .= grpstats_backward.correct_mean
    grpstats_backward.persev_sem .= grpstats_backward.correct_sem
    grpstats = vcat(grpstats_backward, grpstats_forward)
    grpstats = grpstats[.!isnan.(grpstats.correct_sem),:]
    
    return grpstats
end

function grp_stats_hmm(df)
    df.negpresInBlock .+= 1 # So the last trial is 0 and not -1
    # After switch
    sdf = groupby(df, [:subject, :presInBlock, :condition])
    mean_per_sub = combine(sdf, :persev => mean => :persev, :correct => mean => :correct, :zrt => mean=> :rt)
    
    sdf = groupby(mean_per_sub, [:presInBlock, :condition])
    groupstat_after = combine(sdf, :persev => mean, :persev => sem, :correct => mean, :correct => sem, :rt => mean, :rt => sem)
    groupstat_after = groupstat_after[.!isnan.(groupstat_after.persev_sem),:] 

    # Before switch
    sdf = groupby(df, [:subject, :negpresInBlock, :nextCondition])

    mean_per_sub = combine(sdf, :correct_or_persev => mean => :persev,:correct => mean => :correct, :zrt => mean => :rt)
    sdf = groupby(mean_per_sub, [:negpresInBlock, :nextCondition])
    groupstat_before = combine(sdf, :persev => mean, :persev => sem, :correct => mean, :correct => sem, :rt => mean, :rt => sem)
    groupstat_before = groupstat_before[.!isnan.(groupstat_before.persev_sem),:] 
    rename!(groupstat_before, :negpresInBlock => :presInBlock)
    rename!(groupstat_before, :nextCondition => :condition)
    grpstats = vcat(groupstat_before, groupstat_after)
    return grpstats
end


function grp_plot(grpstats, stat="correct", stable="all", cols=eachindex(unique(grpstats.condition)); kwargs...)
    cpal = StatsPlots.palette(:Dark2)[cols]
    if stable == "all"
        gg = groupby(grpstats, [:condition, :presInBlock])
        grpstats = combine(gg, 
            :correct_mean => mean => :correct_mean, 
            :correct_sem => sem => :correct_sem, 
            :persev_mean => mean => :persev_mean, 
            :persev_sem => sem => :persev_sem, 
            :explo_mean => mean => :explo_mean, 
            :explo_sem => sem => :explo_sem, 
            :rt_mean => mean => :rt_mean, 
            :rt_sem => sem => :rt_sem)
        grpstats = grpstats[.!isnan.(grpstats.correct_sem),:]
    elseif stable == "true"
        grpstats = grpstats[grpstats.stableAS,:]
    elseif stable == "false"
        grpstats = grpstats[.!grpstats.stableAS,:]
    end

    before_df = grpstats[grpstats.presInBlock .<= 0,:]
    after_df = grpstats[grpstats.presInBlock .> 0,:]
    if stat == "correct"
        @df before_df plot(:presInBlock, :correct_mean; palette = cpal, group=:condition, ribbon=:correct_sem, linewidth=3, kwargs...)
        @df after_df plot!(:presInBlock, :correct_mean; palette = cpal, group=:condition, ribbon=:correct_sem, linewidth=3,kwargs..., label="")
    elseif stat == "persev"
        @df grpstats plot(:presInBlock, :persev_mean;  palette = cpal, group=:condition, ribbon=:persev_sem, linewidth=3, kwargs...)
    elseif stat == "explo"
        @df before_df plot(:presInBlock, :explo_mean;  palette = cpal, group=:condition, ribbon=:explo_sem, linewidth=3, kwargs...)
        @df after_df plot!(:presInBlock, :explo_mean;  palette = cpal, group=:condition, ribbon=:explo_sem, linewidth=3, kwargs..., label="")
    elseif stat == "rt"
        @df grpstats plot(:presInBlock, :rt_mean;  palette = cpal, group=:condition, ribbon=:rt_sem, linewidth=3, kwargs...)
    end
end

function grp_plot_hmm(grpstats, stat="persev", cols=eachindex(unique(grpstats.condition)); kwargs...)
    cpal = StatsPlots.palette(:Dark2)[cols]

    before_df = grpstats[grpstats.presInBlock .<= 0,:]
    #before_df.presInBlock = before_df.presInBlock .+ 0.9
    after_df = grpstats[grpstats.presInBlock .> 0,:]
    #after_df.presInBlock = after_df.presInBlock .- 0.9
    if stat == "explo"
        @df before_df plot(:presInBlock, :explo_mean; palette = cpal, group=:condition, ribbon=:explo_sem, linewidth=3, kwargs...)
        @df after_df plot!(:presInBlock, :correct_mean; palette = cpal, group=:condition, ribbon=:explo_sem, linewidth=3, kwargs..., label="")

    elseif stat == "persev"
        @df before_df plot(:presInBlock, :persev_mean;  palette = cpal, group=:condition, ribbon=:persev_sem, linewidth=3, kwargs...)
        @df after_df plot!(:presInBlock, :persev_mean;  palette = cpal, group=:condition, ribbon=:persev_sem, linewidth=3, kwargs..., label="")
        #plot!([-0.5, 0.5, 0.5, -0.5],  [0, 0, 1, 1], seriestype=:shape, alpha = 0.6, linewidth = 0, color=:grey, fillstyle = :/, label = "")
    elseif stat == "correct"
        @df before_df plot(:presInBlock, :correct_mean;  palette = cpal, group=:condition, ribbon=:correct_sem, linewidth=3, kwargs...)
        @df after_df plot!(:presInBlock, :correct_mean;  palette = cpal, group=:condition, ribbon=:correct_sem,linewidth=3, kwargs..., label="")
    elseif stat == "rt"
        @df before_df plot(:presInBlock, :rt_mean;  palette = cpal, group=:condition, ribbon=:rt_sem, linewidth=3, kwargs...)
        @df after_df plot!(:presInBlock, :rt_mean;  palette = cpal, group=:condition, ribbon=:rt_sem, linewidth=3, kwargs..., label="")
        #plot!([-0.5, 0.5, 0.5, -0.5], [-1, -1, 1, 1], seriestype=:shape, alpha = 0.6, linewidth = 0, color=:grey, fillstyle = :/, label = "")

    end
end


function grp_plot!(grpstats, stat="correct", stable="all", cols=eachindex(unique(grpstats.condition)), errorstyle = "ribbon";kwargs...)
    cpal = StatsPlots.palette(:Dark2)[cols]

    if stable == "all"
        gg = groupby(grpstats, [:condition, :presInBlock])
        grpstats = combine(gg, 
            :correct_mean => mean => :correct_mean, 
            :correct_sem => sem => :correct_sem, 
            :persev_mean => mean => :persev_mean, 
            :persev_sem => sem => :persev_sem, 
            :explo_mean => mean => :explo_mean, 
            :explo_sem => sem => :explo_sem, 
            :rt_mean => mean => :rt_mean, 
            :rt_sem => sem => :rt_sem)
        grpstats = grpstats[.!isnan.(grpstats.correct_sem),:]
    elseif stable == "true"
        grpstats = grpstats[grpstats.stableAS,:]
    elseif stable == "false"
        grpstats = grpstats[.!grpstats.stableAS,:]
    end


    before_df = grpstats[grpstats.presInBlock .<= 0,:]
    after_df = grpstats[grpstats.presInBlock .> 0,:]
    if errorstyle == "ribbon"
        if stat == "correct"
            @df before_df plot!(:presInBlock, :correct_mean; palette = cpal, group=:condition, ribbon=:correct_sem, linewidth=3, msw=3, msc = cpal', kwargs...)
            @df after_df plot!(:presInBlock, :correct_mean; palette = cpal, group=:condition, ribbon=:correct_sem, linewidth=3, msw=3, msc = cpal', kwargs..., label="")
        elseif stat == "persev"
            @df grpstats plot!(:presInBlock, :persev_mean;  palette = cpal, group=:condition, yerror=:persev_sem, linewidth=3, msw=3, msc = cpal, kwargs...)
        elseif stat == "explo"
            @df before_df plot!(:presInBlock, :explo_mean;  palette = cpal, group=:condition, yerror=:explo_sem, linewidth=3, msw=3, msc = cpal, kwargs...)
            @df after_df plot!(:presInBlock, :explo_mean;  palette = cpal, group=:condition, yerror=:explo_sem, linewidth=3, msw=3, msc = cpal, kwargs..., label="")
        elseif stat == "rt"
            @df grpstats plot!(:presInBlock, :rt_mean;  palette = cpal, group=:condition, yerror=:rt_sem, linewidth=3, msw=3, msc = cpal, kwargs...)
        end
    elseif errorstyle == "bar"
        if stat == "correct"
            @df before_df plot!(:presInBlock, :correct_mean; palette = cpal, group=:condition, yerror=:correct_sem, linewidth=3, kwargs...)
            @df after_df plot!(:presInBlock, :correct_mean; palette = cpal, group=:condition, yerror=:correct_sem, linewidth=3,kwargs..., label="")
            # Replot markers
            @df before_df scatter!(:presInBlock, :correct_mean; palette = cpal, group=:condition, yerror=:correct_sem,markershape=:circle,  markercolor=:white, markersize=4, msw=2, msc = cpal', kwargs...)
            @df after_df scatter!(:presInBlock, :correct_mean; palette = cpal, group=:condition, yerror=:correct_sem, markershape=:circle, markercolor=:white, markersize=4, msw=2, msc = cpal', kwargs..., label="")
        elseif stat == "persev"
            @df grpstats plot!(:presInBlock, :persev_mean;  palette = cpal, group=:condition, yerror=:persev_sem, linewidth=3, msw=3, msc = cpal, kwargs...)
        elseif stat == "explo"
            @df before_df plot!(:presInBlock, :explo_mean;  palette = cpal, group=:condition, yerror=:explo_sem, linewidth=3, msw=3, msc = cpal, kwargs...)
            @df after_df plot!(:presInBlock, :explo_mean;  palette = cpal, group=:condition, yerror=:explo_sem, linewidth=3, msw=3, msc = cpal, kwargs..., label="")
        elseif stat == "rt"
            @df grpstats plot!(:presInBlock, :rt_mean;  palette = cpal, group=:condition, yerror=:rt_sem, linewidth=3, msw=3, msc = cpal, kwargs...)
        end
    end
end
function boxplot2d(X, Y; kwargs...)
    qX = quantile.([X], [0.05, 0.25, 0.5, 0.75, 0.95])
    qY = quantile.([Y], [0.05, 0.25, 0.5, 0.75, 0.95])
    pl = plot([qY[2], qY[4], qY[4], qY[2]], [qX[2], qX[2], qX[4], qX[4]];kwargs..., seriestype=:shape)
    plot!([qY[1], qY[5]], [qX[3], qX[3]]; color=:black, kwargs..., label="")
    plot!([qY[3], qY[3]], [qX[1], qX[5]]; color=:black, kwargs..., label="")
end

function boxplot2d!(X, Y; kwargs...)
    qX = quantile.([X], [0.05, 0.25, 0.5, 0.75, 0.95])
    qY = quantile.([Y], [0.05, 0.25, 0.5, 0.75, 0.95])
    plot!([qY[2], qY[4], qY[4], qY[2]], [qX[2], qX[2], qX[4], qX[4]];kwargs..., seriestype=:shape)
    plot!([qY[1], qY[5]], [qX[3], qX[3]]; color=:black, kwargs..., label="")
    plot!([qY[3], qY[3]], [qX[1], qX[5]]; color=:black, kwargs..., label="")
end