using Revise
includet(joinpath(@__DIR__, "VowelPerception.jl"))

const VP = VowelPerception

using CairoMakie
using DataFrameMacros
using Chain
using DataFrames
using StatsBase

binsize = 52

# enter 1 for precursor data
function plot_binned_temporal_development(arg)
    binned = @chain begin
        if arg ==1
            VP.data_precursors()
        else
            VP.data_blocked()
        end

        groupby([:participant, :context_speaker])
        # @c operates by column 
        @transform :trial_in_block = @c :trial_number .- minimum(:trial_number) .+ 1
        @aside @show maximum(_.trial_in_block)
        @transform :trial_bin = fld1(:trial_in_block, binsize)
        # speaker center, high FF, high F0
        groupby([:context_speaker, :trial_bin, :speaker, :morph])
        @combine begin
            :mean_po = mean(:response_o)
            :std_po = sem(:response_o)
            :var_po = var(:response_o)
            :n_trials = length(:response_o)
        end
        # only take trials without center speaker
        @subset :context_speaker == :speaker

        sort([:speaker, :trial_bin, :morph])
        # sort([:context_speaker, :trial_bin, :morph])
    end


    centered = @chain begin
        if arg ==1
            VP.data_precursors()
        else
            VP.data_blocked()
        end

        groupby([:participant, :context_speaker])
        # @c operates by column 
        @transform :trial_in_block = @c :trial_number .- minimum(:trial_number) .+ 1
        # @aside @show maximum(_.trial_in_block)
        @transform :trial_bin = fld1(:trial_in_block, binsize)
        # speaker center, high FF, high F0
        groupby([:context_speaker, :trial_bin, :speaker, :morph])
        @combine begin
            :mean_po = mean(:response_o)
            :std_po = sem(:response_o)
            :var_po = var(:response_o)
            :n_trials = length(:response_o)
        end
        # only take trials from center speaker in context
        @subset :context_speaker != :speaker
        # filter 0.5 morph
        @subset :morph == 0.5

        # sort([:speaker, :trial_bin, :morph])
        sort([:context_speaker, :trial_bin, :morph])
    end

    bins = sort(unique(binned.trial_bin))
    binedges = map(bins) do b
        (b-1) * binsize + 1, min(207, b * binsize)
    end
    binlabels = map(binedges) do (lo, hi)
        "$lo - $hi"
    end

    set_theme!(VP.base_theme())
    f = Figure(resolution = (800, 600))

    axes = Axis[]

    # iterate over gdf grouped data frame

    # println(binned)
    # binned_ = subset(binned, :context_speaker => ByRow(==(:speaker))) # --> does not work?
    # println(binned_)

    # binned_ = filter(row -> (row.context_speaker == row.speaker),  binned)
    # println(binned_)


    # for (i, gdf) in enumerate(groupby(
    #     filter(row -> (row.context_speaker == row.speaker),  binned), 
    #     :speaker, sort = true))
    line_style = [:dash, :dot, :solid, :dashdotdot, [0.5, 1.0, 1.5, 2.5]]

    # iterate ofer grouped dataframe with coontext speaker + gdf w center speaker
    for (i, (gdf, gdf_c)) in enumerate(zip(groupby(binned, :context_speaker, sort = true), groupby(centered, :context_speaker, sort = true)))

        ax = Axis(f[fldmod1(i, 2)...],
            title = VP.speaker_pretty(gdf.context_speaker[1])*(" context"),
            xticks = (bins, binlabels),
            yticks = ([0, 0.5, 1.0], ["u", "¹/₂", "o"]),
            xticklabelrotation = pi/4,
            xlabel = "Trial bin",
            ylabel = "Mean p(o)"
        )


        for (j, ggdf) in enumerate(groupby(gdf, :morph))
            c = VP.contextcolor(ggdf.context_speaker[1]) 

            if arg == 1
                # precursor condition has only 3 morph levels
                lines!(ggdf.trial_bin, ggdf.mean_po, color = (c, j/3), linestyle = line_style[j+1])
                # standard error of the mean as light background band
                band!(ggdf.trial_bin,(ggdf.mean_po .- ggdf.std_po), (ggdf.mean_po .+ ggdf.std_po), color = (c, (j/3)*0.2))
            else
                lines!(ggdf.trial_bin, ggdf.mean_po, color = (c, j/5), linestyle = line_style[j])
                band!(ggdf.trial_bin,(ggdf.mean_po .- ggdf.std_po), (ggdf.mean_po .+ ggdf.std_po), color = (c, (j/5)*0.2) )
            end

        end

        # # plot 0.5 morph for center speaker - same for both conditions
        c = VP.contextcolor(gdf_c.speaker[1])
        lines!(gdf_c.trial_bin, gdf_c.mean_po, color = (c, 3/5), linestyle = :solid)
        band!(gdf_c.trial_bin,(gdf_c.mean_po .- gdf_c.std_po), (gdf_c.mean_po .+ gdf_c.std_po), color = (c, (3/5)*0.2))

        ylims!(ax, 0, 1)
        push!(axes, ax)
    end

    axes = permutedims(reshape(axes, (2, 2)))
    hideydecorations!.(axes[:, 2])
    hidexdecorations!.(axes[1, :])

    # Legend(f[1:2, 3],
    #     [LineElement(color = (:black, j/5), linestyle = line_style[j]) for j in 1:5],
    #     ["u", "0.25", "0.5", "0.75", "o"],
    #     "Morph",
    #     framevisible = false,
    # )

    glegend = f[1:2, 3] = GridLayout(
        tellheight = false,
        tellwidth = false,
    )
    colsize!(f.layout, 3, Fixed(150))


    if arg ==1
        Legend(
            glegend[1, 3],
            [LineElement(color = (:black, j/3), linestyle = line_style[j+1]) for j in 1:3],
            ["0.25", "0.5", "0.75"],
            "Morph",
            framevisible = false,
            padding = (0, 0, 0, 0),
        )
    else
        Legend(
            glegend[1, 3],
            [LineElement(color = (:black, j/5), linestyle = line_style[j]) for j in 1:5],
            ["u", "0.25", "0.5", "0.75", "o"],
            "Morph",
            framevisible = false,
            padding = (0, 0, 0, 0),
        )
    end

    c = VP.contextcolor("center") 
    Legend(
        glegend[2, 3],
        [
            LineElement(color = VP.contextcolor("low_ff")),
            LineElement(color = VP.contextcolor("high_ff")),
            LineElement(color = VP.contextcolor("low_f0")),
            LineElement(color = VP.contextcolor("high_f0")),
            LineElement(color = (c, 3/5)),
        ],
        ["Low FF", "High FF", "Low F0", "High F0", "Center"],
        "Speaker",
        framevisible = false,
        padding = (0, 0, 0, 0),
    )
    colgap!(glegend, 5)

    if arg==1
        Label(f[0, 1:3], "Mean responses across time in precursor condition",
        font = VP.medium_font(), textsize = 22, color = :gray40)
    else
        Label(f[0, 1:3], "Mean responses across time in blocked condition",
        font = VP.medium_font(), textsize = 22, color = :gray40
    )
    end

    set_theme!()
    f
end

@chain begin
    plot_binned_temporal_development(1)
    save(joinpath(@__DIR__,"..", "figures/binned_temporal_development_center_precursor.pdf"), _)
end

@chain begin
    plot_binned_temporal_development(2)
    save(joinpath(@__DIR__,"..","figures/binned_temporal_development_center_blocked.pdf"), _)
end

#############


