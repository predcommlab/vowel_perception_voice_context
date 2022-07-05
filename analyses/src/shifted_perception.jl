includet(joinpath(@__DIR__, "VowelPerception.jl"))

const VP = VowelPerception

using CairoMakie
using DataFrameMacros
using Chain
using DataFrames
using StatsBase

function speaker_axis_settings!(axdict)
    for (speaker, ax) in axdict
        if speaker != "low_ff"
            hideydecorations!(ax)
        end
        if speaker != "low_f0"
            hidexdecorations!(ax)
        end
    end
end

function shift_plot()
    set_theme!(VP.base_theme())

    df_mixed = VP.data_mixed()
    df_blocked = VP.data_blocked()
    df_precursors = VP.data_precursors()

    f = Figure(resolution = (800, 800))

    g1 = f[1, 1] = GridLayout()
    mixed_axes = Dict()
    for gdf in groupby(df_mixed, :speaker)
        speaker = gdf.speaker[1]
        facet = VP.facet_index(speaker)
        mixed_axes[speaker] = Axis(
            g1[facet...],
            xlabel = "Morph level",
            ylabel = "Response",
            limits = (0, 1, 0, 1),
            xticks = (0:0.5:1, ["uː", "¹/₂", "oː"]),
            yticks = (0:0.5:1, ["uː", "¹/₂", "oː"]),
        )

        responses = @combine(
            groupby(gdf, [:participant, :morph], sort = true),
            :p_o = mean(:response_o)
        )

        for gdf2 in groupby(responses, :participant)
            lines!(gdf2.morph, gdf2.p_o, color = :gray90, linewidth = 0.75)
        end

        responses_mean = @combine(
            groupby(responses, :morph),
            :p_o = mean(:p_o)
        )

        l = scatterlines!(
            responses_mean.morph,
            responses_mean.p_o,
            color = VP.contextcolor("center"),
            markercolor = :white,
            strokecolor = VP.contextcolor("center"),
            strokewidth = 0.75,
            markersize = 6,
            linestyle = VP.exp_linestyle("mixed"),
        )
        translate!(l, 0, 0, 10)

    end



    g2 = f[1, 2] = GridLayout()
    blocked_axes = Dict()
    for gdf in groupby(df_blocked, :speaker)
        speaker = gdf.speaker[1]
        facet = VP.facet_index(speaker)
        blocked_axes[speaker] = Axis(
            g2[facet...],
            xlabel = "Morph level",
            ylabel = "Response",
            limits = (0, 1, 0, 1),
            xticks = (0:0.5:1, ["uː", "¹/₂", "oː"]),
            yticks = (0:0.5:1, ["uː", "¹/₂", "oː"]),
        )

        for ggdf in groupby(gdf, :context_speaker)
            color = VP.contextcolor(ggdf.context_speaker[1])

            responses = @combine(
                groupby(ggdf, [:participant, :morph], sort = true),
                :p_o = mean(:response_o)
            )

            for gdf2 in groupby(responses, :participant)
                lines!(gdf2.morph, gdf2.p_o, color = (color, 0.1), linewidth = 0.75)
            end

            responses_mean = @combine(
                groupby(responses, :morph),
                :p_o = mean(:p_o)
            )

            l = scatterlines!(
                responses_mean.morph,
                responses_mean.p_o,
                color = color,
                markercolor = :white,
                strokecolor = color,
                strokewidth = 0.75,
                markersize = 6,
                linestyle = VP.exp_linestyle("blocked"),
            )
            translate!(l, 0, 0, 10)
        end

    end

    g3 = f[2, 1] = GridLayout()
    precursor_axes = Dict()
    for gdf in groupby(df_precursors, :speaker)
        speaker = gdf.speaker[1]
        facet = VP.facet_index(speaker)
        precursor_axes[speaker] = Axis(
            g3[facet...],
            xlabel = "Morph level",
            ylabel = "Response",
            limits = (0, 1, 0, 1),
            xticks = (0:0.5:1, ["uː", "¹/₂", "oː"]),
            yticks = (0:0.5:1, ["uː", "¹/₂", "oː"]),
        )

        for ggdf in groupby(gdf, :context_speaker)
            color = VP.contextcolor(ggdf.context_speaker[1])

            responses = @combine(
                groupby(ggdf, [:participant, :morph], sort = true),
                :p_o = mean(:response_o)
            )

            for gdf2 in groupby(responses, :participant)
                lines!(gdf2.morph, gdf2.p_o, color = (color, 0.1), linewidth = 0.75)
            end

            responses_mean = @combine(
                groupby(responses, :morph),
                :p_o = mean(:p_o)
            )

            l = scatterlines!(
                responses_mean.morph,
                responses_mean.p_o,
                color = color,
                markercolor = :white,
                strokecolor = color,
                strokewidth = 0.75,
                markersize = 6,
                linestyle = VP.exp_linestyle("precursors"),
            )
            translate!(l, 0, 0, 10)
        end

    end

    g4 = f[2, 2] = GridLayout()
    combined_axes = Dict()
    for gdf in groupby(VP.data_all(), [:experiment, :speaker])
        speaker = gdf.speaker[1]
        facet = VP.facet_index(speaker)
        ax = get!(combined_axes, speaker) do
            Axis(
                g4[facet...],
                xlabel = "Morph level",
                ylabel = "Response",
                limits = (0, 1, 0, 1),
                xticks = (0:0.5:1, ["uː", "¹/₂", "oː"]),
                yticks = (0:0.5:1, ["uː", "¹/₂", "oː"]),
            )
        end
        current_axis!(ax)

        for ggdf in groupby(gdf, :context_speaker)
            context_speaker = ggdf.context_speaker[1]
            color = VP.contextcolor(context_speaker)

            responses_mean = @combine(
                groupby(ggdf, :morph, sort = true),
                :p_o = mean(:response_o)
            )

            l = scatterlines!(
                responses_mean.morph,
                responses_mean.p_o,
                color = color,
                markercolor = :white,
                strokecolor = color,
                strokewidth = 0.75,
                markersize = 6,
                linestyle = VP.exp_linestyle(gdf.experiment[1]),
            )
            translate!(l, 0, 0, 10)
        end
    end

    speaker_axis_settings!(blocked_axes)
    speaker_axis_settings!(mixed_axes)
    speaker_axis_settings!(precursor_axes)
    speaker_axis_settings!(combined_axes)

    for g in [g1, g2, g3, g4]
        colgap!(g, 10)
        rowgap!(g, 10)
    end

    glegend = f[1:2, 1:2] = GridLayout(
        tellheight = false,
        tellwidth = false,
    )
    Legend(
        glegend[1, 1],
        [
            LineElement(color = VP.contextcolor("low_ff")),
            LineElement(color = VP.contextcolor("high_ff")),
            LineElement(color = VP.contextcolor("low_f0")),
            LineElement(color = VP.contextcolor("high_f0")),
            LineElement(color = VP.contextcolor("center")),
        ],
        ["Low FF", "High FF", "Low F0", "High F0", "No context"],
        "Context",
        framevisible = false,
        padding = (0, 0, 0, 0),
    )
    Legend(
        glegend[1, 2],
        [
            LineElement(color = :black, linestyle = VP.exp_linestyle("mixed")),
            LineElement(color = :black, linestyle = VP.exp_linestyle("blocked")),
            LineElement(color = :black, linestyle = VP.exp_linestyle("precursors")),
        ],
        ["Mixed", "Blocked", "Precursors"],
        "Condition",
        framevisible = false,
        padding = (0, 0, 0, 0),
    )
    colgap!(glegend, 5)

    ga = f[1, 1] = GridLayout(tellheight = false, tellwidth = false, halign = :left, valign = :top)
    Label(ga[1, 1], "A", font = VP.bold_font(), textsize = 36,
        halign = :left)
    Label(ga[2, 1], "Exp. 1\nMixed", halign = :left)

    gb = f[1, 2] = GridLayout(tellheight = false, tellwidth = false, halign = :right, valign = :top)
    Label(gb[1, 1], "B", font = VP.bold_font(), textsize = 36,
        halign = :right, valign = :top)
    lb = Label(gb[2, 1], "Exp. 2\nBlocked", halign = :right)
    lb.elements[:text].justification = :right

    gc = f[2, 1] = GridLayout(tellheight = false, tellwidth = false, halign = :left, valign = :top)
    Label(gc[1, 1], "C", font = VP.bold_font(), textsize = 36,
        halign = :left, valign = :top)
    lc = Label(gc[2, 1], "Exp. 3\nPrecursors", halign = :left)
    lc.elements[:text].justification = :left

    gd = f[2, 2] = GridLayout(tellheight = false, tellwidth = false, halign = :right, valign = :top)
    Label(gd[1, 1], "D", font = VP.bold_font(), textsize = 36,
        halign = :right, valign = :top)
    Label(gd[2, 1], "Comparison", halign = :right)

    rowgap!.([ga, gb, gc, gd], 0)

    set_theme!()
    f
end

@chain begin
    shift_plot()
    save("figures/perceptual_shifts.pdf", _)
end

