module VowelPerception
    using CSV
    using DataFrames
    using Chain
    using CairoMakie

    datadir(args...) = joinpath(@__DIR__, "..", "data", args...)
    data_mixed() = CSV.read(datadir("mixed_valid.csv"), DataFrame)
    data_blocked() = CSV.read(datadir("blocked_valid.csv"), DataFrame)
    data_precursors() = CSV.read(datadir("precursors_valid.csv"), DataFrame)
    data_all() = CSV.read(datadir("all_valid.csv"), DataFrame)

    function base_theme()
        Theme(
            font = normal_font(),
            Axis = (
                xgridvisible = false,
                ygridvisible = false,
                topspinevisible = false,
                rightspinevisible = false,
                xlabelfont = medium_font(),
                ylabelfont = medium_font(),
                titlefont = medium_font(),
            ),
            Legend = (
                titlefont = medium_font(),
                framevisible = false,
            )
        )
    end

    function facet_index(speaker)
        Dict(
            "center" => (2, 2),
            "high_ff" => (2, 3),
            "high_f0" => (1, 2),
            "low_f0" => (3, 2),
            "low_ff" => (2, 1),
        )[speaker]
    end

    function speaker_pretty(speaker)
        Dict(
            missing => "Center",
            "center" => "Center",
            "high_ff" => "High FF",
            "high_f0" => "High F0",
            "low_f0" => "Low F0",
            "low_ff" => "Low FF",
        )[speaker]
    end

    # # blue / red dark / light
    function contextcolor(speaker)
        Dict(
            "center" => "gray30",
            missing => "gray30",
            "high_ff" => "#f33f51",
            "high_f0" => "#3196aa",
            "low_f0" => "#3a5d83",
            "low_ff" => "#9b3156",
        )[speaker]
    end

    function exp_linestyle(experiment)
        Dict(
            "mixed" => :solid,
            "blocked" => :dash,
            "precursors" => :dot,
        )[experiment]
    end

    # fonts may need to be downloaded 
    # should work fine without specifying font even with provided code 
    normal_font() = raw"C:\Users\~\AppData\Local\Microsoft\Windows\Fonts\HelveticaNeueLTStd-Lt.otf"
    medium_font() = raw"C:\Users\~\AppData\Local\Microsoft\Windows\Fonts\HelveticaNeueLTStd-Md.otf"
    bold_font() = raw"C:\Users\~\AppData\Local\Microsoft\Windows\Fonts\HelveticaNeueLTStd-Bd.otf"

    scale_morph(m) = 4m - 2
    unscale_morph(s) = (s + 2) / 4
    scale_f0(f0) = f0 == 1 ? 0.0 : f0 > 1 ? 1.0 : -1.0
    scale_ff(ff) = ff == 1 ? 0.0 : ff > 1 ? 1.0 : -1.0

    reportpath(args...) = joinpath(@__DIR__, "..", "reports", args...)

    function contextspeaker_to_ff_f0(s)
        Dict(
            missing => (f0 = 1.0, ff = 1.0),
            "center" => (f0 = 1.0, ff = 1.0),
            "high_f0" => (f0 = 1.24, ff = 1.0),
            "high_ff" => (f0 = 1.0, ff = 1.146),
            "low_f0" => (f0 = 0.806, ff = 1.0),
            "low_ff" => (f0 = 1.0, ff = 0.873),
        )[s]
    end

    hz_to_bark(hz) = (26.81hz / (1960 + hz)) - 0.53
    bark_to_hz(b) = -(490 * (53 + 100 * b)) / (-657 + 25 * b)
end
