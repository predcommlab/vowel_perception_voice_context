includet(joinpath(@__DIR__, "VowelPerception.jl"))

using MixedModels
using CategoricalArrays
#using StatsModels


df = @chain begin
    VP.data_all()
    @transform begin
        :ffc = VP.scale_ff(:ff)
        :f0c = VP.scale_f0(:f0)
        :morphc = VP.scale_morph(:morph)
    end
    ## include rows from all mixed trials + from all blocked and precurse that are NOT center
    @subset :experiment == "mixed" || :speaker != "center"
    ## select columns
    select([:experiment, :participant, :wordpair, :morphc, :ffc, :f0c, :response_o])
    transform(:participant => categorical => :participant)
    transform(:experiment => categorical => :experiment)
    @aside levels!(_.experiment, ["mixed", "blocked", "precursors"])
end

# first model with default DummyCoding
model = fit(
    MixedModel,
    @formula(response_o ~ experiment * morphc + experiment * ffc + experiment * f0c + (1 + morphc | participant) + (1 + morphc | wordpair)),
    df,
    Bernoulli(),
)

function latexed_table(model)
    io = IOBuffer()
    show(io, MIME"text/markdown"(), model)
    s = String(take!(io))
    s = replace(s, r"^\|:\-.*$"m => "\\hline")
    s = replace(s, r"^\| "m => "")
    s = replace(s, r"\|$"m => "\\\\")
    s = replace(s, "&" => "\\&")
    s = replace(s, "<" => "\\textless")
    s = replace(s, r"σ_(\w*)" => s"$\\sigma_\\textrm{\1}$")
    s = replace(s, r"\|"m => "&")
end

latexed_table(model) |> println

open(VP.reportpath("mixedmodel_shifts_all.txt"), "w") do file
    print(file, model)
end

open(VP.reportpath("mixedmodel_shifts_all_latex.txt"), "w") do file
    print(file, latexed_table(model))
end

# first model with sequential difference coding
model_sequDiffCod = fit(
    MixedModel,
    @formula(response_o ~ experiment * morphc + experiment * ffc + experiment * f0c + (1 + morphc | participant) + (1 + morphc | wordpair)),
    df,
    Bernoulli(),
    contrasts = Dict(:experiment => SeqDiffCoding())
)

latexed_table(model) |> println

open(VP.reportpath("mixedmodel_shifts_sequDiffCod_all.txt"), "w") do file
    print(file, model_sequDiffCod)
end

open(VP.reportpath("mixedmodel_shifts_sequDiffCod_all_latex.txt"), "w") do file
    print(file, latexed_table(model_sequDiffCod))
end


## model of center speaker in blocked only

df_center_blocked = @chain begin
    VP.data_all()
    @subset :speaker == "center" && :experiment == "blocked"
    @transform AsTable = VP.contextspeaker_to_ff_f0(:context_speaker)
    @transform begin
        :ffc = VP.scale_ff(:ff)
        :f0c = VP.scale_f0(:f0)
        :morphc = VP.scale_morph(:morph)
    end
    select([:experiment, :participant, :wordpair, :morphc, :ffc, :f0c, :response_o])
    transform(:participant => categorical => :participant)
end

model_center_blocked = fit(
    MixedModel,
    @formula(response_o ~ morphc + ffc + f0c + (1 + morphc | participant) + (1 + morphc | wordpair)),
    df_center_blocked,
    Bernoulli(),
)

model_center_blocked

open(VP.reportpath("mixedmodel_shifts_center_blocked.txt"), "w") do file
    print(file, model_center_blocked)
end

open(VP.reportpath("mixedmodel_shifts_center_blocked_latex.txt"), "w") do file
    print(file, latexed_table(model_center_blocked))
end

## model of center speaker after precursors only

df_center_precursors = @chain begin
    VP.data_all()
    @subset :speaker == "center" && :experiment == "precursors"
    @transform AsTable = VP.contextspeaker_to_ff_f0(:context_speaker)
    @transform begin
        :ffc = VP.scale_ff(:ff)
        :f0c = VP.scale_f0(:f0)
    end
    select([:experiment, :participant, :wordpair, :ffc, :f0c, :response_o])
    transform(:participant => categorical => :participant)
end

model_center_precursors = fit(
    MixedModel,
    @formula(response_o ~ ffc + f0c + (1 | participant) + (1 | wordpair)),
    df_center_precursors,
    Bernoulli(),
)

model_center_precursors

open(VP.reportpath("mixedmodel_shifts_center_precursors.txt"), "w") do file
    print(file, model_center_precursors)
end

open(VP.reportpath("mixedmodel_shifts_center_precursors_latex.txt"), "w") do file
    print(file, latexed_table(model_center_precursors))
end

## compare center speaker in blocked + precurser context
df_center_mixed_blocked_precursors = @chain begin
    VP.data_all()
    @subset :speaker == "center"  && :morph == 0.5 && :experiment != "mixed"
    @transform AsTable = VP.contextspeaker_to_ff_f0(:context_speaker)
    @transform begin
        :ffc = VP.scale_ff(:ff)
        :f0c = VP.scale_f0(:f0)
    end
    transform(:experiment => categorical => :experiment)

    select([:experiment, :participant, :wordpair, :ffc, :f0c, :response_o])
    transform(:participant => categorical => :participant)
end

model_center_blocked_precursors = fit(
    MixedModel,
    @formula(response_o ~ experiment*ffc + experiment*f0c + (1 | participant) + (1 | wordpair)),
    df_center_mixed_blocked_precursors,
    Bernoulli(),
    #contrasts = Dict(:experiment => SeqDiffCoding())
)

model_center_blocked_precursors

open(VP.reportpath("mixedmodel_shifts_center_blocked_precursors.txt"), "w") do file
    print(file, model_center_blocked_precursors)
end

open(VP.reportpath("mixedmodel_shifts_center_blocked_precursors_latex.txt"), "w") do file
    print(file, latexed_table(model_center_blocked_precursors))
end