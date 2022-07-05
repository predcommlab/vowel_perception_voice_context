using Revise
includet(joinpath(@__DIR__, "VowelPerception.jl"))

using MixedModels
using CategoricalArrays
using Chain
using DataFrameMacros
using StatsBase
using DataFrames

const VP = VowelPerception

df = @chain begin
    VP.data_all()
    @transform begin
        :ffc = VP.scale_ff(:ff)
        :f0c = VP.scale_f0(:f0)
        :morphc = VP.scale_morph(:morph)
    end
    # @subset :experiment == "blocked" || :experiment == "precursor"
    # @subset :experiment == ("blocked" ||  "precursor")
    # @subset( (:experiment=="blocked") || (:experiment =="precursor"))
    # filter((:experiment => (:experiment .== "blocked") .| (:experiment .== "precursor")), _)
    # filter((:experiment => ( == "blocked") || ( == "precursor")), _)
    # filtervals = ["blocked", "precursor"]
    # filter(:experiment -> contains(==, filtervals, :experiment), _)
    @aside println(combine(groupby(_, [:experiment]), nrow => :count))

    # filter(experiment -> experiment == "blocked" || experiment == "precursors" , _)
    # filter!(experiment -> experiment in ["blocked", "precursors"], _)
    filter(:experiment => experiment -> experiment == "blocked" || experiment == "precursors", _)
    # filter([:x, :y] => (x, y) -> x == 1 || y == "b", df)

    # filter(:id => >(6), _)
    # filter(!"mixed",_)
    @aside println(combine(groupby(_, [:experiment]), nrow => :count))

    select([:experiment, :participant, :wordpair, :morphc, :ffc, :f0c, :response_o])
    transform(:participant => categorical => :participant)
    transform(:experiment => categorical => :experiment)
    @aside levels!(_.experiment, ["blocked", "precursors"])
end

println(first(df,6))
println(combine(groupby(df, [:experiment]), nrow => :count))

# println(unique(df, :experiment))


model = fit(
    MixedModel,
    @formula(response_o ~ experiment * morphc + experiment * ffc + experiment * f0c + (1 + morphc | participant) + (1 + morphc | wordpair)),
    df,
    Bernoulli(),
)

model2 = fit(
    MixedModel,
    @formula(response_o ~ 1+ experiment + experiment * ffc + experiment * f0c + (1 | participant) + (1 | wordpair)),
    df,
    Bernoulli(),

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

open(VP.reportpath("mixedmodel_shifts_block_prec.txt"), "w") do file
    print(file, model)
end

open(VP.reportpath("mixedmodel_shifts_all_block_prec_latex.txt"), "w") do file
    print(file, latexed_table(model))
end

##############
# Dummy code #
##############

df = @chain begin
    VP.data_all()
    @transform begin
        :ffc = VP.scale_ff(:ff)
        :f0c = VP.scale_f0(:f0)
        :morphc = VP.scale_morph(:morph)
    end
    @transform :dummy_c = 

    @aside println(combine(groupby(_, [:experiment]), nrow => :count))

    filter(:experiment => experiment -> experiment == "blocked" || experiment == "precursors", _)

    @aside println(combine(groupby(_, [:experiment]), nrow => :count))

    select([:experiment, :participant, :wordpair, :morphc, :ffc, :f0c, :response_o])
    transform(:participant => categorical => :participant)
    transform(:experiment => categorical => :experiment)
    @aside levels!(_.experiment, ["blocked", "precursors"])
end

println(first(df,6))
println(combine(groupby(df, [:experiment]), nrow => :count))

# println(unique(df, :experiment))


model = fit(
    MixedModel,
    @formula(response_o ~ experiment * morphc + experiment * ffc + experiment * f0c + (1 + morphc | participant) + (1 + morphc | wordpair)),
    df,
    Bernoulli(),
)

model2 = fit(
    MixedModel,
    @formula(response_o ~ experiment * morphc + experiment * ffc + experiment * f0c + (1 + morphc | participant) + (1 + morphc | wordpair)),
    df,
    Bernoulli(),

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

open(VP.reportpath("mixedmodel_shifts_block_prec.txt"), "w") do file
    print(file, model)
end

open(VP.reportpath("mixedmodel_shifts_all_block_prec_latex.txt"), "w") do file
    print(file, latexed_table(model))
end
