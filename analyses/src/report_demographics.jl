using CSV
using DataFrames
using Chain
using DataFrameMacros
using StatsBase

df = CSV.read(joinpath(@__DIR__, "..", "data", "demographics.csv"), DataFrame)
   
combine(groupby(df, :experiment), :age .=> [minimum, maximum, StatsBase.median])
##
combine(groupby(df, :experiment), :sex => countmap)
