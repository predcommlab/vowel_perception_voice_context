using SimpleWeave


simpleweave(
    joinpath(@__DIR__, "report_demographics.jl"),
    joinpath(@__DIR__, "..", "reports")
)