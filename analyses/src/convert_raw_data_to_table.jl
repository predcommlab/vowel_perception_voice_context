using CSV
using DataFrames
using DataFrameMacros
using Chain
using StatsBase


# set up raw data directories
dir_mixed_raw = raw"C:\Users\Krumbiegel\projects\vowel perception\data\2021-04-22 experiment 1a\accepted"
dir_blocked_raw = raw"C:\Users\Krumbiegel\projects\vowel perception\data\2021-04-22 experiment 1b\accepted"
dir_precursors_raw = raw"C:\Users\Krumbiegel\projects\vowel perception\data\2021-04-26 experiment 2\accepted"


function timestamp(filename)
    s = match(r"SESSION_(.*)\.csv", filename)[1]
    String(s)
end

function read_csvs(folder)
    files = filter(endswith(".csv"), readdir(folder, join = true))
    mapreduce(vcat, files) do file
        df = CSV.read(file, DataFrame)
        insertcols!(df, :timestamp => timestamp(basename(file)))
    end
end

# mixed
@chain begin
    read_csvs(dir_mixed_raw)
    insertcols!(:experiment => "mixed")
    transform(_, :timestamp => function(ts)
        uts = sort(unique(ts))
        d = Dict(uts .=> 1:length(uts))
        [d[t] for t in ts]
    end => :participant)
    @subset(isequal(:task_segment, "word-classification"))
    select([:experiment, :participant, :trial_number, :f0, :ff, :morph, :response, :rt,
        :speaker, :time_elapsed, :wordpair])
    transform(:response => ByRow(x -> x ∉ ("o", "u") ? missing : x), renamecols = false)
    transform(:rt => ByRow(x -> all(isdigit, x) ? parse(Int, x) : missing), renamecols = false)
    disallowmissing(Not([:response, :rt]))
    transform(
        :response => (x -> x .== "o") => :response_o,
        :rt => (x -> x .> 400) => :rt_valid,
        :time_elapsed => (x -> x ./ 1000) => :time_elapsed,
    )
    @transform(:trial_valid = !ismissing(:response) && :rt_valid)
    sort([:participant, :trial_number])
    @aside CSV.write(joinpath(@__DIR__, "..", "data", "mixed_all.csv"), _)
    subset(:trial_valid)
    select(Not([:rt_valid, :trial_valid]))
    disallowmissing
    @aside CSV.write(joinpath(@__DIR__, "..", "data", "mixed_valid.csv"), _)
end

## blocked
function add_blockspeaker!(df)
    @chain df begin
        groupby([:PROLIFIC_PID, :block])
        transform!(:speaker => function(bs)
            c = countmap(bs)
            ma = findmax(last, c)
            ma[2]
        end => :context_speaker)
    end
end

@chain begin
    read_csvs(dir_blocked_raw)
    insertcols!(:experiment => "blocked")
    filter(:task_segment => isequal("word-classification"), _)
    add_blockspeaker!
    transform(_, :timestamp => function(ts)
        uts = sort(unique(ts))
        d = Dict(uts .=> 1:length(uts))
        [d[t] for t in ts] .+ 30 # offset for unique participant numbers
    end => :participant)
    select([:experiment, :participant, :trial_number, :f0, :ff, :morph, :response, :rt,
        :speaker, :context_speaker, :time_elapsed, :wordpair])
    @transform(:response = :response ∉ ("o", "u") ? missing : :response)
    @transform(:rt = all(isdigit, :rt) ? parse(Int, :rt) : missing)
    disallowmissing!(Not([:response, :rt]))
    @transform(
        :response_o = :response == "o",
        :rt_valid = :rt > 400,
        :time_elapsed = :time_elapsed / 1000,
    )
    @transform(:trial_valid = !ismissing(:response) && :rt_valid)
    sort([:participant, :trial_number])
    @aside CSV.write(joinpath(@__DIR__, "..", "data", "blocked_all.csv"), _)
    subset(:trial_valid)
    select(Not([:rt_valid, :trial_valid]))
    disallowmissing
    @aside CSV.write(joinpath(@__DIR__, "..", "data", "blocked_valid.csv"), _)
    _
end


## precursors
@chain begin
    read_csvs(dir_precursors_raw)
    insertcols!(:experiment => "precursors")
    transform(_, :timestamp => function(ts)
        uts = sort(unique(ts))
        d = Dict(uts .=> 1:length(uts))
        [d[t] for t in ts] .+ 60 # offset for unique participant numbers
    end => :participant)
    @subset(isequal(:task_segment, "word-classification") && !ismissing(:precursor_speaker))
    select([:experiment, :participant, :trial_number, :f0, :ff, :morph, :response, :rt,
        :speaker, :precursor_speaker, :time_elapsed, :wordpair, :precursor])
    rename!(:precursor_speaker => :context_speaker) # same name as in blocked
    @transform(:response = :response ∉ ("o", "u") ? missing : :response)
    @transform(:rt = all(isdigit, :rt) ? parse(Int, :rt) : missing)
    disallowmissing(Not([:response, :rt]))
    @transform(
        :response_o = :response == "o",
        :rt_valid = :rt > 400,
        :time_elapsed = :time_elapsed / 1000,
    )
    @transform(:trial_valid = !ismissing(:response) && :rt_valid)
    sort([:participant, :trial_number])
    @aside CSV.write(joinpath(@__DIR__, "..", "data", "precursors_all.csv"), _)
    filter(:trial_valid => identity, _)
    select(Not([:rt_valid, :trial_valid]))
    disallowmissing
    @aside CSV.write(joinpath(@__DIR__, "..", "data", "precursors_valid.csv"), _)
    _
end

## combined
@chain begin
    vcat(
        CSV.read(joinpath(@__DIR__, "..", "data", "mixed_valid.csv"), DataFrame),
        CSV.read(joinpath(@__DIR__, "..", "data", "blocked_valid.csv"), DataFrame),
        CSV.read(joinpath(@__DIR__, "..", "data", "precursors_valid.csv"), DataFrame),
        cols = :union)
    CSV.write(joinpath(@__DIR__, "..", "data", "all_valid.csv"), _)
end 