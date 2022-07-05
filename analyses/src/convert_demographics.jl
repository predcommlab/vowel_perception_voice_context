using CSV
using DataFrames
using DataFrameMacros
using Chain
using Dates


path_mixed = raw"C:\Users\Krumbiegel\projects\vowel perception\data\2021-04-22 experiment 1a"
path_blocked = raw"C:\Users\Krumbiegel\projects\vowel perception\data\2021-04-22 experiment 1b"
path_precursors = raw"C:\Users\Krumbiegel\projects\vowel perception\data\2021-04-26 experiment 2"

accepted_ids = mapreduce(vcat, [path_mixed, path_blocked, path_precursors]) do path
    csvpaths = filter(endswith(".csv"), readdir(joinpath(path, "accepted"), join = true))
    map(csvpaths) do csvp
        match(r"ID_(.*?)_", csvp)[1]
    end
end


paths_and_experiments = zip([path_mixed, path_blocked, path_precursors], ["mixed", "blocked", "precursors"])

df = mapreduce((x, y) -> vcat(x, y, cols = :union), enumerate(paths_and_experiments)) do (i, (path, experiment))
    @chain begin
        CSV.read(joinpath(path, "demographics.csv"), DataFrame)
        @subset :status == "APPROVED"
        @transform :completed_date_time = DateTime(:completed_date_time, "y-m-d H:M:S.s")
        @transform :started_datetime = DateTime(:started_datetime, "y-m-d H:M:S.s")
        @subset :participant_id ∈ accepted_ids
        sort!(:completed_date_time)
        @transform :participant = @c eachindex(:participant_id) .+ (i - 1) * 30
        @transform :experiment = experiment
        select(["experiment", "participant", "age", "Sex", "Nationality"])
        rename!(_, names(_) .=> lowercase.(names(_)))
        @transform :sex = lowercase(:sex)[1]
        @transform :nationality = lowercase(:nationality)
        disallowmissing!
    end
end

CSV.write("data/demographics.csv", df)