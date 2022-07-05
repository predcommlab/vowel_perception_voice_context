using Revise
includet(joinpath(@__DIR__, "VowelPerception.jl"))

const VP = VowelPerception

using CSV
using CairoMakie
using DataFrameMacros
using Chain
using DataFrames
using StatsBase
# using StatsPlots
using Turing
using CategoricalArrays
using MCMCChains
using ArviZ

function compute_hillenbrand_slope_variance()
    df = CSV.read(VP.datadir("hillenbrand.csv"), DataFrame)

    barkmeans = @chain df begin
        select(Not(r"f\d_\d\d"))
        select(Not([:filename, :duration_msec]))
        transform([:f0_steady, :f1_steady, :f2_steady, :f3_steady]
            .=> ByRow(VP.hz_to_bark)
            .=> [:f0b, :f1b, :f2b, :f3b]
        )
        # we just call the mean of f1, f2 and f3 "ffb" as a proxy
        @transform :ffb = mean((:f1b, :f2b, :f3b))
        dropmissing(:ffb)
        groupby([:type, :number])
        combine([:f0b, :ffb] .=> mean)
    end

    @model function ff_prediction(f0s, ffs)

        intercept ~ Uniform(0, 20)
        slope ~ Normal(0, 5)
        σ² ~ truncated(Normal(0, 100), 0, Inf)
    
        for i in eachindex(f0s)
            ffs[i] ~ Normal(intercept + slope * f0s[i], sqrt(σ²))
        end
    end
    
    iterations = 4000
    n_chains = 4
    
    chain = sample(
        ff_prediction(barkmeans.f0b_mean, barkmeans.ffb_mean),
        NUTS(0.65),
        MCMCThreads(),
        iterations,
        n_chains)

    # describe(chain, showall=true)
    # summarize(chain)
    x = summarize(chain; sections=[:parameters, :internals])
    # println(chain[:lp])
    open("reports/hillenbrand_model_summarize.txt", "w") do file
        show(file, MIME"text/plain"(), x)
    end
    
    intercept = mean(chain[:intercept])
    slope = mean(chain[:slope])
    variance = mean(chain[:σ²])
    (; intercept, slope, variance)
end

hb = compute_hillenbrand_slope_variance()
const hb_slope = hb.slope
const hb_intercept = hb.intercept
const hb_variance = hb.variance

function predict_ff_from_f0(hz)
    intercept = hb_intercept
    slope = hb_slope
    VP.bark_to_hz(VP.hz_to_bark(hz) * slope + intercept)
end
# other way around? iterative?

guessed_ff = predict_ff_from_f0(170)

# function: combination of distributions for later (prior with context distribution)
# math?
function combine_distributions(n1::Normal, n2::Normal)
    t = 1 / n1.σ ^ 2
    u = 1 / n2.σ ^ 2
    m = n1.μ
    n = n2.μ
    Normal( (t * m + u * n) / (t + u) , sqrt(1 / (t + u)))
end

# data
n_o_df = @chain VP.data_all() begin
    @subset :experiment != "precursors"
    @transform(:f0_bark = VP.hz_to_bark(170 .* :f0))
    @transform(:ff_bark = VP.hz_to_bark(guessed_ff * :ff))
    @transform(:exp_index = :experiment == "mixed" ? 1 : 2)
    @transform(AsTable = begin
        cff, cf0 = if :exp_index == 1
            guessed_ff, 170.0
        else
            if :context_speaker == "high_ff"
                1.146 * guessed_ff, 170.0
            elseif :context_speaker == "low_ff"
                0.873 * guessed_ff, 170.0
            elseif :context_speaker == "high_f0"
                guessed_ff, 170 * 1.24
            elseif :context_speaker == "low_f0"
                guessed_ff, 170 * 0.806
            end
        end
        (context_ff = VP.hz_to_bark(cff), context_f0 = VP.hz_to_bark(cf0))
    end)
    select(Not([:wordpair, :time_elapsed, :experiment, :rt, :trial_number,]))
    @transform(:morph_scaled = :morph * 2 - 1)
    # now convert to number of o responses
    groupby([:exp_index, :speaker, :context_speaker, :ff_bark, :f0_bark, :context_ff, :context_f0, :morph_scaled], sort = true)
    @combine(:n_o = sum(:response_o), :n = length(:response_o))
end

open("reports/data.txt", "w") do file
    show(file, MIME"text/plain"(), n_o_df)
end



@model function pred_error_binomial(f0s, ffs, context_ffs, context_f0s, morphs,
    n_os, ns, experiment_indices, prior_variance)

    intercept = hb_intercept
    slope = hb_slope

    # truncated because o² can only be positive
    σ²_ff ~ truncated(Normal(0, 100), 0, Inf)
    pred_err_scale_1a ~ Normal(0, 2)
    # why two scales ? --> one for mixed one for blocked 
    pred_err_scale_1b ~ Normal(0, 2)
    log_morphscale ~ Normal(0, 3)

    for i in eachindex(f0s)
        stimulus_ff = ffs[i]
        # the distribution for the guess of the current formant frequency mean
        # given the stimulus pitch
        ei = experiment_indices[i]

        # this is the mean formant value from the preceding trials for blocked, or of the center
        # speaker for mixed
        # in bark
        context_ff = context_ffs[i]
        # this is the pitch value from the preceding trials for blocked, or of the center
        # speaker for mixed
        context_f0 = context_f0s[i]
        # what is the expectation of ff given the context f0?
        context_ff_from_f0_prediction = intercept + slope * context_f0
        # this is the learned error between what the f0 values of the previous trials predict and
        # what the ff values of the previous trials were
        # which is actually always zero for the mixed condition if we always use the center speaker
        # as the context, because the prediction for the ff of the center speaker is exactly what its
        # f0 value would predict (because the prediction itself is derived from this value)
        # no flexibility in f0 here 
        learned_ff_prediction_error = context_ff - context_ff_from_f0_prediction

        # this value controls how strongly the learned_ff_prediction_error is neutralized
        pred_err_scale = ei == 1 ? pred_err_scale_1a : pred_err_scale_1b

        # the formant value that would be expected given the stimulus f0
        stimulus_ff_from_f0_prediction = intercept + slope * f0s[i]
        # this is combined with the scaled learned_ff_prediction_error because the subjectsusing 
        # should have learned (at least in the blocked condition)
        # to counter the prediction error between
        # stimulus_ff_from_f0_prediction and stimulus_ff
        context_biased_ff_prediction = stimulus_ff_from_f0_prediction + pred_err_scale * learned_ff_prediction_error

        # the full prior's mean is stimulus_ff_from_f0_prediction biased by the
        # learned_ff_prediction_error, the uncertainty of that estimatio is given by
        # prior_variance
        # local prior
        prior = Normal(context_biased_ff_prediction, sqrt(prior_variance))

        # the distribution representing the mean formant frequency
        # estimated from the stimulus
        stimulus_ff_distribution = Normal(stimulus_ff, sqrt(σ²_ff))

        # the final estimation of formant mean is made by combining
        # prior and stimulus distributions
        combined = combine_distributions(prior, stimulus_ff_distribution)

        # the probability of responding o depends on the location where the morph
        # formant pattern lies on the pdf or cdf of the combined distribution
        # the impact of the morphs is scaled by morphscale, which works in combination
        # with the variance of the combined distribution
        p_o = cdf(combined, exp(log_morphscale) * morphs[i] + stimulus_ff)



        # why stimulus_ff and not predicted stimulus_ff? -> cdf verschoben ist schon template prediction
        

        n_os[i] ~ Binomial(ns[i], p_o)
    end
end

@model function pred_error_binomial_null(f0s, ffs, context_ffs, context_f0s, morphs,
    n_os, ns, experiment_indices, prior_variance)

    intercept = hb_intercept
    slope = hb_slope

    # truncated because o² can only be positive
    σ²_ff ~ truncated(Normal(0, 100), 0, Inf)
    # pred_err_scale_1a ~ Normal(0, 2)
    # pred_err_scale_1b ~ Normal(0, 2)
    log_morphscale ~ Normal(0, 3)

    for i in eachindex(f0s)
        stimulus_ff = ffs[i]
        # the distribution for the guess of the current formant frequency mean
        # given the stimulus pitch
        ei = experiment_indices[i]

       
        # this is the learned error between what the f0 values of the previous trials predict and
        # what the ff values of the previous trials were
        # which is actually always zero for the mixed condition if we always use the center speaker
        # as the context, because the prediction for the ff of the center speaker is exactly what its
        # f0 value would predict (because the prediction itself is derived from this value)
        # no flexibility in f0 here 
        learned_ff_prediction_error = 0

        # this value controls how strongly the learned_ff_prediction_error is neutralized
        # pred_err_scale = ei == 1 ? pred_err_scale_1a : pred_err_scale_1b

        # the formant value that would be expected given the stimulus f0
        stimulus_ff_from_f0_prediction = intercept + slope * f0s[i]
        # this is combined with the scaled learned_ff_prediction_error because the subjectsusing 
        # should have learned (at least in the blocked condition)
        # to counter the prediction error between
        # stimulus_ff_from_f0_prediction and stimulus_ff
        context_biased_ff_prediction = stimulus_ff_from_f0_prediction + learned_ff_prediction_error

        # the full prior's mean is stimulus_ff_from_f0_prediction biased by the
        # learned_ff_prediction_error, the uncertainty of that estimatio is given by
        # prior_variance
        # local prior
        prior = Normal(context_biased_ff_prediction, sqrt(prior_variance))

        # the distribution representing the mean formant frequency
        # estimated from the stimulus
        stimulus_ff_distribution = Normal(stimulus_ff, sqrt(σ²_ff))

        # the final estimation of formant mean is made by combining
        # prior and stimulus distributions
        combined = combine_distributions(prior, stimulus_ff_distribution)

        # the probability of responding o depends on the location where the morph
        # formant pattern lies on the pdf or cdf of the combined distribution
        # the impact of the morphs is scaled by morphscale, which works in combination
        # with the variance of the combined distribution
        p_o = cdf(combined, exp(log_morphscale) * morphs[i] + stimulus_ff)

        # why stimulus_ff and not predicted stimulus_ff? -> cdf verschoben ist schon template prediction
        n_os[i] ~ Binomial(ns[i], p_o)
    end
end

param_mod = pred_error_binomial(n_o_df.f0_bark, n_o_df.ff_bark, n_o_df.context_ff, n_o_df.context_f0,n_o_df.morph_scaled, n_o_df.n_o, n_o_df.n, n_o_df.exp_index,hb_variance )
param_mod_null = pred_error_binomial_null(n_o_df.f0_bark, n_o_df.ff_bark, n_o_df.context_ff, n_o_df.context_f0,n_o_df.morph_scaled, n_o_df.n_o, n_o_df.n, n_o_df.exp_index,hb_variance )


chains = let
    iterations = 8000
    n_chains = 4

    sample(
        param_mod,
        NUTS(0.65),
        MCMCThreads(),
        iterations,
        n_chains,
        progress = false,
    )
end


# model evidence - calculate loglikelihood 
ℓ = Turing.pointwise_loglikelihoods(param_mod, chains)
ℓ_mat = reduce(hcat, values(ℓ));
ℓ_arr = reshape(ℓ_mat, 1, size(ℓ_mat)...); # (chain_idx, sample_idx, parameter_idx)
# println(size(ℓ_arr))

data = ArviZ.from_mcmcchains(
    chains,
    library = "Turing",
    log_likelihood = Dict("y" => ℓ_arr)
)

loo_result = ArviZ.arviz.loo(data)

# x = summarize(chains; sections=[:parameters, :internals])
# println(chain[:lp])
# open("reports/bayesian_model_summarize_blocked.csv", "w") do file
#     show(file, MIME"text/csv"(), x)
# end



open("reports/bayesian_model_output.txt", "w") do file
    show(file, MIME"text/plain"(), chains)
end

open("reports/bayesian_model_loo.txt", "w") do file
    show(file, MIME"text/plain"(), loo_result)
end

# for null model
chains_mod_null = let
    iterations = 8000
    n_chains = 4

    sample(
        param_mod_null,
        NUTS(0.65),
        MCMCThreads(),
        iterations,
        n_chains,
        progress = false,
    )
end

ℓ = Turing.pointwise_loglikelihoods(param_mod_null, chains_mod_null)
ℓ_mat = reduce(hcat, values(ℓ));
ℓ_arr = reshape(ℓ_mat, 1, size(ℓ_mat)...); # (chain_idx, sample_idx, parameter_idx)

data_null = ArviZ.from_mcmcchains(
    chains_mod_null,
    library = "Turing",
    log_likelihood = Dict("y" => ℓ_arr)
)

loo_result_null = ArviZ.arviz.loo(data_null)

# x = summarize(chains; sections=[:parameters, :internals])
# # println(chain[:lp])
# open("reports/bayesian_model_3_summarize_blocked.csv", "w") do file
#     show(file, MIME"text/csv"(), x)
# end

open("reports/bayesian_model_null_output.txt", "w") do file
    show(file, MIME"text/plain"(), chains_mod_null)
end

open("reports/bayesian_model_null_loo.txt", "w") do file
    show(file, MIME"text/plain"(), loo_result_null)
end


### compare models
compare_dict = Dict(
           "model" => loo_result,
           "null model" => loo_result_null,
       );

results = compare(compare_dict; ic="loo")

open("reports/bayesian_model_comparison_latex.txt", "w") do file
    show(file, MIME"text/plain"(), latexed_table(results))
end

# describe(chains, showall=true) # Lists statistics of the samples.
# plot(chains, colordim = :parameters)


### prior and posterior predictive check
prior = Turing.sample(param_mod, Prior(), 8000);

# vector of missing number of o-responses to predict
n_os_missing = Vector{Missing}(missing, nrow(n_o_df)) 

# feed into model
model_missing = pred_error_binomial(
    n_o_df.f0_bark, # f0s 
    n_o_df.ff_bark, # ffs
    n_o_df.context_ff, # context ff
    n_o_df.context_f0, # context f0
    n_o_df.morph_scaled, # :morph_scaled = :morph * 2 - 1
    n_os_missing,
    n_o_df.n, # ns
    n_o_df.exp_index,
    hb_variance, # prior_variance
    )

prior_predictive = Turing.predict(model_missing, prior)
posterior_predictive = Turing.predict(model_missing, chains)

open("reports/prior_predictive_raw.txt", "w") do file
    show(file, MIME"text/plain"(), prior_predictive)
end

open("reports/posterior_predictive_raw.txt", "w") do file
    show(file, MIME"text/plain"(), posterior_predictive)
end

# plot prior and posterior predictive check
# @chain begin
#     StatsPlots.plot(prior_predictive; colordim = :parameter)
#     save("figures/prior_pred_para.png", _)
# end

# @chain begin
#     StatsPlots.plot(posterior_predictive; colordim = :parameter)
#     save("figures/post_pred_para.png", _)
# end

# @chain begin
#     StatsPlots.plot(prior_predictive; colordim = :chain)
#     save("figures/prior_pred_chain.png", _)
# end

# @chain begin
#     StatsPlots.plot(posterior_predictive; colordim = :chain)
#     save("figures/post_pred_chain.png", _)
# end

# @show summarystats(prior_predictive[:, 1:5, :])
# @show summarystats(posterior_predictive[:, 1:5, :])


function pred_error_binomial_prediction(f0s, ffs, context_ffs, context_f0s, morphs,
    experiment_indices, prior_variance, σ²_ff, pred_err_scale_1a,
    pred_err_scale_1b, log_morphscale)
    # pred_err_scale_1a

    intercept = hb_intercept
    slope = hb_slope

    map(eachindex(f0s)) do i
        # the distribution for the guess of the current formant frequency mean
        # given the stimulus pitch
        ei = experiment_indices[i]
        stimulus_ff = ffs[i]

        context_ff = context_ffs[i]
        context_f0 = context_f0s[i]
        context_ff_from_f0_prediction = intercept + slope * context_f0
        # actually for 1a this is always 0 currently, so there's no use for pred_err_scale_1a
        learned_ff_prediction_error = context_ff - context_ff_from_f0_prediction

        pred_err_scale = ei == 1 ? pred_err_scale_1a : pred_err_scale_1b

        stimulus_ff_from_f0_prediction = intercept + slope * f0s[i]
        context_biased_ff_prediction = stimulus_ff_from_f0_prediction + pred_err_scale * learned_ff_prediction_error

        prior = Normal(context_biased_ff_prediction, sqrt(prior_variance))

        # the distribution representing the guess of formant frequency mean from
        # the raw input itself
        ff_input = Normal(stimulus_ff, sqrt(σ²_ff))

        # combine these two into one normal distribution
        combined = combine_distributions(prior, ff_input)

        p_o = cdf(combined, exp(log_morphscale) * morphs[i] + stimulus_ff)
        (; pred = p_o, pred_err_scale, learned_ff_prediction_error, context_biased_ff_prediction, mu_combined = combined.μ)
    end |> DataFrame
end

function pred_error_binomial_prediction_null(f0s, ffs, context_ffs, context_f0s, morphs,
    experiment_indices, prior_variance, σ²_ff,  log_morphscale)
    # pred_err_scale_1a,
    # pred_err_scale_1b,

    intercept = hb_intercept
    slope = hb_slope

    map(eachindex(f0s)) do i
        # the distribution for the guess of the current formant frequency mean
        # given the stimulus pitch
        ei = experiment_indices[i]
        stimulus_ff = ffs[i]

        # actually for 1a this is always 0 currently, so there's no use for pred_err_scale_1a
        learned_ff_prediction_error = 0

        # pred_err_scale = ei == 1 ? pred_err_scale_1a : pred_err_scale_1b

        stimulus_ff_from_f0_prediction = intercept + slope * f0s[i]
        context_biased_ff_prediction = stimulus_ff_from_f0_prediction + learned_ff_prediction_error

        prior = Normal(context_biased_ff_prediction, sqrt(prior_variance))

        # the distribution representing the guess of formant frequency mean from
        # the raw input itself
        ff_input = Normal(stimulus_ff, sqrt(σ²_ff))

        # combine these two into one normal distribution
        combined = combine_distributions(prior, ff_input)

        p_o = cdf(combined, exp(log_morphscale) * morphs[i] + stimulus_ff)
        # pred_err_scale
        (; pred = p_o,  learned_ff_prediction_error, context_biased_ff_prediction, mu_combined = combined.μ)
    end |> DataFrame
end

predictions = pred_error_binomial_prediction(
    n_o_df.f0_bark,
    n_o_df.ff_bark,
    n_o_df.context_ff,
    n_o_df.context_f0,
    n_o_df.morph_scaled,
    n_o_df.exp_index,
    hb_variance,
    mean(chains["σ²_ff"]),
    mean(chains["pred_err_scale_1a"]),
    mean(chains["pred_err_scale_1b"]),
    mean(chains["log_morphscale"]),
)

predictions_null = pred_error_binomial_prediction_null(
    n_o_df.f0_bark,
    n_o_df.ff_bark,
    n_o_df.context_ff,
    n_o_df.context_f0,
    n_o_df.morph_scaled,
    n_o_df.exp_index,
    hb_variance,
    mean(chains_mod_null["σ²_ff"]),
    # mean(chains_mod_null["pred_err_scale_1a"]),
    # mean(chains_mod_null["pred_err_scale_1b"]),
    mean(chains_mod_null["log_morphscale"]),
)

errors = @chain n_o_df begin
    hcat(predictions)
    @transform :p_o = :n_o / :n  # from nr of total o responses to probability
    @transform(:err = abs(:p_o - :pred))
end

# open("reports/bayesian_model_errors.csv", "w") do file
#     show(file, MIME"text/csv"(), errors)
# end

errors_null = @chain n_o_df begin
    hcat(predictions_null)
    @transform :p_o = :n_o / :n
    @transform(:err = abs(:p_o - :pred))
end

# open("reports/bayesian_model_errors_null.csv", "w") do file
#     show(file, MIME"text/csv"(), errors_null)
# end

function visualize_model(errors)

    errors = @chain errors begin
        @transform(:context_speaker = @c categorical(:context_speaker))
        @transform(:speaker = @c categorical(:speaker))
    end
    levels!(errors.context_speaker, ["low_ff", "high_ff", "low_f0", "high_f0", "center"])
    levels!(errors.speaker, ["low_ff", "high_ff", "low_f0", "high_f0", "center"])
    sort!(errors, [:context_speaker, :speaker])

    # @show errors

    f = with_theme(VP.base_theme()) do
        Figure(resolution = (800, 500))
    end

    axs = Matrix{Axis}(undef, 2, 5)

    for (row, gdf) in enumerate(groupby(errors, :exp_index))
        for (col, ggdf) in enumerate(groupby(gdf, :speaker))
            ax = axs[row, col] = Axis(f[row, col], titlefont = VP.normal_font())
            ax.xticks = (-1:1, ["u", "¹/₂", "o"])
            ax.yticks = (0:0.5:1, ["u", "¹/₂", "o"])
            xlims!(ax, -1, 1)
            ylims!(ax, 0, 1)

            if row == 1
                ax.title = VP.speaker_pretty(ggdf.speaker[1])
            end

            for gggdf in groupby(ggdf, :context_speaker)
                c = VP.contextcolor(gggdf.context_speaker[1])
                # measured
                lines!(gggdf.morph_scaled, gggdf.p_o, color = c)
                # predicted
                lines!(gggdf.morph_scaled, gggdf.pred, linestyle = :dash, color = c)
            end
        end
        Label(f[row, 5, Right()],
            gdf.exp_index[1] == 1 ? "Mixed" : "Blocked",
            rotation = pi/2,
            font = VP.medium_font()
        )
    end

    hideydecorations!.(axs[:, 2:end])
    hidexdecorations!.(axs[1, :])

    Label(f[1:2, 1, Left()], "Response p(o)", rotation = pi/2,
        padding = (0, 30, 0, 0),
        font = VP.medium_font()
    )

    Label(f[2, 3, Bottom()], "Morph level",
        padding = (0, 0, 0, 30),
        font = VP.medium_font()
    )

    Label(f[1, 3, Top()], "Speaker",
        padding = (0, 0, 30, 0),
        font = VP.medium_font()
    )

    colororder = unique(errors.context_speaker)

    Legend(f[3, :],
        [
            [LineElement(color = c) for c in VP.contextcolor.(colororder)],
            [LineElement(linestyle = ls) for ls in [:solid, :dash]],
        ],
        [
            VP.speaker_pretty.(colororder),
            ["Measured", "Predicted"]
        ],
        ["Context", "Type"],
        orientation = :horizontal,
        groupgap = 50,
        titleposition = :left,
        nbanks = 2,
        titlegap = 15,
    )

    f
end



@chain begin
    visualize_model(errors)
    save("figures/bayesian_model.pdf", _)
end

@chain begin
    visualize_model(errors_null)
    save("figures/bayesian_model_null.pdf", _)
end




