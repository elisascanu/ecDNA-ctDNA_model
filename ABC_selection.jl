include("model.jl")

using Random
using Statistics
using Distributions

const R_PRIOR = Uniform(0.02, 0.06)
const LOGN_PRIOR = Uniform(7.0, 9.0)
const S_POSITIVE_PRIOR = Uniform(1.05, 3.0)

function abc_distance(sim_ecDNA, sim_chr, obs_ecDNA, obs_chr)
    length(sim_ecDNA) == length(obs_ecDNA) || error("ecDNA series have different lengths.")
    length(sim_chr) == length(obs_chr) || error("Chromosomal series have different lengths.")

    d_ecDNA = mean(abs.(log1p.(sim_ecDNA) .- log1p.(obs_ecDNA)))
    d_chr = mean(abs.(log1p.(sim_chr) .- log1p.(obs_chr)))
    return sqrt((d_ecDNA^2 + d_chr^2) / 2)
end

function build_selection_reference(
    n_particles;
    sample_times = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0],
    r_prior = R_PRIOR,
    logN_prior = LOGN_PRIOR,
    s_positive_prior = S_POSITIVE_PRIOR
)
    reference = NamedTuple[]

    while length(reference) < n_particles
        r = rand(r_prior)
        logN = rand(logN_prior)
        N_first = 10.0^logN

        positive = rand(Bool)
        s_max = positive ? rand(s_positive_prior) : 1.0

        sim = simulate_longitudinal_blood_from_Nfirst(
            r, s_max, N_first;
            sample_times = sample_times
        )

        sim === nothing && continue

        push!(reference, (
            r = r,
            logN = logN,
            N_first = N_first,
            s_max = s_max,
            positive = positive,
            ecDNA = sim.ecDNA,
            chromosomal = sim.chromosomal
        ))
    end

    return reference
end

function selection_abc(
    obs_ecDNA,
    obs_chromosomal,
    reference;
    n_accept = 50
)
    n_accept <= length(reference) || error("n_accept exceeds reference-table size.")

    distances = [
        abc_distance(
            p.ecDNA, p.chromosomal,
            obs_ecDNA, obs_chromosomal
        )
        for p in reference
    ]

    accepted_idx = sortperm(distances)[1:n_accept]
    accepted = reference[accepted_idx]

    positive_support = mean(p.positive for p in accepted)
    logN_values = [p.logN for p in accepted]
    r_values = [p.r for p in accepted]

    return (
        positive_support = positive_support,
        classified_positive = positive_support > 0.5,
        median_logN = median(logN_values),
        logN_q25 = quantile(logN_values, 0.25),
        logN_q75 = quantile(logN_values, 0.75),
        median_r = median(r_values),
        accepted = accepted,
        distances = distances[accepted_idx]
    )
end

function simulate_selection_patient(
    r, s_max, N_first;
    sample_times = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0]
)
    return simulate_longitudinal_blood_from_Nfirst(
        r, s_max, N_first;
        sample_times = sample_times
    )
end

function validate_selection(
    reference;
    s_values = [1.0, 1.05, 1.10, 1.25, 1.50, 1.75, 2.0, 2.5],
    n_patients = 20,
    n_accept = 50,
    r_prior = R_PRIOR,
    logN_prior = LOGN_PRIOR,
    sample_times = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0]
)
    results = NamedTuple[]

    for s_true in s_values
        for patient in 1:n_patients
            obs = nothing
            r_true = NaN
            logN_true = NaN

            while obs === nothing
                r_true = rand(r_prior)
                logN_true = rand(logN_prior)
                obs = simulate_selection_patient(
                    r_true, s_true, 10.0^logN_true;
                    sample_times = sample_times
                )
            end

            fit = selection_abc(
                obs.ecDNA,
                obs.chromosomal,
                reference;
                n_accept = n_accept
            )

            push!(results, (
                patient = patient,
                s_true = s_true,
                r_true = r_true,
                logN_true = logN_true,
                positive_support = fit.positive_support,
                classified_positive = fit.classified_positive,
                median_logN = fit.median_logN,
                median_r = fit.median_r
            ))
        end
    end

    return results
end

function summarise_selection_validation(results)
    s_values = sort(unique(r.s_true for r in results))
    summary = NamedTuple[]

    for s in s_values
        rows = [r for r in results if r.s_true == s]
        weights = [r.positive_support for r in rows]
        classified = [r.classified_positive for r in rows]

        push!(summary, (
            s_true = s,
            median_positive_support = median(weights),
            q25_positive_support = quantile(weights, 0.25),
            q75_positive_support = quantile(weights, 0.75),
            fraction_classified_positive = mean(classified)
        ))
    end

    return summary
end
