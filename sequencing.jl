include("ABC_treatment.jl")

using Random
using Statistics
using Distributions

const DEPTHS = [Inf, 100.0, 30.0, 10.0, 5.0]
const THRESHOLDS = [1, 2, 3]
const TRUE_D_VALUES = [0.001, 0.005, 0.010]

# Effective molecular-depth model applied to one 15 ml blood sample
function sequencing_observation(ecDNA::Integer, chromosomal::Integer,
                                N_target::Real, L::Int, rng)
    total = ecDNA + chromosomal
    total == 0 && return (ecDNA = 0, chromosomal = 0)

    if isinf(N_target)
        Y_e = Int(ecDNA)
        Y_c = Int(chromosomal)
    else
        p_seq = min(1.0, N_target / total)
        N_seq = rand(rng, Binomial(total, p_seq))

        if N_seq == 0
            return (ecDNA = 0, chromosomal = 0)
        end

        p_ec = ecDNA / total
        p_chr = chromosomal / total

        sampled = rand(
            rng,
            Multinomial(N_seq, [p_ec, p_chr])
        )

        Y_e = sampled[1]
        Y_c = sampled[2]
    end

    return (
        ecDNA = Y_e >= L ? Y_e : 0,
        chromosomal = Y_c
    )
end

function observe_trajectory(ecDNA, chromosomal, N_target, L, rng)
    obs_ec = Int[]
    obs_chr = Int[]

    for i in eachindex(ecDNA)
        x = sequencing_observation(ecDNA[i], chromosomal[i], N_target, L, rng)
        push!(obs_ec, x.ecDNA)
        push!(obs_chr, x.chromosomal)
    end

    return (ecDNA = obs_ec, chromosomal = obs_chr)
end

function observe_reference(reference, N_target, L, rng)
    return [
        begin
            obs = observe_trajectory(p.ecDNA, p.chromosomal, N_target, L, rng)
            (
                r = p.r,
                logN = p.logN,
                N_first = p.N_first,
                D = p.D,
                ecDNA = obs.ecDNA,
                chromosomal = obs.chromosomal
            )
        end
        for p in reference
    ]
end

# Paired validation cohort: same patient background tested at three D values
function build_sequencing_validation(
    n_patients = 50;
    D_values = TRUE_D_VALUES,
    r_prior = R_PRIOR,
    logN_prior = LOGN_PRIOR,
    seed = 271828
)
    rng = MersenneTwister(seed)
    validation = NamedTuple[]
    patient_id = 0
    attempts = 0

    while patient_id < n_patients
        attempts += 1
        r_true = rand(rng, r_prior)
        logN_true = rand(rng, logN_prior)
        patient_seed = seed + attempts
        sims = NamedTuple[]
        evaluable = true

        for D_true in D_values
            Random.seed!(patient_seed)
            sim = simulate_evaluable_treatment_patient(r_true, logN_true, D_true)

            if sim === nothing || !sim.treatment_started || sim.detection_time != 0.0
                evaluable = false
                break
            end

            push!(sims, (
                D_true = D_true,
                ecDNA = copy(sim.sampled_ecDNA),
                chromosomal = copy(sim.sampled_chromosomal)
            ))
        end

        evaluable || continue
        patient_id += 1

        for sim in sims
            push!(validation, (
                patient_id = patient_id,
                r_true = r_true,
                logN_true = logN_true,
                true_D = sim.D_true,
                observed_ecDNA = sim.ecDNA,
                observed_chromosomal = sim.chromosomal
            ))
        end
    end

    return validation, attempts
end

function observe_validation(validation, N_target, L, rng)
    return [
        begin
            obs = observe_trajectory(
                row.observed_ecDNA,
                row.observed_chromosomal,
                N_target,
                L,
                rng
            )
            (
                patient_id = row.patient_id,
                true_D = row.true_D,
                ecDNA = obs.ecDNA,
                chromosomal = obs.chromosomal
            )
        end
        for row in validation
    ]
end

function analyse_sequencing_condition(
    reference,
    validation,
    N_target,
    L;
    replicate = 1,
    acceptance_fraction = 0.01,
    seed = 1_100_000
)
    depth_index = findfirst(==(N_target), DEPTHS)
    threshold_index = findfirst(==(L), THRESHOLDS)

    rng_ref = MersenneTwister(
        seed + 100_000 * depth_index + 10_000 * threshold_index + replicate
    )
    rng_val = MersenneTwister(
        seed + 600_000 + 100_000 * depth_index + 10_000 * threshold_index + replicate
    )

    observed_reference = observe_reference(reference, N_target, L, rng_ref)
    observed_validation = observe_validation(validation, N_target, L, rng_val)

    patient_results = NamedTuple[]

    for row in observed_validation
        fit = treatment_abc(
            row.ecDNA,
            row.chromosomal,
            observed_reference;
            acceptance_fraction = acceptance_fraction
        )

        push!(patient_results, (
            patient_id = row.patient_id,
            true_D = row.true_D,
            true_class = treatment_class(row.true_D),
            median_D = fit.median_D,
            predicted_class = fit.predicted_class,
            p_weak = fit.corrected_weak,
            p_intermediate = fit.corrected_intermediate,
            p_strong = fit.corrected_strong
        ))
    end

    performance = treatment_confusion_matrix(patient_results)

    recovered = Dict(
        D => median(r.median_D for r in patient_results if isapprox(r.true_D, D; atol=1e-12))
        for D in TRUE_D_VALUES
    )

    patient_ids = sort(unique(r.patient_id for r in patient_results))
    full_ordering = 0
    strong_above_weak = 0

    for id in patient_ids
        rows = [r for r in patient_results if r.patient_id == id]
        weak = only(r for r in rows if isapprox(r.true_D, 0.001; atol=1e-12))
        intermediate = only(r for r in rows if isapprox(r.true_D, 0.005; atol=1e-12))
        strong = only(r for r in rows if isapprox(r.true_D, 0.010; atol=1e-12))

        weak.median_D < intermediate.median_D < strong.median_D &&
            (full_ordering += 1)
        strong.median_D > weak.median_D &&
            (strong_above_weak += 1)
    end

    all_ecDNA = vcat((r.ecDNA for r in observed_validation)...)
    baseline_detection = mean(r.ecDNA[1] > 0 for r in observed_validation)

    return (
        depth = N_target,
        threshold = L,
        replicate = replicate,
        balanced_accuracy = performance.balanced_accuracy,
        accuracy = performance.accuracy,
        full_ordering = full_ordering / length(patient_ids),
        strong_above_weak = strong_above_weak / length(patient_ids),
        recovered_weak = recovered[0.001],
        recovered_intermediate = recovered[0.005],
        recovered_strong = recovered[0.010],
        detection_fraction = mean(all_ecDNA .> 0),
        baseline_detection_fraction = baseline_detection,
        confusion = performance.confusion
    )
end

function run_sequencing_sensitivity(
    reference,
    validation;
    depths = DEPTHS,
    thresholds = THRESHOLDS,
    n_tech_reps = 10,
    acceptance_fraction = 0.01
)
    results = NamedTuple[]

    for depth in depths
        for L in thresholds
            n_reps = isinf(depth) ? 1 : n_tech_reps

            for rep in 1:n_reps
                push!(
                    results,
                    analyse_sequencing_condition(
                        reference,
                        validation,
                        depth,
                        L;
                        replicate = rep,
                        acceptance_fraction = acceptance_fraction
                    )
                )
            end
        end
    end

    return results
end

function summarise_sequencing(results)
    summary = NamedTuple[]

    for depth in sort(unique(r.depth for r in results); rev = true)
        for L in sort(unique(r.threshold for r in results))
            rows = [r for r in results if r.depth == depth && r.threshold == L]

            push!(summary, (
                depth = depth,
                threshold = L,
                balanced_accuracy = median(r.balanced_accuracy for r in rows),
                full_ordering = median(r.full_ordering for r in rows),
                strong_above_weak = median(r.strong_above_weak for r in rows),
                baseline_detection = median(r.baseline_detection_fraction for r in rows),
                recovered_weak = median(r.recovered_weak for r in rows),
                recovered_intermediate = median(r.recovered_intermediate for r in rows),
                recovered_strong = median(r.recovered_strong for r in rows)
            ))
        end
    end

    return summary
end
