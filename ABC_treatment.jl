include("treatment.jl")

using Random
using Statistics
using Distributions
using Serialization

const R_PRIOR = Uniform(0.02, 0.06)
const LOGN_PRIOR = Uniform(7.0, 9.0)
const D_PRIOR = Uniform(0.0, 0.012)
const S_MAX_TREATMENT = 1.2
const EARLY_BLOOD_TIMES = [0.0, 0.125, 0.25, 0.5, 1.0, 2.0]

function blood_distance(sim_ecDNA, sim_chr, obs_ecDNA, obs_chr)
    d_ec = mean(abs.(log1p.(sim_ecDNA) .- log1p.(obs_ecDNA)))
    d_chr = mean(abs.(log1p.(sim_chr) .- log1p.(obs_chr)))
    return sqrt((d_ec^2 + d_chr^2) / 2)
end

function treatment_class(D)
    D < 0.003 && return "Weak"
    D < 0.0075 && return "Intermediate"
    return "Strong"
end

function build_treatment_reference(
    n_particles;
    r_prior = R_PRIOR,
    logN_prior = LOGN_PRIOR,
    D_prior = D_PRIOR,
    s_max = S_MAX_TREATMENT,
    blood_times = EARLY_BLOOD_TIMES,
    seed = 20261024
)
    Random.seed!(seed)
    reference = NamedTuple[]
    attempts = 0

    while length(reference) < n_particles
        attempts += 1

        r = rand(r_prior)
        logN = rand(logN_prior)
        N_first = 10.0^logN
        D = rand(D_prior)

        sim = simulate_chemo_from_Nfirst(
            r, s_max, N_first, D;
            blood_draw_times = blood_times,
            track_times = blood_times
        )

        sim === nothing && continue
        (!sim.treatment_started || sim.detection_time != 0.0) && continue

        push!(reference, (
            r = r,
            logN = logN,
            N_first = N_first,
            D = D,
            ecDNA = copy(sim.sampled_ecDNA),
            chromosomal = copy(sim.sampled_chromosomal)
        ))
    end

    return reference, attempts
end

function treatment_abc(
    obs_ecDNA,
    obs_chromosomal,
    reference;
    acceptance_fraction = 0.01
)
    n_accept = max(1, round(Int, acceptance_fraction * length(reference)))

    distances = [
        blood_distance(
            p.ecDNA, p.chromosomal,
            obs_ecDNA, obs_chromosomal
        )
        for p in reference
    ]

    accepted_idx = sortperm(distances)[1:n_accept]
    accepted = reference[accepted_idx]

    accepted_D = [p.D for p in accepted]
    accepted_r = [p.r for p in accepted]
    accepted_logN = [p.logN for p in accepted]

    classes = ("Weak", "Intermediate", "Strong")

    ref_counts = [
        count(p -> treatment_class(p.D) == c, reference)
        for c in classes
    ]
    accepted_counts = [
        count(D -> treatment_class(D) == c, accepted_D)
        for c in classes
    ]

    ref_mass = ref_counts ./ sum(ref_counts)
    posterior_mass = accepted_counts ./ sum(accepted_counts)

    scores = posterior_mass ./ ref_mass
    corrected_mass = scores ./ sum(scores)

    return (
        median_D = median(accepted_D),
        D_q25 = quantile(accepted_D, 0.25),
        D_q75 = quantile(accepted_D, 0.75),
        median_r = median(accepted_r),
        median_logN = median(accepted_logN),
        posterior_weak = posterior_mass[1],
        posterior_intermediate = posterior_mass[2],
        posterior_strong = posterior_mass[3],
        corrected_weak = corrected_mass[1],
        corrected_intermediate = corrected_mass[2],
        corrected_strong = corrected_mass[3],
        predicted_class = classes[argmax(corrected_mass)],
        accepted = accepted,
        accepted_distances = distances[accepted_idx]
    )
end

function simulate_evaluable_treatment_patient(
    r, logN, D;
    s_max = S_MAX_TREATMENT,
    blood_times = EARLY_BLOOD_TIMES
)
    return simulate_chemo_from_Nfirst(
        r, s_max, 10.0^logN, D;
        blood_draw_times = blood_times,
        track_times = blood_times
    )
end

function validate_treatment(
    reference;
    D_values = [0.001, 0.005, 0.010],
    n_patients = 50,
    acceptance_fraction = 0.01,
    r_prior = R_PRIOR,
    logN_prior = LOGN_PRIOR,
    seed = 314159
)
    rng = MersenneTwister(seed)
    results = NamedTuple[]
    attempts = 0
    patient = 0

    while patient < n_patients
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
        patient += 1

        for sim in sims
            fit = treatment_abc(
                sim.ecDNA,
                sim.chromosomal,
                reference;
                acceptance_fraction = acceptance_fraction
            )

            push!(results, (
                patient = patient,
                r_true = r_true,
                logN_true = logN_true,
                D_true = sim.D_true,
                true_class = treatment_class(sim.D_true),
                median_D = fit.median_D,
                D_q25 = fit.D_q25,
                D_q75 = fit.D_q75,
                predicted_class = fit.predicted_class,
                p_weak = fit.corrected_weak,
                p_intermediate = fit.corrected_intermediate,
                p_strong = fit.corrected_strong
            ))
        end
    end

    return results, attempts
end

function treatment_confusion_matrix(results)
    classes = ("Weak", "Intermediate", "Strong")
    confusion = zeros(Int, 3, 3)

    for row in results
        i = findfirst(==(row.true_class), classes)
        j = findfirst(==(row.predicted_class), classes)
        confusion[i, j] += 1
    end

    recalls = [
        sum(confusion[i, :]) == 0 ? NaN :
        confusion[i, i] / sum(confusion[i, :])
        for i in 1:3
    ]

    accuracy = sum(confusion[i, i] for i in 1:3) / sum(confusion)
    balanced_accuracy = mean(filter(x -> !isnan(x), recalls))

    return (
        confusion = confusion,
        accuracy = accuracy,
        balanced_accuracy = balanced_accuracy,
        recalls = recalls
    )
end

function save_reference(filename, reference)
    serialize(filename, reference)
end

function load_reference(filename)
    return deserialize(filename)
end
