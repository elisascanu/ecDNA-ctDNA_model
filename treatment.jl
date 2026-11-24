include("model.jl")

chemo_death_rate(copies, d0, D, chemo_active) =
    chemo_active ? d0 + D * copies : d0

function treatment_event_rates(population, b, d0, D, s_max, k_e, k_m, chemo_active)
    events = Tuple{Int,Float64,Float64}[]
    total_rate = 0.0

    for (copies, count) in population
        birth_rate = b * count * s_t(copies, s_max, k_e, k_m)
        death_rate = chemo_death_rate(copies, d0, D, chemo_active) * count
        total_rate += birth_rate + death_rate
        push!(events, (copies, birth_rate, death_rate))
    end

    return events, total_rate
end

function simulate_chemo_from_Nfirst(
    r, s_max, N_first, D;
    q_s = 1.4e-4,
    b = 0.14,
    epsilon = log(2) / 30 * 60 * 24,
    k_e = 8,
    k_m = 4,
    initial_population = Dict(1 => 30),
    max_population = 100_000,
    resample_size = 10_000,
    blood_draw_times = collect(0.0:10.0:80.0),
    track_times = collect(0.0:1.0:80.0),
    sample_fraction = 15 / 5000,
    treatment_duration = 30.0
)
    d0 = b - r
    d0 <= 0 && return nothing

    population = copy(initial_population)
    ctDNA = Dict{Int,Int}(0 => 0, 1 => 0, 2 => 0)
    particle_weight = 1.0
    pre_time = 0.0

    pre_track_times = Float64[]
    pre_effective_population = Float64[]
    pre_mean_copy_number = Float64[]
    pre_fraction_ecDNA_negative = Float64[]

    function record_pre_state!()
        total_cells = sum(values(population))
        push!(pre_track_times, pre_time)
        push!(pre_effective_population, particle_weight * total_cells)
        push!(pre_mean_copy_number, mean_copy_number(population))
        push!(pre_fraction_ecDNA_negative, fraction_ecDNA_negative(population))
    end

    record_pre_state!()
    next_pre_record = 1.0

    # Latent tumour growth before the first blood draw
    while true
        total_cells = sum(values(population))
        total_cells == 0 && return nothing
        particle_weight * total_cells >= N_first && break

        if total_cells >= max_population
            old_population_size = total_cells
            population = resample_population(population, resample_size)
            particle_weight *= old_population_size / resample_size
            continue
        end

        events, total_rate = event_rates(population, b, d0, s_max, k_e, k_m)
        total_rate <= 0 && return nothing

        dt = rand(Exponential(1 / total_rate))
        clear_blood!(ctDNA, dt, epsilon)
        pre_time += dt

        event_type, copies = draw_event(events, total_rate)
        apply_event!(population, ctDNA, event_type, copies, q_s, particle_weight)

        if pre_time >= next_pre_record
            record_pre_state!()
            next_pre_record = pre_time + 1.0
        end
    end

    record_pre_state!()

    track_times = sort(Float64.(track_times))
    blood_draw_times = sort(Float64.(blood_draw_times))
    isempty(track_times) && error("track_times cannot be empty.")
    track_times[1] != 0.0 && error("The first track time must be 0.")

    observation_time = 0.0
    treatment_started = false
    detection_time = NaN
    chemo_start = NaN
    chemo_end = NaN

    track_time_out = Float64[]
    mean_copy_out = Float64[]
    fraction_negative_out = Float64[]
    effective_population_out = Float64[]
    tumour_ecDNA_total_out = Float64[]
    whole_ecDNA_out = Int[]
    whole_chromosomal_out = Int[]
    chemo_active_out = Bool[]

    draw_time_out = Float64[]
    sampled_ecDNA_out = Int[]
    sampled_chromosomal_out = Int[]

    function record_state!(t)
        total_cells = sum(values(population))
        effective_population = particle_weight * total_cells
        mean_copy = mean_copy_number(population)
        fraction_negative = fraction_ecDNA_negative(population)
        chromosomal = ctDNA[0] + ctDNA[1]
        ecDNA = ctDNA[2]
        chemo_active = treatment_started && t >= chemo_start && t < chemo_end

        push!(track_time_out, t)
        push!(mean_copy_out, mean_copy)
        push!(fraction_negative_out, fraction_negative)
        push!(effective_population_out, effective_population)
        push!(tumour_ecDNA_total_out, effective_population * mean_copy)
        push!(whole_ecDNA_out, ecDNA)
        push!(whole_chromosomal_out, chromosomal)
        push!(chemo_active_out, chemo_active)
    end

    function perform_blood_draw!(t)
        chromosomal = ctDNA[0] + ctDNA[1]
        ecDNA = ctDNA[2]

        sampled_chr = rand(Binomial(chromosomal, sample_fraction))
        sampled_ec = rand(Binomial(ecDNA, sample_fraction))

        push!(draw_time_out, t)
        push!(sampled_chromosomal_out, sampled_chr)
        push!(sampled_ecDNA_out, sampled_ec)

        if !treatment_started && sampled_ec > 0
            treatment_started = true
            detection_time = t
            chemo_start = t
            chemo_end = t + treatment_duration
        end
    end

    draw_time_set = Set(blood_draw_times)

    record_state!(0.0)
    0.0 in draw_time_set && perform_blood_draw!(0.0)

    for target_time in track_times[2:end]
        while observation_time < target_time
            boundary = target_time
            chemo_active_now = treatment_started &&
                               observation_time >= chemo_start &&
                               observation_time < chemo_end

            if chemo_active_now && chemo_end < boundary
                boundary = chemo_end
            end

            while observation_time < boundary
                total_cells = sum(values(population))
                total_cells == 0 && return nothing

                if total_cells >= max_population
                    old_population_size = total_cells
                    population = resample_population(population, resample_size)
                    particle_weight *= old_population_size / resample_size
                    continue
                end

                chemo_active_now = treatment_started &&
                                   observation_time >= chemo_start &&
                                   observation_time < chemo_end

                events, total_rate = treatment_event_rates(
                    population, b, d0, D, s_max, k_e, k_m, chemo_active_now
                )
                total_rate <= 0 && return nothing

                dt_event = rand(Exponential(1 / total_rate))
                remaining = boundary - observation_time

                if dt_event >= remaining
                    clear_blood!(ctDNA, remaining, epsilon)
                    observation_time = boundary
                    break
                end

                clear_blood!(ctDNA, dt_event, epsilon)
                observation_time += dt_event

                event_type, copies = draw_event(events, total_rate)
                apply_event!(
                    population, ctDNA, event_type, copies, q_s, particle_weight
                )
            end
        end

        record_state!(target_time)
        target_time in draw_time_set && perform_blood_draw!(target_time)
    end

    return (
        r = r,
        s_max = s_max,
        D = D,
        requested_N_first = N_first,
        actual_N_first = effective_population_out[1],
        latent_tumour_age = pre_time,
        treatment_started = treatment_started,
        detection_time = detection_time,
        chemo_start = chemo_start,
        chemo_end = chemo_end,
        treatment_duration = treatment_duration,
        pre_track_times = pre_track_times,
        pre_effective_population = pre_effective_population,
        pre_mean_copy_number = pre_mean_copy_number,
        pre_fraction_ecDNA_negative = pre_fraction_ecDNA_negative,
        track_times = track_time_out,
        effective_population = effective_population_out,
        mean_copy_number = mean_copy_out,
        fraction_ecDNA_negative = fraction_negative_out,
        tumour_ecDNA_total = tumour_ecDNA_total_out,
        whole_ecDNA = whole_ecDNA_out,
        whole_chromosomal = whole_chromosomal_out,
        chemo_active = chemo_active_out,
        blood_draw_times = draw_time_out,
        sampled_ecDNA = sampled_ecDNA_out,
        sampled_chromosomal = sampled_chromosomal_out
    )
end
