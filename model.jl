using Random
using Distributions

# Copy-number-dependent selection
function s_t(n, s_max, k_e, k_m)
    n == 0 && return 1.0
    n >= k_e && return s_max
    return 1 + (s_max - 1) *
           (1 + (k_e - n) / (k_e - k_m)) *
           (n / k_e)^(k_e / (k_e - k_m))
end

# Cell division with random ecDNA segregation
function birth!(population::Dict{Int,Int}, copies::Int)
    population[copies] -= 1
    population[copies] == 0 && delete!(population, copies)

    total_copies = 2 * copies
    daughter1 = rand(Binomial(total_copies, 0.5))
    daughter2 = total_copies - daughter1

    population[daughter1] = get(population, daughter1, 0) + 1
    population[daughter2] = get(population, daughter2, 0) + 1
    return population
end

# Sample cells without replacement across copy-number classes
function resample_population(population::Dict{Int,Int}, n_sample::Int)
    total_cells = sum(values(population))
    n_sample > total_cells && error("n_sample cannot exceed population size.")

    copy_classes = sort(collect(keys(population)))
    sampled_population = Dict{Int,Int}()
    remaining_total = total_cells
    remaining_sample = n_sample

    for copies in copy_classes[1:end-1]
        count = population[copies]
        sampled = rand(Hypergeometric(count, remaining_total - count, remaining_sample))
        sampled > 0 && (sampled_population[copies] = sampled)
        remaining_sample -= sampled
        remaining_total -= count
    end

    if remaining_sample > 0
        last_class = copy_classes[end]
        sampled_population[last_class] =
            get(sampled_population, last_class, 0) + remaining_sample
    end

    return sampled_population
end

function clear_blood!(ctDNA::Dict{Int,Int}, dt, epsilon)
    dt <= 0 && return ctDNA
    p_survival = clamp(exp(-epsilon * dt), 0.0, 1.0)

    for i in 0:2
        ctDNA[i] > 0 && (ctDNA[i] = rand(Binomial(ctDNA[i], p_survival)))
    end
    return ctDNA
end

# ctDNA[0]: chromosomal DNA from ecDNA-negative cells
# ctDNA[1]: chromosomal DNA from ecDNA-positive cells
# ctDNA[2]: circulating ecDNA copies
function weighted_shedding!(ctDNA::Dict{Int,Int}, copies::Int, q_s, particle_weight)
    biological_deaths = max(1, round(Int, particle_weight))
    chromosomal_shed = rand(Binomial(biological_deaths, q_s))

    if copies == 0
        ctDNA[0] += chromosomal_shed
    else
        ctDNA[1] += chromosomal_shed
        ctDNA[2] += rand(Binomial(biological_deaths * copies, q_s))
    end
    return ctDNA
end

function mean_copy_number(population::Dict{Int,Int})
    N = sum(values(population))
    N == 0 && return NaN
    return sum(copies * count for (copies, count) in population) / N
end

function fraction_ecDNA_negative(population::Dict{Int,Int})
    N = sum(values(population))
    N == 0 && return NaN
    return get(population, 0, 0) / N
end

function event_rates(population, b, d, s_max, k_e, k_m)
    events = Tuple{Int,Float64,Float64}[]
    total_rate = 0.0

    for (copies, count) in population
        birth_rate = b * count * s_t(copies, s_max, k_e, k_m)
        death_rate = d * count
        total_rate += birth_rate + death_rate
        push!(events, (copies, birth_rate, death_rate))
    end

    return events, total_rate
end

function draw_event(events, total_rate)
    x = rand() * total_rate
    cumulative = 0.0

    for (copies, birth_rate, death_rate) in events
        cumulative += birth_rate
        x <= cumulative && return :birth, copies

        cumulative += death_rate
        x <= cumulative && return :death, copies
    end

    error("No event selected.")
end

function apply_event!(population, ctDNA, event_type, copies, q_s, particle_weight)
    if event_type == :birth
        birth!(population, copies)
    else
        population[copies] -= 1
        population[copies] == 0 && delete!(population, copies)
        weighted_shedding!(ctDNA, copies, q_s, particle_weight)
    end
end

function simulate_longitudinal_blood_from_Nfirst(
    r, s_max, N_first;
    q_s = 1.4e-4,
    b = 0.14,
    epsilon = log(2) / 30 * 60 * 24,
    k_e = 8,
    k_m = 4,
    initial_population = Dict(1 => 30),
    max_population = 100_000,
    resample_size = 10_000,
    sample_times = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0],
    sample_fraction = 15 / 5000
)
    d = b - r
    d <= 0 && return nothing

    population = copy(initial_population)
    ctDNA = Dict{Int,Int}(0 => 0, 1 => 0, 2 => 0)
    particle_weight = 1.0
    pre_time = 0.0

    # Latent growth up to the first blood draw
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

        events, total_rate = event_rates(population, b, d, s_max, k_e, k_m)
        total_rate <= 0 && return nothing

        dt = rand(Exponential(1 / total_rate))
        clear_blood!(ctDNA, dt, epsilon)
        pre_time += dt

        event_type, copies = draw_event(events, total_rate)
        apply_event!(population, ctDNA, event_type, copies, q_s, particle_weight)
    end

    sample_times = sort(Float64.(sample_times))
    isempty(sample_times) && error("sample_times cannot be empty.")
    sample_times[1] != 0.0 && error("The first sample time must be 0.")

    observation_time = 0.0
    whole_chromosomal = Int[]
    whole_ecDNA = Int[]
    sampled_chromosomal = Int[]
    sampled_ecDNA = Int[]
    effective_population_samples = Float64[]
    mean_copy_samples = Float64[]
    fraction_zero_samples = Float64[]
    particle_weight_samples = Float64[]

    function record_sample!()
        chromosomal = ctDNA[0] + ctDNA[1]
        ecDNA = ctDNA[2]

        push!(whole_chromosomal, chromosomal)
        push!(whole_ecDNA, ecDNA)
        push!(sampled_chromosomal, rand(Binomial(chromosomal, sample_fraction)))
        push!(sampled_ecDNA, rand(Binomial(ecDNA, sample_fraction)))

        total_cells = sum(values(population))
        push!(effective_population_samples, particle_weight * total_cells)
        push!(mean_copy_samples, mean_copy_number(population))
        push!(fraction_zero_samples, fraction_ecDNA_negative(population))
        push!(particle_weight_samples, particle_weight)
    end

    record_sample!()

    for target_time in sample_times[2:end]
        while observation_time < target_time
            total_cells = sum(values(population))
            total_cells == 0 && return nothing

            if total_cells >= max_population
                old_population_size = total_cells
                population = resample_population(population, resample_size)
                particle_weight *= old_population_size / resample_size
                continue
            end

            events, total_rate = event_rates(population, b, d, s_max, k_e, k_m)
            total_rate <= 0 && return nothing

            dt_event = rand(Exponential(1 / total_rate))
            remaining_time = target_time - observation_time

            if dt_event >= remaining_time
                clear_blood!(ctDNA, remaining_time, epsilon)
                observation_time = target_time
                break
            end

            clear_blood!(ctDNA, dt_event, epsilon)
            observation_time += dt_event

            event_type, copies = draw_event(events, total_rate)
            apply_event!(population, ctDNA, event_type, copies, q_s, particle_weight)
        end

        record_sample!()
    end

    return (
        sample_times = sample_times,
        chromosomal = sampled_chromosomal,
        ecDNA = sampled_ecDNA,
        whole_chromosomal = whole_chromosomal,
        whole_ecDNA = whole_ecDNA,
        effective_population = effective_population_samples,
        mean_copy_number = mean_copy_samples,
        fraction_ecDNA_negative = fraction_zero_samples,
        particle_weight = particle_weight_samples,
        latent_tumour_age = pre_time,
        actual_N_first = effective_population_samples[1]
    )
end
