using QuasiMonteCarlo
using Random
using Statistics
using Printf
using CPUTime

include(joinpath(@__DIR__, "integrands.jl"))


function run_rqmc_replications(integrand, k, s, m, pad, first_seed)
    n = 2^k
    
    elapsed_time = @CPUelapsed begin
        estimates = Vector{Float64}(undef, m)
        rng = MersenneTwister(first_seed)
        sampler = SobolSample(
            R = OwenScramble(
                base = 2,
                pad = pad,
                rng = rng,
            ),
        )
        randomized_point_sets = DesignMatrix(n, s, sampler, m)

        for (replication, points) in enumerate(randomized_point_sets)
            estimates[replication] = integrand(points)
        end
    end

    return estimates, elapsed_time
end


function corrected_skewness(values::AbstractVector{<:Real})
    n = length(values)
    n >= 3 || return NaN

    mean_value = mean(values)
    moment2 = 0.0
    moment3 = 0.0

    @inbounds for value in values
        deviation = value - mean_value
        deviation2 = deviation * deviation
        moment2 += deviation2
        moment3 += deviation2 * deviation
    end

    moment2 /= n
    moment3 /= n
    moment2 > 0.0 || return NaN

    raw_skewness = moment3 / moment2^(3 / 2)
    return sqrt(n * (n - 1)) / (n - 2) * raw_skewness
end


function corrected_excess_kurtosis(values::AbstractVector{<:Real})
    n = length(values)
    n >= 4 || return NaN

    mean_value = mean(values)
    moment2 = 0.0
    moment4 = 0.0

    @inbounds for value in values
        deviation = value - mean_value
        deviation2 = deviation * deviation
        moment2 += deviation2
        moment4 += deviation2 * deviation2
    end

    moment2 /= n
    moment4 /= n
    moment2 > 0.0 || return NaN

    raw_excess_kurtosis = moment4 / moment2^2 - 3.0
    return ((n - 1) * ((n + 1) * raw_excess_kurtosis + 6.0)) /
           ((n - 2) * (n - 3))
end


function run_timing_experiments(
    models,
    s_values,
    k_values,
    timing_m_values,
    pad,
    seed,
    num_runs,
    output_dir,
)
    mkpath(output_dir)

    for (model_name, integrand) in models
        results = Vector{Tuple{String, Int, Int, Int, Float64}}()

        for s in s_values
            for k in k_values
                for m in timing_m_values
                    elapsed_times = Vector{Float64}(undef, num_runs)

                    for run in 1:num_runs
                        run_first_seed = seed + (run - 1)
                        _, elapsed_times[run] = run_rqmc_replications(
                            integrand,
                            k,
                            s,
                            m,
                            pad,
                            run_first_seed,
                        )
                    end

                    median_time = median(elapsed_times)
                    push!(results, ("SciMl", s, k, m, median_time))

                    @printf(
                        "model=%s, s=%d, k=%d, m=%d, median time=%.6f seconds\n",
                        model_name,
                        s,
                        k,
                        m,
                        median_time,
                    )
                end
            end
        end

        output_file = joinpath(
            output_dir,
            "sciMl_$(model_name)_time_by_rep-$(num_runs).csv",
        )

        open(output_file, "w") do io
            println(io, "method,s,k,m,median_time")
            for (method, s, k, m, median_time) in results
                @printf(io, "%s,%d,%d,%d,%.12f\n", method, s, k, m, median_time)
            end
        end

        println("Results written to ", output_file)
    end
end


function run_statistics_experiments(
    models,
    s_values,
    k_values,
    statistics_m,
    pad,
    seed,
    method,
    output_dir,
)
    mkpath(output_dir)

    for (model_name, integrand) in models
        for s in s_values
            for k in k_values
                estimates, elapsed_time = run_rqmc_replications(
                    integrand,
                    k,
                    s,
                    statistics_m,
                    pad,
                    seed,
                )
                mean_estimate = mean(estimates)
                variance = statistics_m > 1 ? var(estimates) : NaN
                skewness = corrected_skewness(estimates)
                excess_kurtosis = corrected_excess_kurtosis(estimates)

                @printf(
                    "model=%s, s=%d, k=%d, m=%d, mean=%.16e, variance=%.16e, skewness=%.16e, excess kurtosis=%.16e, time=%.6f seconds\n",
                    model_name,
                    s,
                    k,
                    statistics_m,
                    mean_estimate,
                    variance,
                    skewness,
                    excess_kurtosis,
                    elapsed_time,
                )

                output_file = joinpath(
                    output_dir,
                    "$(model_name)-$(s)-$(method)-$(k)-$(statistics_m).dat",
                )

                open(output_file, "w") do io
                    for estimate in sort(estimates)
                        @printf(io, "%.17g\n", estimate)
                    end
                end

                println("Results written to ", output_file)
            end
        end
    end
end


function main()
    TIMING = "TIMING"
    STATISTICS = "STATISTICS"
    mode = TIMING

    method = "NUS-SciML"
    models = [
        ("MC2", mc2),
        ("SumUeU", sumueu),
        ("SmoothPerB4", points -> smoothperb4(points, 1.0)),
        ("Polynomial", polynomial),
        ("Oscillatory", oscillatory),
        ("Gaussian", gaussian),
        ("SmoothGauss", smoothgauss),
        ("PieceLinGauss", piecelingauss),
        ("IndSumNormal", indsumnormal),
    ]

    modelsTiming = [
        ("MC2", mc2)
    ]
    s_values = [2, 4, 8, 16, 32]
    k_values = [8, 10, 12, 14, 16]
    timing_m_values = [1, 2, 5, 10, 20, 30, 50]
    statistics_m = 1000
    pad = 30
    seed = 20171215
    num_runs = 11

    # timing_output_dir = joinpath(@__DIR__, "jl-results")
    timing_output_dir =  joinpath(@__DIR__, "/mnt/c/Users/cherrato/Dropbox/Nus-comparisons/all-methods/") #joinpath(@__DIR__, "jl-results")
    statistics_output_dir = joinpath(timing_output_dir, "dat")

    if mode == TIMING
        run_timing_experiments(
            modelsTiming,
            s_values,
            k_values,
            timing_m_values,
            pad,
            seed,
            num_runs,
            timing_output_dir,
        )
    elseif mode == STATISTICS
        run_statistics_experiments(
            models,
            s_values,
            k_values,
            statistics_m,
            pad,
            seed,
            method,
            statistics_output_dir,
        )
    else
        error("mode must be TIMING or STATISTICS")
    end
end


main()
