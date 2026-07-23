using QuasiMonteCarlo
using Random
using Statistics
using Printf

include(joinpath(@__DIR__, "integrands.jl"))

# Experiment parameters
const s_values = [2, 4]
const k_values = [8, 10]

const num_reps = 10_000
const pad = 30
const seed = 12345


"""
Corrected excess kurtosis corresponding to:

e1071::kurtosis(x, type = 2)
"""
function corrected_excess_kurtosis(values::AbstractVector{<:Real})
    n = length(values)

    n >= 4 || error("At least four estimates are required.")

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

    moment2 > 0.0 || error("Kurtosis is undefined when all estimates are equal.")

    g2 = moment4 / (moment2 * moment2) - 3.0

    return ((n - 1) * ((n + 1) * g2 + 6.0)) /
           ((n - 2) * (n - 3))
end


function run_experiments()
    output_path = joinpath(@__DIR__, "..", "results", "scMl_res.csv")
    mkpath(dirname(output_path))

    rows = Vector{Tuple{Int, Int, String, Float64, Float64, Float64}}()

    println("s,k,method,variance,kurtosis,cpu_time")

    for s in s_values
        for k in k_values
            n = 2^k

            experiment_elapsed = @elapsed begin
                # Restart the random stream for each (s, k) experiment.
                rng = MersenneTwister(seed)

                sampler = SobolSample(
                    R = OwenScramble(
                        base = 2,
                        pad = pad,
                        rng = rng
                    )
                )

                # SciML generates the deterministic Sobol point set once and
                # applies num_reps independent NUS randomizations to that set.
                randomized_point_sets =
                    DesignMatrix(n, s, sampler, num_reps)

                estimates = Vector{Float64}(undef, num_reps)

                for (rep, points) in enumerate(randomized_point_sets)
                    estimates[rep] = mc2(points)
                end
            end

            variance = var(estimates)
            kurtosis = corrected_excess_kurtosis(estimates)

            push!(rows, (s, k, "SciML", variance, kurtosis, experiment_elapsed))

            @printf(
                "%d,%d,%s,%.15e,%.15e,%.6f\n",
                s,
                k,
                "SciML",
                variance,
                kurtosis,
                experiment_elapsed
            )
        end
    end

    open(output_path, "w") do io
        println(io, "s,k,method,variance,kurtosis,cpu_time")
        for (s, k, method, variance, kurtosis, cpu_time) in rows
            @printf(
                io,
                "%d,%d,%s,%.15e,%.15e,%.6f\n",
                s,
                k,
                method,
                variance,
                kurtosis,
                cpu_time
            )
        end
    end

    println("Wrote results to ", output_path)
end


run_experiments()