using QuasiMonteCarlo
using Random
using Statistics
using Printf

include(joinpath(@__DIR__, "integrands.jl"))

const s_values = [2, 4, 8, 16, 32]
const k_values = [8, 10, 12, 14, 16]
const m_values = [1, 2, 5, 10, 20, 30, 50]

const pad = 30
const seed = 20171215
const num_runs = 11

# output_dir = joinpath(@__DIR__, "jl-results") 
output_dir = joinpath(@__DIR__, "/home/otman/Documents/GitHub/Data/o-test/nus/nus-timing")
mkpath(output_dir)

results = Vector{Tuple{String, Int, Int, Int, Float64}}()

for s in s_values
    for k in k_values
        for m in m_values
            cpu_times = Vector{Float64}(undef, num_runs)

            for run in 1:num_runs
                elapsed = @elapsed begin
                    rng = MersenneTwister(seed + run)

                    sampler = SobolSample(
                        R = OwenScramble(
                            base = 2,
                            pad = pad,
                            rng = rng
                        )
                    )

                    n = 2^k
                    randomized_point_sets = DesignMatrix(n, s, sampler, m)
                    estimates = Vector{Float64}(undef, m)

                    for (rep, points) in enumerate(randomized_point_sets)
                        estimates[rep] = mc2(points)
                    end
                end

                cpu_times[run] = elapsed
            end

            median_time = median(cpu_times)
            push!(results, ("SciMl", s, k, m, median_time))

            @printf(
                "s=%d, k=%d, m=%d, median time=%.6f seconds\n",
                s,
                k,
                m,
                median_time
            )
        end
    end
end

output_file = joinpath(output_dir, "sciMl_time_by_rep-$(num_runs).csv")

open(output_file, "w") do io
    println(io, "method,s,k,m,median_time")
    for (method, s, k, m, median_time) in results
        @printf(io, "%s,%d,%d,%d,%.12f\n", method, s, k, m, median_time)
    end
end

println("Results written to ", output_file)
