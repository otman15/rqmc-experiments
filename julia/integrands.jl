"""
Compute one MC2 RQMC estimate from an s × n matrix of points.

This is equivalent to:
mean(prod((s - x) / (s - 0.5), dimensions) - 1)
"""
function mc2(points::AbstractMatrix{<:Real})
    s, n = size(points)
    scale = 1.0 / (s - 0.5)
    sum_values = 0.0

    @inbounds for i in 1:n
        product = 1.0

        for j in 1:s
            product *= (s - points[j, i]) * scale
        end

        sum_values += product - 1.0
    end

    return sum_values / n
end