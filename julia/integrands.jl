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


"""
Compute one SumUeU RQMC estimate from an s × n matrix of points.
"""
function sumueu(points::AbstractMatrix{<:Real})
    s, n = size(points)
    sum_values = 0.0

    @inbounds for i in 1:n
        value = -s

        for j in 1:s
            u = points[j, i]
            value += u * exp(u)
        end

        sum_values += value
    end

    return sum_values / n
end


"""
Compute one SmoothPerB4 RQMC estimate from an s × n matrix of points.
"""
function smoothperb4(points::AbstractMatrix{<:Real}, expw::Real = 1.0)
    s, n = size(points)
    sum_values = 0.0

    @inbounds for i in 1:n
        product = 1.0

        for j in 1:s
            u = points[j, i]
            omega = j^(-expw)
            product *= 1.0 + (30.0 * u^2 * (1.0 - u)^2 - 1.0) * omega
        end

        sum_values += product - 1.0
    end

    return sum_values / n
end


"""
Compute one Polynomial RQMC estimate from an s × n matrix of points.
"""
function polynomial(points::AbstractMatrix{<:Real})
    s, n = size(points)
    sum_values = 0.0

    @inbounds for i in 1:n
        product = 1.0

        for j in 1:s
            product *= 1.0 + (j / s) * (points[j, i] - 0.5)
        end

        sum_values += product - 1.0
    end

    return sum_values / n
end
