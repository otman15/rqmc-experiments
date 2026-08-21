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


const _SQRT2 = sqrt(2.0)
const _SQRT2PI = sqrt(2.0 * pi)


function _erf(x::Float64)
    return ccall(:erf, Float64, (Float64,), x)
end


function _normal_cdf(x::Float64)
    return 0.5 * (1.0 + _erf(x / _SQRT2))
end


function _normal_quantile(p::Float64)
    p = clamp(p, nextfloat(0.0), prevfloat(1.0))

    a1 = -3.969683028665376e+01
    a2 = 2.209460984245205e+02
    a3 = -2.759285104469687e+02
    a4 = 1.383577518672690e+02
    a5 = -3.066479806614716e+01
    a6 = 2.506628277459239e+00
    b1 = -5.447609879822406e+01
    b2 = 1.615858368580409e+02
    b3 = -1.556989798598866e+02
    b4 = 6.680131188771972e+01
    b5 = -1.328068155288572e+01
    c1 = -7.784894002430293e-03
    c2 = -3.223964580411365e-01
    c3 = -2.400758277161838e+00
    c4 = -2.549732539343734e+00
    c5 = 4.374664141464968e+00
    c6 = 2.938163982698783e+00
    d1 = 7.784695709041462e-03
    d2 = 3.224671290700398e-01
    d3 = 2.445134137142996e+00
    d4 = 3.754408661907416e+00
    p_low = 0.02425
    p_high = 1.0 - p_low

    if p < p_low
        q = sqrt(-2.0 * log(p))
        x = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
            ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0)
    elseif p <= p_high
        q = p - 0.5
        r = q * q
        x = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
            (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0)
    else
        q = sqrt(-2.0 * log1p(-p))
        x = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
             ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0)
    end

    error = _normal_cdf(x) - p
    density = exp(-0.5 * x * x) / _SQRT2PI
    if density > 0.0
        x -= (error / density) / (1.0 + 0.5 * x * error / density)
    end

    return x
end


function _normal_score(u::Real)
    return _normal_quantile(Float64(u))
end


"""
Compute one centered Oscillatory RQMC estimate.

This matches the SSJ Oscillatory model, centered to match the Genz
oscillatory form.
"""
function oscillatory(points::AbstractMatrix{<:Real})
    s, n = size(points)
    mean_product = 1.0
    sum_a = 0.0
    sum_values = 0.0

    @inbounds for j in 1:s
        a = j / s
        mean_product *= 2.0 * sin(a / 2.0) / a
        sum_a += a
    end

    exact_mean = cos(sum_a / 2.0) * mean_product

    @inbounds for i in 1:n
        angle = 0.0

        for j in 1:s
            angle += (j / s) * points[j, i]
        end

        sum_values += cos(angle) - exact_mean
    end

    return sum_values / n
end


"""
Compute one centered Gaussian RQMC estimate.

This matches the SSJ Gaussian model with the Genz peak sign convention.
"""
function gaussian(points::AbstractMatrix{<:Real})
    s, n = size(points)
    exact_mean = (sqrt(pi) / 2.0 * _erf(1.0))^s
    sum_values = 0.0

    @inbounds for i in 1:n
        sum_squares = 0.0

        for j in 1:s
            u = points[j, i]
            sum_squares += u * u
        end

        sum_values += exp(-sum_squares) - exact_mean
    end

    return sum_values / n
end


"""
Compute one SmoothGauss RQMC estimate matching the SSJ benchmark.
"""
function smoothgauss(points::AbstractMatrix{<:Real})
    s, n = size(points)
    inv_sqrt_s = 1.0 / sqrt(s)
    shift = _normal_cdf(1.0 / _SQRT2)
    sum_values = 0.0

    @inbounds for i in 1:n
        normal_sum = 0.0

        for j in 1:s
            normal_sum += _normal_score(points[j, i])
        end

        sum_values += _normal_cdf(1.0 + normal_sum * inv_sqrt_s) - shift
    end

    return sum_values / n
end


"""
Compute one PieceLinGauss RQMC estimate matching the SSJ benchmark.
"""
function piecelingauss(points::AbstractMatrix{<:Real})
    s, n = size(points)
    inv_sqrt_s = 1.0 / sqrt(s)
    shift = _normal_cdf(-1.0) - exp(-0.5) / _SQRT2PI
    sum_values = 0.0

    @inbounds for i in 1:n
        normal_sum = 0.0

        for j in 1:s
            normal_sum += _normal_score(points[j, i])
        end

        sum_values += max(normal_sum * inv_sqrt_s - 1.0, 0.0) + shift
    end

    return sum_values / n
end


"""
Compute one IndSumNormal RQMC estimate matching the SSJ benchmark.
"""
function indsumnormal(points::AbstractMatrix{<:Real})
    s, n = size(points)
    inv_sqrt_s = 1.0 / sqrt(s)
    shift = _normal_cdf(-1.0)
    sum_values = 0.0

    @inbounds for i in 1:n
        normal_sum = 0.0

        for j in 1:s
            normal_sum += _normal_score(points[j, i])
        end

        sum_values += normal_sum * inv_sqrt_s >= 1.0 ? 1.0 - shift : -shift
    end

    return sum_values / n
end
