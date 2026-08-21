#include "Integrands.h"

#include <cmath>

namespace
{
constexpr double sqrt2 = 1.41421356237309504880;
constexpr double sqrt2Pi = 2.50662827463100050242;
constexpr double sqrtPi = 1.77245385090551602730;

double normalCdf(double x)
{
    return 0.5 * (1.0 + std::erf(x / sqrt2));
}

double normalQuantile(double p)
{
    const double low = std::nextafter(0.0, 1.0);
    const double high = std::nextafter(1.0, 0.0);

    if (p < low)
        p = low;
    else if (p > high)
        p = high;

    const double a1 = -3.969683028665376e+01;
    const double a2 = 2.209460984245205e+02;
    const double a3 = -2.759285104469687e+02;
    const double a4 = 1.383577518672690e+02;
    const double a5 = -3.066479806614716e+01;
    const double a6 = 2.506628277459239e+00;
    const double b1 = -5.447609879822406e+01;
    const double b2 = 1.615858368580409e+02;
    const double b3 = -1.556989798598866e+02;
    const double b4 = 6.680131188771972e+01;
    const double b5 = -1.328068155288572e+01;
    const double c1 = -7.784894002430293e-03;
    const double c2 = -3.223964580411365e-01;
    const double c3 = -2.400758277161838e+00;
    const double c4 = -2.549732539343734e+00;
    const double c5 = 4.374664141464968e+00;
    const double c6 = 2.938163982698783e+00;
    const double d1 = 7.784695709041462e-03;
    const double d2 = 3.224671290700398e-01;
    const double d3 = 2.445134137142996e+00;
    const double d4 = 3.754408661907416e+00;
    const double pLow = 0.02425;
    const double pHigh = 1.0 - pLow;

    double x;

    if (p < pLow)
    {
        const double q = std::sqrt(-2.0 * std::log(p));
        x =
            (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6)
            / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
    }
    else if (p <= pHigh)
    {
        const double q = p - 0.5;
        const double r = q * q;
        x =
            (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6)
            * q
            / (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r
               + 1.0);
    }
    else
    {
        const double q = std::sqrt(-2.0 * std::log1p(-p));
        x =
            -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q
              + c6)
            / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
    }

    const double error = normalCdf(x) - p;
    const double density = std::exp(-0.5 * x * x) / sqrt2Pi;

    if (density > 0.0)
    {
        const double ratio = error / density;
        x -= ratio / (1.0 + 0.5 * x * ratio);
    }

    return x;
}
}

double sumueu(const std::vector<double>& x)
{
    double value = 0.0;

    for (double xi : x)
        value += xi * std::exp(xi);

    return value - static_cast<double>(x.size());
}

double smoothperb4(
    const std::vector<double>& x,
    double expw)
{
    double value = 1.0;

    for (std::size_t j = 0; j < x.size(); ++j)
    {
        double omega =
            std::pow(static_cast<double>(j + 1), -expw);

        double xi = x[j];

        double term =
            1.0
            + (30.0 * xi * xi
               * (1.0 - xi) * (1.0 - xi)
               - 1.0)
              * omega;

        value *= term;
    }

    return value - 1.0;
}

double oscillatory(const std::vector<double>& x)
{
    const double s = static_cast<double>(x.size());
    double meanProduct = 1.0;
    double sumA = 0.0;
    double angle = 0.0;

    for (std::size_t j = 0; j < x.size(); ++j)
    {
        const double a = static_cast<double>(j + 1) / s;
        meanProduct *= 2.0 * std::sin(a / 2.0) / a;
        sumA += a;
        angle += a * x[j];
    }

    const double exactMean =
        std::cos(sumA / 2.0) * meanProduct;

    return std::cos(angle) - exactMean;
}

double gaussian(const std::vector<double>& x)
{
    const double s = static_cast<double>(x.size());
    const double exactMean =
        std::pow(sqrtPi / 2.0 * std::erf(1.0), s);
    double sumSquares = 0.0;

    for (double xi : x)
        sumSquares += xi * xi;

    return std::exp(-sumSquares) - exactMean;
}

double smoothgauss(const std::vector<double>& x)
{
    const double invSqrtS =
        1.0 / std::sqrt(static_cast<double>(x.size()));
    const double shift = normalCdf(1.0 / sqrt2);
    double sumNormals = 0.0;

    for (double xi : x)
        sumNormals += normalQuantile(xi);

    return normalCdf(1.0 + sumNormals * invSqrtS) - shift;
}

double piecelingauss(const std::vector<double>& x)
{
    const double invSqrtS =
        1.0 / std::sqrt(static_cast<double>(x.size()));
    const double shift =
        normalCdf(-1.0) - std::exp(-0.5) / sqrt2Pi;
    double sumNormals = 0.0;

    for (double xi : x)
        sumNormals += normalQuantile(xi);

    const double value = sumNormals * invSqrtS - 1.0;
    return (value < 0.0 ? 0.0 : value) + shift;
}

double indsumnormal(const std::vector<double>& x)
{
    const double invSqrtS =
        1.0 / std::sqrt(static_cast<double>(x.size()));
    const double shift = normalCdf(-1.0);
    double sumNormals = 0.0;

    for (double xi : x)
        sumNormals += normalQuantile(xi);

    return sumNormals * invSqrtS >= 1.0
        ? 1.0 - shift
        : -shift;
}

double mc2(const std::vector<double>& x)
{
    double s = static_cast<double>(x.size());
    double value = 1.0;

    for (double xi : x)
        value *= (s - xi) / (s - 0.5);

    return value - 1.0;
}

double polynomial(const std::vector<double>& x)
{
    double s = static_cast<double>(x.size());
    double value = 1.0;

    for (std::size_t j = 0; j < x.size(); ++j)
    {
        double a = static_cast<double>(j + 1) / s;
        value *= 1.0 + (x[j] - 0.5) * a;
    }

    return value - 1.0;
}
