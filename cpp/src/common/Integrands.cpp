#include "Integrands.h"

#include <cmath>

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
