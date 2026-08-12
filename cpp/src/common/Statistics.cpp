#include "Statistics.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

double mean(const std::vector<double>& values)
{
    if (values.empty())
        return std::numeric_limits<double>::quiet_NaN();

    double sum = 0.0;

    for (double value : values)
        sum += value;

    return sum / static_cast<double>(values.size());
}

double median(const std::vector<double>& values)
{
    if (values.empty())
        return std::numeric_limits<double>::quiet_NaN();

    std::vector<double> sorted = values;
    std::sort(sorted.begin(), sorted.end());

    const std::size_t n = sorted.size();

    if (n % 2 == 1)
        return sorted[n / 2];

    return 0.5 * (sorted[n / 2 - 1] + sorted[n / 2]);
}

double sampleVariance(const std::vector<double>& values)
{
    const std::size_t n = values.size();

    if (n < 2)
        return std::numeric_limits<double>::quiet_NaN();

    const double avg = mean(values);

    double sumSquaredDeviations = 0.0;

    for (double value : values)
    {
        const double deviation = value - avg;
        sumSquaredDeviations += deviation * deviation;
    }

    return sumSquaredDeviations /
           static_cast<double>(n - 1);
}

double skewness(const std::vector<double>& values)
{
    const std::size_t n = values.size();

    if (n < 3)
        return std::numeric_limits<double>::quiet_NaN();

    const double avg = mean(values);

    double m2 = 0.0;
    double m3 = 0.0;

    for (double value : values)
    {
        const double deviation = value - avg;
        const double deviation2 = deviation * deviation;

        m2 += deviation2;
        m3 += deviation2 * deviation;
    }

    const double nd = static_cast<double>(n);

    m2 /= nd;
    m3 /= nd;

    if (m2 == 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    const double g1 =
        m3 / std::pow(m2, 1.5);

    return std::sqrt(nd * (nd - 1.0)) /
           (nd - 2.0) *
           g1;
}

double excessKurtosis(const std::vector<double>& values)
{
    const std::size_t n = values.size();

    if (n < 4)
        return std::numeric_limits<double>::quiet_NaN();

    const double avg = mean(values);

    double m2 = 0.0;
    double m4 = 0.0;

    for (double value : values)
    {
        const double deviation = value - avg;
        const double deviation2 = deviation * deviation;

        m2 += deviation2;
        m4 += deviation2 * deviation2;
    }

    const double nd = static_cast<double>(n);

    m2 /= nd;
    m4 /= nd;

    if (m2 == 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    const double g2 =
        m4 / (m2 * m2) - 3.0;

    return
        ((nd - 1.0) /
         ((nd - 2.0) * (nd - 3.0))) *
        ((nd + 1.0) * g2 + 6.0);
}
