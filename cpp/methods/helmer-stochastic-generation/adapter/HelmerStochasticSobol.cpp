#include "HelmerStochasticSobol.h"

#include <stdexcept>
#include <vector>

#include "rng.h"
#include "xor_values.h"

HelmerStochasticSobol::HelmerStochasticSobol(
    unsigned int s,
    unsigned int k,
    long seed)
    : dimension_(s),
      numPoints_(std::uint64_t{1} << k),
      currentIndex_(0),
      points_(numPoints_ * s)
{
    if (s == 0 ||
        s > static_cast<unsigned int>(sampling::MAX_SOBOL_DIM))
    {
        throw std::out_of_range(
            "Helmer stochastic Sobol supports dimensions 1 to 64");
    }

    if (k >= 31)
    {
        throw std::out_of_range(
            "Helmer stochastic Sobol requires k < 31");
    }

    sampling::RNG rng(
        static_cast<std::uint64_t>(seed));

    for (unsigned int d = 0; d < dimension_; ++d)
        points_[d] = rng.GetUniformFloat();

    for (unsigned int logN = 0;
         (std::uint64_t{1} << logN) < numPoints_;
         ++logN)
    {
        const std::uint64_t previousLength =
            std::uint64_t{1} << logN;

        const std::uint64_t numberOfStrata =
            previousLength * 2;

        const double inverseStrata =
            1.0 /
            static_cast<double>(numberOfStrata);

        for (std::uint64_t i = 0;
             i < previousLength &&
             previousLength + i < numPoints_;
             ++i)
        {
            const std::uint64_t newIndex =
                previousLength + i;

            for (unsigned int d = 0;
                 d < dimension_;
                 ++d)
            {
                const std::uint64_t previousIndex =
                    i ^
                    sampling::sobol_xors[d][logN];

                const int previousStratum =
                    static_cast<int>(
                        points_[
                            previousIndex * dimension_ + d]
                        * numberOfStrata);

                const int stratum =
                    previousStratum ^ 1;

                points_[
                    newIndex * dimension_ + d] =
                    (rng.GetUniformFloat() + stratum) *
                    inverseStrata;
            }
        }
    }
}

double HelmerStochasticSobol::coordinate(
    unsigned int j) const
{
    return points_[
        currentIndex_ * dimension_ + j];
}

void HelmerStochasticSobol::next()
{
    if (currentIndex_ + 1 < numPoints_)
        ++currentIndex_;
}

std::uint64_t HelmerStochasticSobol::numPoints() const
{
    return numPoints_;
}
