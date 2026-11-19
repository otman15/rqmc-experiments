#include "HelmerNativeSobol.h"

#include <limits>
#include <stdexcept>

#include "ssobol.h"
#include "xor_values.h"

namespace
{
std::uint64_t checkedPointCount(unsigned int s, unsigned int k)
{
    if (s == 0 || s > static_cast<unsigned int>(sampling::MAX_SOBOL_DIM))
        throw std::out_of_range("Native Helmer Sobol supports dimensions 1 to 64");
    if (k >= 31)
        throw std::out_of_range("Native Helmer Sobol requires k < 31");

    const auto n = std::uint64_t{1} << k;
    // Upstream uses signed int arithmetic for point * dimension offsets.
    if (n * s > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
        throw std::length_error("Native Helmer Sobol requires n * s to fit int");
    if (n * s > std::vector<double>().max_size())
        throw std::length_error("Native Helmer Sobol output buffer is too large");
    return n;
}
}

HelmerNativeSobol::HelmerNativeSobol(unsigned int s, unsigned int k)
    : dimension_(s),
      numPoints_(checkedPointCount(s, k)),
      points_(numPoints_ * s)
{
}

void HelmerNativeSobol::randomize()
{
    sampling::GetStochasticSobolSamples(
        static_cast<int>(numPoints_), static_cast<int>(dimension_),
        false, 1, false, points_.data());
    reset();
}

void HelmerNativeSobol::reset()
{
    currentIndex_ = 0;
}

double HelmerNativeSobol::coordinate(unsigned int j) const
{
    return points_[currentIndex_ * dimension_ + j];
}

void HelmerNativeSobol::next()
{
    if (currentIndex_ + 1 < numPoints_)
        ++currentIndex_;
}

std::uint64_t HelmerNativeSobol::numPoints() const
{
    return numPoints_;
}
