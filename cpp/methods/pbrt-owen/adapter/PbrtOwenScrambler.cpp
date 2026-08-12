#include "PbrtOwenScrambler.h"

#include "PbrtSobolCore.h"

PbrtOwenScrambler::PbrtOwenScrambler(
    unsigned int s,
    unsigned int k,
    long seed)
    : dimension_(s),
      numPoints_(std::uint64_t{1} << k),
      currentIndex_(0),
      seed_(seed),
      currentPoint_(s)
{
    updateCurrentPoint();
}

std::uint32_t PbrtOwenScrambler::scramble(
    std::uint32_t value,
    std::uint32_t seed)
{
    if (seed & 1u)
        value ^= 1u << 31;

    for (int b = 1; b < 32; ++b)
    {
        const std::uint32_t mask =
            (~0u) << (32 - b);

        if (static_cast<std::uint32_t>(
                PbrtSobolCore::mixBits(
                    (value & mask) ^ seed)) &
            (1u << b))
        {
            value ^= 1u << (31 - b);
        }
    }

    return value;
}

void PbrtOwenScrambler::updateCurrentPoint()
{
    for (unsigned int j = 0; j < dimension_; ++j)
    {
        const std::uint32_t value =
            PbrtSobolCore::sobolInteger(
                currentIndex_,
                j);

        const std::uint32_t scrambleSeed =
            PbrtSobolCore::dimensionSeed(
                j,
                seed_);

        currentPoint_[j] =
            PbrtSobolCore::toUnitInterval(
                scramble(value, scrambleSeed));
    }
}

double PbrtOwenScrambler::coordinate(unsigned int j) const
{
    return currentPoint_[j];
}

void PbrtOwenScrambler::next()
{
    if (currentIndex_ + 1 < numPoints_)
    {
        ++currentIndex_;
        updateCurrentPoint();
    }
}

std::uint64_t PbrtOwenScrambler::numPoints() const
{
    return numPoints_;
}
