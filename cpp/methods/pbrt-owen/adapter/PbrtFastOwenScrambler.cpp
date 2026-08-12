#include "PbrtFastOwenScrambler.h"

#include "PbrtSobolCore.h"

PbrtFastOwenScrambler::PbrtFastOwenScrambler(
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

std::uint32_t PbrtFastOwenScrambler::scramble(
    std::uint32_t value,
    std::uint32_t seed)
{
    value = PbrtSobolCore::reverseBits32(value);

    value ^= value * 0x3d20adeau;
    value += seed;
    value *= (seed >> 16) | 1u;
    value ^= value * 0x05526c56u;
    value ^= value * 0x53a22864u;

    return PbrtSobolCore::reverseBits32(value);
}

void PbrtFastOwenScrambler::updateCurrentPoint()
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

double PbrtFastOwenScrambler::coordinate(unsigned int j) const
{
    return currentPoint_[j];
}

void PbrtFastOwenScrambler::next()
{
    if (currentIndex_ + 1 < numPoints_)
    {
        ++currentIndex_;
        updateCurrentPoint();
    }
}

std::uint64_t PbrtFastOwenScrambler::numPoints() const
{
    return numPoints_;
}
