#ifndef PBRT_SOBOL_CORE_H
#define PBRT_SOBOL_CORE_H

#include <cstdint>

namespace PbrtSobolCore
{

std::uint32_t sobolInteger(
    std::uint64_t index,
    unsigned int dimension);

std::uint32_t dimensionSeed(
    unsigned int dimension,
    long seed);

std::uint64_t mixBits(std::uint64_t value);

std::uint32_t reverseBits32(std::uint32_t value);

double toUnitInterval(std::uint32_t value);

}

#endif
