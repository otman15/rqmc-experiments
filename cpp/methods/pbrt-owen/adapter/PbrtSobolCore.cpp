#include "PbrtSobolCore.h"

#include <algorithm>
#include <limits>
#include <stdexcept>

#include <pbrt/util/hash.h>
#include <pbrt/util/sobolmatrices.h>

namespace PbrtSobolCore
{

std::uint32_t sobolInteger(
    std::uint64_t index,
    unsigned int dimension)
{
    if (dimension >=
        static_cast<unsigned int>(pbrt::NSobolDimensions))
    {
        throw std::out_of_range(
            "PBRT Sobol dimension exceeds available dimensions");
    }

    if (index >=
        (std::uint64_t{1} << pbrt::SobolMatrixSize))
    {
        throw std::out_of_range(
            "PBRT Sobol index exceeds 52-bit matrix size");
    }

    std::uint32_t value = 0;

    int matrixIndex =
        static_cast<int>(dimension) *
        pbrt::SobolMatrixSize;

    while (index != 0)
    {
        if (index & 1)
            value ^= pbrt::SobolMatrices32[matrixIndex];

        index >>= 1;
        ++matrixIndex;
    }

    return value;
}

std::uint32_t dimensionSeed(
    unsigned int dimension,
    long seed)
{
    if (seed < static_cast<long>(
                   std::numeric_limits<int>::min()) ||
        seed > static_cast<long>(
                   std::numeric_limits<int>::max()))
    {
        throw std::out_of_range(
            "PBRT seed must fit in int");
    }

    const int pbrtDimension =
        static_cast<int>(dimension);

    const int pbrtSeed =
        static_cast<int>(seed);

    return static_cast<std::uint32_t>(
        pbrt::Hash(pbrtDimension, pbrtSeed));
}

std::uint64_t mixBits(std::uint64_t value)
{
    return pbrt::MixBits(value);
}

std::uint32_t reverseBits32(std::uint32_t value)
{
    value = (value << 16) | (value >> 16);

    value =
        ((value & 0x00ff00ffu) << 8) |
        ((value & 0xff00ff00u) >> 8);

    value =
        ((value & 0x0f0f0f0fu) << 4) |
        ((value & 0xf0f0f0f0u) >> 4);

    value =
        ((value & 0x33333333u) << 2) |
        ((value & 0xccccccccu) >> 2);

    value =
        ((value & 0x55555555u) << 1) |
        ((value & 0xaaaaaaaau) >> 1);

    return value;
}

double toUnitInterval(std::uint32_t value)
{
    const float sample =
        std::min(
            value * 0x1p-32f,
            0x1.fffffep-1f);

    return static_cast<double>(sample);
}

}
