#include "AdaptiveRegularTiles.h"

#include "PbrtSobolCore.h"

#include <limits>
#include <stdexcept>

#include <pbrt/util/rng.h>

namespace
{

struct Production
{
    std::uint16_t left;
    std::uint16_t right;

    std::uint16_t operator[](int i) const
    {
        return i == 0 ? left : right;
    }
};

constexpr unsigned int idCount = 32;
constexpr unsigned int maxDimensions = 1024;

constexpr Production production32[idCount] = {
    {30, 31},
    {22, 23},
    {24, 25},
    {26, 27},
    {28, 29},
    {10, 11},
    {12, 13},
    {14, 15},
    {16, 17},
    { 6,  7},
    {18, 19},
    { 2,  3},
    {20, 21},
    {22, 23},
    {24, 25},
    {26, 27},
    {28, 29},
    {10, 11},
    {30, 31},
    {22, 23},
    { 0,  1},
    { 2,  3},
    { 4,  5},
    { 6,  7},
    { 8,  9},
    {10, 11},
    {12, 13},
    {14, 15},
    {16, 17},
    { 6,  7},
    { 0,  1},
    { 2,  3}
};

}

AdaptiveRegularTiles::AdaptiveRegularTiles(
    unsigned int s,
    unsigned int k,
    long seed)
    : dimension_(s),
      numPoints_(std::uint64_t{1} << k),
      currentIndex_(0),
      dna_(maxDimensions * idCount),
      currentPoint_(s)
{
    if (s == 0 || s > maxDimensions)
        throw std::out_of_range(
            "ART-Owen supports dimensions 1 to 1024");

    if (seed < static_cast<long>(std::numeric_limits<int>::min()) ||
        seed > static_cast<long>(std::numeric_limits<int>::max()))
    {
        throw std::out_of_range(
            "ART-Owen seed must fit in int");
    }

    const int artSeed =
        static_cast<int>(seed);

    pbrt::RNG rng(
        static_cast<std::uint64_t>(artSeed));

    for (std::uint32_t& value : dna_)
        value = rng.Uniform<std::uint32_t>();

    updateCurrentPoint();
}

std::uint32_t AdaptiveRegularTiles::scramble(
    std::uint32_t value,
    const std::uint32_t* dna)
{
    int id = 0;
    std::uint32_t result = value;

    for (int i = 0; i < 32; ++i)
    {
        result ^= dna[id] >> i;

        id = production32[id][
            (value >> (31 - i)) & 1u];
    }

    return result;
}

void AdaptiveRegularTiles::updateCurrentPoint()
{
    for (unsigned int j = 0; j < dimension_; ++j)
    {
        const std::uint32_t value =
            PbrtSobolCore::sobolInteger(
                currentIndex_,
                j);

        const std::uint32_t scrambled =
            scramble(
                value,
                dna_.data() + j * idCount);

        currentPoint_[j] =
            PbrtSobolCore::toUnitInterval(
                scrambled);
    }
}

double AdaptiveRegularTiles::coordinate(
    unsigned int j) const
{
    return currentPoint_[j];
}

void AdaptiveRegularTiles::next()
{
    if (currentIndex_ + 1 < numPoints_)
    {
        ++currentIndex_;
        updateCurrentPoint();
    }
}

std::uint64_t AdaptiveRegularTiles::numPoints() const
{
    return numPoints_;
}
