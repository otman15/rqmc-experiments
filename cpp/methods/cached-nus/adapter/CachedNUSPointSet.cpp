/*
 * Adapted from sscpp/src/hups64/cached_nus_point_set.cpp, a port of
 * umontreal.ssj.hups64.CachedNUSPointSet.
 * Copyright (C) 2008-2021 The SSJ Authors, supervised by Pierre L'Ecuyer,
 * Université de Montréal.
 * Licensed under the Apache License, Version 2.0:
 * https://www.apache.org/licenses/LICENSE-2.0
 * Distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied. See the License for permissions and limits.
 *
 * The sorted-order preparation and scrambling are retained. This adapter uses
 * the project's C++17 sequence interface, packed digital nets, and libseq's
 * native drand48 stream instead of sscpp's PointSet/RandomStream interfaces.
 */
#include "CachedNUSPointSet.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <utility>

#include "PbrtSobolCore.h"

namespace
{
constexpr std::uint64_t javaIntMax = 2147483647;

std::size_t checkedProduct(std::size_t a, std::size_t b)
{
    if (b != 0 && a > std::numeric_limits<std::size_t>::max() / b)
        throw std::length_error("CachedNUSPointSet size product overflows");
    return a * b;
}

template <typename T>
void checkArraySize(std::size_t count)
{
    if (count > std::vector<T>().max_size())
        throw std::length_error("CachedNUSPointSet array is too large");
}

unsigned int trailingZeros(std::uint32_t value)
{
    // C++17 equivalent of sscpp's std::countr_zero, for positive indices.
    unsigned int count = 0;
    while ((value & 1u) == 0)
    {
        value >>= 1;
        ++count;
    }
    return count;
}

std::uint64_t highestOneBit(std::uint64_t value)
{
    // C++17 equivalent of sscpp's std::bit_floor (also returns 0 for 0).
    value |= value >> 1;
    value |= value >> 2;
    value |= value >> 4;
    value |= value >> 8;
    value |= value >> 16;
    value |= value >> 32;
    return value - (value >> 1);
}
}

CachedNUSPointSet::CachedNUSPointSet(const DigitalNet& net, std::uint64_t seed)
{
    static_assert(std::numeric_limits<double>::is_iec559 &&
                  std::numeric_limits<double>::digits == 53,
                  "Java parity requires IEEE binary64 doubles");

    if (net.numPoints == 0 || net.numPoints > javaIntMax ||
        net.dimension == 0 || net.dimension > javaIntMax)
        throw std::invalid_argument("Point count and dimension must be positive Java ints");
    if (net.outputBits < 1 || net.outputBits > 62)
        throw std::invalid_argument("outputBits must be between 1 and 62");

    unsigned int requiredColumns = 0;
    for (std::uint64_t last = net.numPoints - 1; last != 0; last >>= 1)
        ++requiredColumns;
    if (net.numColumns < requiredColumns || net.numColumns > javaIntMax)
        throw std::invalid_argument("Invalid generating column count");
    if (net.numPoints > std::numeric_limits<std::size_t>::max() ||
        net.dimension > std::numeric_limits<std::size_t>::max() ||
        net.numColumns > std::numeric_limits<std::size_t>::max())
        throw std::length_error("Digital net sizes do not fit size_t");

    const auto n = static_cast<std::size_t>(net.numPoints);
    const auto s = static_cast<std::size_t>(net.dimension);
    const auto columnCount = checkedProduct(s, net.numColumns);
    const auto coordinateCount = checkedProduct(s, n);
    const auto maskCount = checkedProduct(s, n - 1);
    checkArraySize<std::uint64_t>(columnCount);
    checkArraySize<std::uint32_t>(coordinateCount);
    checkArraySize<double>(coordinateCount);
    checkArraySize<std::uint64_t>(maskCount);
    checkArraySize<std::uint64_t>(n);
    checkArraySize<std::uint32_t>(n);
    if (net.columns.size() != columnCount)
        throw std::invalid_argument("Generating matrix size does not match the net");
    const auto limit = std::uint64_t{1} << net.outputBits;
    for (const auto column : net.columns)
        if (column >= limit)
            throw std::invalid_argument("Generating column exceeds outputBits");

    numPoints_ = n;
    dimension_ = net.dimension;
    outputBits_ = net.outputBits;
    normFactor_ = std::ldexp(1.0, -static_cast<int>(outputBits_));
    originalIndices_.resize(coordinateCount);
    firstDiffMasks_.resize(maskCount);
    coordinates_.resize(coordinateCount, 0.0);
    reseed(seed);
    prepare(net);
}

CachedNUSPointSet::DigitalNet CachedNUSPointSet::sobolNet(unsigned int s, unsigned int k)
{
    if (k >= 31)
        throw std::invalid_argument("Sobol point count requires k < 31");
    if (s == 0)
        throw std::invalid_argument("Sobol dimension must be positive");
    (void) PbrtSobolCore::sobolInteger(0, s - 1);
    const auto columnCount = checkedProduct(s, k);
    checkArraySize<std::uint64_t>(columnCount);
    DigitalNet net{std::uint64_t{1} << k, s, k, 32,
                   std::vector<std::uint64_t>(columnCount)};
    for (unsigned int j = 0; j < s; ++j)
        for (unsigned int c = 0; c < k; ++c)
            net.columns[static_cast<std::size_t>(j) * k + c] =
                PbrtSobolCore::sobolInteger(std::uint64_t{1} << c, j);
    return net;
}

CachedNUSPointSet::CachedNUSPointSet(unsigned int s, unsigned int k, long seed)
    : CachedNUSPointSet(sobolNet(s, k), static_cast<std::uint64_t>(seed))
{
}

void CachedNUSPointSet::prepare(const DigitalNet& net)
{
    std::vector<std::uint64_t> bits(numPoints_), sortedBits(numPoints_);
    std::vector<std::uint32_t> tempIndices(numPoints_);
    std::array<std::uint32_t, 256> counts{}, positions{};
    const unsigned int numPasses = (outputBits_ + 7) / 8;

    for (unsigned int j = 0; j < dimension_; ++j)
    {
        auto* dimensionIndices = originalIndices_.data() +
            static_cast<std::size_t>(j) * numPoints_;
        // Choose the starting buffer so the final pass ends in originalIndices_.
        auto* indices = (numPasses & 1u) == 0 ? dimensionIndices : tempIndices.data();
        auto* sortedIndices = (numPasses & 1u) == 0 ? tempIndices.data() : dimensionIndices;
        const auto columnOffset = static_cast<std::size_t>(j) * net.numColumns;
        bits[0] = 0;
        for (std::size_t i = 1; i < numPoints_; ++i)
            bits[i] = bits[i - 1] ^
                net.columns[columnOffset + trailingZeros(static_cast<std::uint32_t>(i))];
        for (std::size_t i = 0; i < numPoints_; ++i)
            indices[i] = static_cast<std::uint32_t>(i);

        // Stable base-256 radix sorting preserves original order for ties.
        for (unsigned int shift = 0; shift < outputBits_; shift += 8)
        {
            counts.fill(0);
            for (const auto value : bits)
                ++counts[(value >> shift) & 0xffu];
            positions[0] = 0;
            for (std::size_t bin = 1; bin < counts.size(); ++bin)
                positions[bin] = positions[bin - 1] + counts[bin - 1];
            for (std::size_t i = 0; i < numPoints_; ++i)
            {
                const auto bin = (bits[i] >> shift) & 0xffu;
                const auto destination = positions[bin]++;
                sortedBits[destination] = bits[i];
                sortedIndices[destination] = indices[i];
            }
            bits.swap(sortedBits);
            std::swap(indices, sortedIndices);
        }

        const auto maskOffset = static_cast<std::size_t>(j) * (numPoints_ - 1);
        for (std::size_t i = 1; i < numPoints_; ++i)
        {
            const auto firstDifference = highestOneBit(bits[i - 1] ^ bits[i]);
            firstDiffMasks_[maskOffset + i - 1] =
                firstDifference == 0 ? 0 : (firstDifference << 1) - 1;
        }
    }
}

std::uint64_t CachedNUSPointSet::randomWord()
{
    const std::uint64_t first = static_cast<std::uint32_t>(drand48() * 0x1p32);
    if (outputBits_ <= 32)
        return first >> (32 - outputBits_);
    const std::uint64_t second = static_cast<std::uint32_t>(drand48() * 0x1p32);
    return ((first << 32) | second) >> (64 - outputBits_);
}

void CachedNUSPointSet::randomize()
{
    for (unsigned int j = 0; j < dimension_; ++j)
    {
        const auto offset = static_cast<std::size_t>(j) * numPoints_;
        const auto maskOffset = static_cast<std::size_t>(j) * (numPoints_ - 1);
        const auto* indices = originalIndices_.data() + offset;
        auto* output = coordinates_.data() + offset;
        auto z = randomWord();
        output[indices[0]] = static_cast<double>(z) * normFactor_ + 0x1p-55;
        for (std::size_t i = 1; i < numPoints_; ++i)
        {
            const auto word = randomWord();
            const auto mask = firstDiffMasks_[maskOffset + i - 1];
            const auto lowerMask = mask >> 1;
            const auto firstDiffBit = mask ^ lowerMask;
            // Duplicates keep z unchanged, but still consume a random word.
            z = ((z ^ firstDiffBit) & ~lowerMask) | (word & lowerMask);
            output[indices[i]] = static_cast<double>(z) * normFactor_ + 0x1p-55;
        }
    }
    reset();
}

void CachedNUSPointSet::reseed(std::uint64_t seed)
{
    srand48(static_cast<long>(static_cast<std::uint32_t>(seed)));
}

double CachedNUSPointSet::getCoordinate(std::uint64_t i, unsigned int j) const
{
    if (i >= numPoints_ || j >= dimension_)
        throw std::out_of_range("CachedNUSPointSet coordinate index out of range");
    return coordinates_[static_cast<std::size_t>(j) * numPoints_ +
                        static_cast<std::size_t>(i)];
}

double CachedNUSPointSet::coordinate(unsigned int j) const
{
    return getCoordinate(currentIndex_, j);
}

void CachedNUSPointSet::next()
{
    if (currentIndex_ + 1 < numPoints_)
        ++currentIndex_;
}

void CachedNUSPointSet::reset()
{
    currentIndex_ = 0;
}

std::uint64_t CachedNUSPointSet::numPoints() const
{
    return numPoints_;
}
