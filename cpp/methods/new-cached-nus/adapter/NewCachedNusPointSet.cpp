/*
 * Port of umontreal.ssj.hups64.NewCachedNusPointSet, licensed under the
 * Apache License, Version 2.0: https://www.apache.org/licenses/LICENSE-2.0
 * Distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied. See the License for permissions and limits.
 */
#include "NewCachedNusPointSet.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <utility>

#include "PbrtSobolCore.h"

namespace
{

constexpr std::uint32_t noParent = std::numeric_limits<std::uint32_t>::max();
constexpr std::uint64_t javaIntMax = 2147483647;

std::size_t checkedProduct(std::size_t a, std::size_t b)
{
    if (b != 0 && a > std::numeric_limits<std::size_t>::max() / b)
        throw std::length_error("NewCachedNusPointSet size product overflows");
    return a * b;
}

std::size_t checkedSum(std::size_t a, std::size_t b)
{
    if (a > std::numeric_limits<std::size_t>::max() - b)
        throw std::length_error("NewCachedNusPointSet storage size overflows");
    return a + b;
}

template <typename T>
void checkArraySize(std::size_t count)
{
    if (count > std::vector<T>().max_size())
        throw std::length_error("NewCachedNusPointSet array is too large");
}

unsigned int trailingZeros(std::uint32_t value)
{
    // Called only for positive point indices.
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
    value |= value >> 1;
    value |= value >> 2;
    value |= value >> 4;
    value |= value >> 8;
    value |= value >> 16;
    value |= value >> 32;
    return value - (value >> 1);
}

}

struct NewCachedNusPointSet::RandomState
{
    // Match libseq's process-global drand48 stream and low-32-bit seeds.
    explicit RandomState(std::uint64_t seed) { reseed(seed); }

    void reseed(std::uint64_t seed)
    {
        srand48(static_cast<long>(static_cast<std::uint32_t>(seed)));
    }

    std::uint32_t nextWord()
    {
        // Same extraction as libseq's get_32_random_bits().
        return static_cast<std::uint32_t>(drand48() * 0x1p32);
    }
};

NewCachedNusPointSet::NewCachedNusPointSet(
    const DigitalNet& net,
    std::uint64_t seed)
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

    // Check peak storage arithmetic as well as every individual vector before
    // allocating: permanent arrays, work array, and five preparation buffers.
    const auto cacheBytes = checkedProduct(coordinateCount,
        sizeof(std::uint32_t) + sizeof(double));
    const auto maskBytes = checkedProduct(maskCount, sizeof(std::uint64_t));
    const auto workAndTemporaryBytes = checkedProduct(n,
        3 * sizeof(std::uint64_t) + 3 * sizeof(std::uint32_t));
    (void) checkedSum(checkedSum(cacheBytes, maskBytes), workAndTemporaryBytes);

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
    parents_.resize(coordinateCount);
    firstDiffMasks_.resize(maskCount);
    coordinates_.resize(coordinateCount, 0.0);
    scrambledBits_.resize(n);
    random_ = std::make_unique<RandomState>(seed);
    prepare(net);
}

NewCachedNusPointSet::DigitalNet NewCachedNusPointSet::sobolNet(
    unsigned int s,
    unsigned int k)
{
    if (k >= 31)
        throw std::invalid_argument("Sobol point count requires k < 31");
    if (s == 0)
        throw std::invalid_argument("Sobol dimension must be positive");
    // Validate even when k == 0, before allocation. The core owns PBRT's limits.
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

NewCachedNusPointSet::NewCachedNusPointSet(
    unsigned int s,
    unsigned int k,
    long seed)
    : NewCachedNusPointSet(sobolNet(s, k), static_cast<std::uint64_t>(seed))
{
}

NewCachedNusPointSet::~NewCachedNusPointSet() = default;

void NewCachedNusPointSet::prepare(const DigitalNet& net)
{
    std::vector<std::uint64_t> bits(numPoints_), sortedBits(numPoints_);
    std::vector<std::uint32_t> indices(numPoints_), sortedIndices(numPoints_);
    std::vector<std::uint32_t> stack(numPoints_);
    std::array<std::uint32_t, 256> counts{}, positions{};

    for (unsigned int j = 0; j < dimension_; ++j)
    {
        bits[0] = 0;
        const auto columnOffset = static_cast<std::size_t>(j) * net.numColumns;
        for (std::size_t i = 1; i < numPoints_; ++i)
            bits[i] = bits[i - 1] ^
                net.columns[columnOffset + trailingZeros(static_cast<std::uint32_t>(i))];
        for (std::size_t i = 0; i < numPoints_; ++i)
            indices[i] = static_cast<std::uint32_t>(i);

        // Stable base-256 radix sort; ties keep their original point order.
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
            indices.swap(sortedIndices);
        }

        const auto parentOffset = static_cast<std::size_t>(j) * numPoints_;
        const auto maskOffset = static_cast<std::size_t>(j) * (numPoints_ - 1);
        std::size_t stackSize = 0;
        // Port of the active Java stack construction. The two candidates are
        // the nearest sorted neighbors with earlier original indices. Keep the
        // strict comparisons, including the left candidate on equal XORs.
        for (std::size_t next = 0; next <= numPoints_; ++next)
        {
            while (stackSize > 0 && (next == numPoints_ ||
                   indices[stack[stackSize - 1]] > indices[next]))
            {
                const auto rank = stack[--stackSize];
                const auto i = indices[rank];
                if (i == 0)
                    continue;
                auto parent = noParent;
                auto difference = static_cast<std::uint64_t>(
                    std::numeric_limits<std::int64_t>::max());
                if (stackSize != 0)
                {
                    const auto previous = stack[stackSize - 1];
                    parent = indices[previous];
                    difference = bits[rank] ^ bits[previous];
                }
                if (next < numPoints_)
                {
                    const auto nextDifference = bits[rank] ^ bits[next];
                    if (nextDifference < difference)
                    {
                        parent = indices[next];
                        difference = nextDifference;
                    }
                }
                parents_[parentOffset + i] = parent;
                const auto firstDifference = highestOneBit(difference);
                firstDiffMasks_[maskOffset + i - 1] =
                    firstDifference == 0 ? 0 : (firstDifference << 1) - 1;
            }
            if (next < numPoints_)
                stack[stackSize++] = static_cast<std::uint32_t>(next);
        }
        parents_[parentOffset] = noParent;
    }
    // All sorting buffers and the stack are released here; no matrices remain.
}

std::uint64_t NewCachedNusPointSet::randomWord()
{
    const std::uint64_t first = random_->nextWord();
    if (outputBits_ <= 32)
        return first >> (32 - outputBits_);
    const std::uint64_t second = random_->nextWord();
    return ((first << 32) | second) >> (64 - outputBits_);
}

void NewCachedNusPointSet::randomize()
{
    for (unsigned int j = 0; j < dimension_; ++j)
    {
        const auto offset = static_cast<std::size_t>(j) * numPoints_;
        const auto maskOffset = static_cast<std::size_t>(j) * (numPoints_ - 1);
        auto z = randomWord();
        scrambledBits_[0] = z;
        coordinates_[offset] = static_cast<double>(z) * normFactor_ + 0x1p-55;
        for (std::size_t i = 1; i < numPoints_; ++i)
        {
            const auto word = randomWord(); // Also consumed for duplicates.
            const auto mask = firstDiffMasks_[maskOffset + i - 1];
            const auto lowerMask = mask >> 1;
            const auto firstDiffBit = mask ^ lowerMask;
            const auto parent = scrambledBits_[parents_[offset + i]];
            z = ((parent ^ firstDiffBit) & ~lowerMask) | (word & lowerMask);
            scrambledBits_[i] = z;
            // Convert the integer to binary64 first, as Java does. Keeping
            // exact integers above avoids recovering rounded parents from doubles.
            coordinates_[offset + i] = static_cast<double>(z) * normFactor_ + 0x1p-55;
        }
    }
    reset();
}

void NewCachedNusPointSet::reseed(std::uint64_t seed)
{
    random_->reseed(seed);
}

double NewCachedNusPointSet::getCoordinate(std::uint64_t i, unsigned int j) const
{
    if (i >= numPoints_ || j >= dimension_)
        throw std::out_of_range("NewCachedNusPointSet coordinate index out of range");
    return coordinates_[static_cast<std::size_t>(j) * numPoints_ +
                        static_cast<std::size_t>(i)];
}

double NewCachedNusPointSet::coordinate(unsigned int j) const
{
    return getCoordinate(currentIndex_, j);
}

void NewCachedNusPointSet::next()
{
    if (currentIndex_ + 1 < numPoints_)
        ++currentIndex_;
}

void NewCachedNusPointSet::reset()
{
    currentIndex_ = 0;
}

std::uint64_t NewCachedNusPointSet::numPoints() const
{
    return numPoints_;
}
