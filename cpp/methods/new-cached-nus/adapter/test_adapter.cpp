#include "NewCachedNusPointSet.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>
#include <new>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <pbrt/util/sobolmatrices.h>

#include "rng.h"

namespace
{
bool forbidAllocations = false;
}

// Check the no-allocation contract, including temporary allocations that would
// be invisible to vector capacity checks. The flag is scoped around hot calls.
void* operator new(std::size_t size)
{
    if (forbidAllocations)
        throw std::bad_alloc();
    if (void* memory = std::malloc(size == 0 ? 1 : size))
        return memory;
    throw std::bad_alloc();
}

void* operator new[](std::size_t size) { return ::operator new(size); }
void operator delete(void* memory) noexcept { std::free(memory); }
void operator delete[](void* memory) noexcept { std::free(memory); }
void operator delete(void* memory, std::size_t) noexcept { std::free(memory); }
void operator delete[](void* memory, std::size_t) noexcept { std::free(memory); }

struct NewCachedNusPointSetTestAccess
{
    static const auto& parents(const NewCachedNusPointSet& p) { return p.parents_; }
    static const auto& masks(const NewCachedNusPointSet& p) { return p.firstDiffMasks_; }
    static const auto& integers(const NewCachedNusPointSet& p) { return p.scrambledBits_; }
    static std::uint64_t word(NewCachedNusPointSet& p) { return p.randomWord(); }
};

namespace
{

using Net = NewCachedNusPointSet::DigitalNet;
using Access = NewCachedNusPointSetTestAccess;

// Exceptions, rather than assert(), keep every check active in Release builds.
void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

template <typename Exception, typename Action>
void expectThrows(Action action, const std::string& message)
{
    try
    {
        action();
    }
    catch (const Exception&)
    {
        return;
    }
    throw std::runtime_error("Expected exception: " + message);
}

std::uint64_t doubleBits(double value)
{
    static_assert(sizeof(double) == sizeof(std::uint64_t), "binary64 required");
    std::uint64_t bits;
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
}

double javaCoordinate(std::uint64_t z, unsigned int precision)
{
    // Force binary64 conversion before scaling and adding epsilon. In
    // particular, do not compute the expression in long double and round once.
    volatile double converted = static_cast<double>(static_cast<std::int64_t>(z));
    volatile double scaled = converted * std::ldexp(1.0, -static_cast<int>(precision));
    return scaled + 0x1p-55;
}

std::uint64_t referenceWord(sampling::RNG& rng, unsigned int precision)
{
    std::uint64_t value = rng.GetUniformInt();
    unsigned int width = 32;
    if (precision > 32)
    {
        value = value * 0x100000000ULL + rng.GetUniformInt();
        width = 64;
    }
    return value >> (width - precision);
}

std::vector<std::uint64_t> inputIntegers(const Net& net)
{
    std::vector<std::uint64_t> result(net.dimension * net.numPoints);
    for (unsigned int j = 0; j < net.dimension; ++j)
        for (std::uint64_t i = 0; i < net.numPoints; ++i)
        {
            // Direct matrix multiplication by Gray(i), independent of the
            // production trailing-zero recurrence.
            std::uint64_t gray = i ^ (i >> 1);
            std::uint64_t value = 0;
            for (unsigned int c = 0; gray != 0; ++c, gray >>= 1)
                if (gray & 1)
                    value ^= net.columns[static_cast<std::size_t>(j) * net.numColumns + c];
            result[static_cast<std::size_t>(j) * net.numPoints + i] = value;
        }
    return result;
}

std::vector<std::uint64_t> lazyTree(
    const Net& net,
    const std::vector<std::uint64_t>& input,
    sampling::RNG& rng)
{
    std::vector<std::uint64_t> output(input.size());
    for (unsigned int j = 0; j < net.dimension; ++j)
    {
        // Each (bit position, original binary prefix) owns one permutation of
        // {0,1}. The first visit chooses it using this point's random word.
        // Existing nodes reuse it, including the branching node. This reference
        // has no sorting, parent links, first-difference masks, or stack sweep.
        std::map<std::pair<unsigned int, std::uint64_t>, unsigned int> flips;
        for (std::uint64_t i = 0; i < net.numPoints; ++i)
        {
            const auto word = referenceWord(rng, net.outputBits);
            const auto offset = static_cast<std::size_t>(j) * net.numPoints + i;
            const auto x = input[offset];
            std::uint64_t z = 0;
            for (int bit = static_cast<int>(net.outputBits) - 1; bit >= 0; --bit)
            {
                const auto inputBit = static_cast<unsigned int>((x >> bit) & 1);
                const auto randomBit = static_cast<unsigned int>((word >> bit) & 1);
                const auto key = std::make_pair(static_cast<unsigned int>(bit), x >> (bit + 1));
                const auto node = flips.emplace(key, inputBit ^ randomBit).first;
                z |= std::uint64_t{inputBit ^ node->second} << bit;
            }
            output[offset] = z;
        }
    }
    return output;
}

std::vector<std::uint64_t> cacheBits(const NewCachedNusPointSet& p, unsigned int s)
{
    std::vector<std::uint64_t> result;
    for (unsigned int j = 0; j < s; ++j)
        for (std::uint64_t i = 0; i < p.numPoints(); ++i)
            result.push_back(doubleBits(p.getCoordinate(i, j)));
    return result;
}

Net makeNet(std::uint64_t n, unsigned int s, unsigned int precision, unsigned int pattern)
{
    unsigned int columns = 0;
    for (auto last = n - 1; last != 0; last >>= 1)
        ++columns;
    Net net{n, s, columns, precision, std::vector<std::uint64_t>(s * columns)};
    sampling::RNG rng(0x123456789abcdef0ULL + n + precision);
    for (unsigned int j = 0; j < s; ++j)
        for (unsigned int c = 0; c < columns; ++c)
        {
            auto value = referenceWord(rng, precision);
            if (pattern == 1)
                value = c % 2 == 0 ? 0 : value; // Singular matrices and duplicates.
            if (pattern == 2)
                value &= 7; // Deep shared prefixes even at 62-bit precision.
            net.columns[static_cast<std::size_t>(j) * columns + c] = value;
        }
    return net;
}

void checkAgainstTree(const Net& net, std::uint64_t seed, int repetitions = 3)
{
    NewCachedNusPointSet actual(net, seed);
    sampling::RNG rng(seed);
    const auto input = inputIntegers(net);
    const auto parents = Access::parents(actual);
    const auto masks = Access::masks(actual);
    require(actual.numPoints() == net.numPoints, "point count");
    for (const auto bits : cacheBits(actual, net.dimension))
        require(bits == 0, "construction must leave a zero cache");

    for (int rep = 0; rep < repetitions; ++rep)
    {
        const auto expected = lazyTree(net, input, rng);
        actual.next();
        actual.randomize();
        for (unsigned int j = 0; j < net.dimension; ++j)
            for (std::uint64_t i = 0; i < net.numPoints; ++i)
            {
                const auto offset = static_cast<std::size_t>(j) * net.numPoints + i;
                require(doubleBits(actual.getCoordinate(i, j)) ==
                        doubleBits(javaCoordinate(expected[offset], net.outputBits)),
                        "lazy permutation tree mismatch at precision " + std::to_string(net.outputBits));
                if (j + 1 == net.dimension)
                    require(Access::integers(actual)[i] == expected[offset],
                            "exact integer work array mismatch");
            }
        for (unsigned int j = 0; j < net.dimension; ++j)
            require(doubleBits(actual.coordinate(j)) == doubleBits(actual.getCoordinate(0, j)),
                    "randomize must reset iteration");
        require(Access::parents(actual) == parents, "randomization changed prepared parents");
        require(Access::masks(actual) == masks, "randomization changed prepared masks");
    }
}

void testReference()
{
    for (unsigned int precision = 1; precision <= 62; ++precision)
        checkAgainstTree(makeNet(7, 3, precision, precision % 3), 1234, 2);
    for (const unsigned int precision : {1u, 3u, 8u, 9u, 31u, 32u, 33u, 53u, 55u, 62u})
        for (const std::uint64_t n : {1u, 2u, 3u, 17u, 65u})
            for (unsigned int pattern = 0; pattern < 3; ++pattern)
                checkAgainstTree(makeNet(n, 4, precision, pattern), 0xfedcba9876543210ULL);

    // All duplicate coordinates, and n > 2^outputBits: every draw must still
    // be consumed, across both dimensions and subsequent randomizations.
    checkAgainstTree(Net{19, 3, 5, 1, std::vector<std::uint64_t>(15, 0)}, 0);
    // Extra, unused generating columns are legal and must not affect the points.
    checkAgainstTree(Net{3, 2, 4, 4, {8, 4, 15, 7, 3, 5, 12, 11}}, 42);
}

void testWords()
{
    for (unsigned int precision = 1; precision <= 62; ++precision)
        for (const std::uint64_t seed : {0ULL, 1ULL, 1234ULL, 0xffffffffffffffffULL})
        {
            NewCachedNusPointSet p(Net{1, 1, 0, precision, {}}, seed);
            sampling::RNG rng(seed);
            for (int i = 0; i < 16; ++i)
                require(Access::word(p) == referenceWord(rng, precision), "PCG word assembly");
        }
}

void testLifecycleAndSnapshot()
{
    auto net = makeNet(17, 3, 62, 0);
    NewCachedNusPointSet p(net, 1234), control(net, 1234);
    // No pointer to the input structure or its generating columns may survive.
    net.columns.assign(1, 0);
    net.numPoints = 0;
    net.dimension = 0;
    net.numColumns = 0;
    net.outputBits = 0;
    p.next();
    p.reset();
    p.reseed(1234);
    for (const auto bits : cacheBits(p, 3))
        require(bits == 0, "reset/reseed must preserve initial zeros");
    p.randomize();
    control.randomize();
    const auto first = cacheBits(p, 3);
    require(first == cacheBits(control, 3), "input snapshot independence");
    p.next();
    p.next();
    p.reseed(9876);
    require(cacheBits(p, 3) == first, "reseed changed cached coordinates");
    require(doubleBits(p.coordinate(1)) == doubleBits(p.getCoordinate(2, 1)),
            "reseed changed iteration position");
    p.reset();
    require(cacheBits(p, 3) == first, "reset changed cached coordinates");
    for (std::uint64_t i = 0; i < p.numPoints(); ++i)
    {
        for (unsigned int j = 0; j < 3; ++j)
            require(doubleBits(p.coordinate(j)) == doubleBits(p.getCoordinate(i, j)),
                    "iteration order");
        p.next();
    }
    p.next();
    require(doubleBits(p.coordinate(0)) == doubleBits(p.getCoordinate(p.numPoints() - 1, 0)),
            "next must stay at the final point");
    control.reseed(9876);
    for (int rep = 0; rep < 3; ++rep)
    {
        p.reset();
        p.randomize();
        control.randomize();
        require(cacheBits(p, 3) == cacheBits(control, 3), "reseed/reset RNG continuity");
        require(doubleBits(p.coordinate(0)) == doubleBits(p.getCoordinate(0, 0)),
                "randomize must rewind iteration");
    }
    p.reseed(1234);
    p.randomize();
    require(cacheBits(p, 3) == first, "reseed must replay the original scrambling");
    expectThrows<std::out_of_range>([&] { p.getCoordinate(17, 0); }, "point bounds");
    expectThrows<std::out_of_range>([&] { p.getCoordinate(0, 3); }, "dimension bounds");
    expectThrows<std::out_of_range>([&] { p.coordinate(3); }, "current coordinate bounds");
    expectThrows<std::out_of_range>([&] {
        p.getCoordinate(std::numeric_limits<std::uint64_t>::max(), 0);
    }, "large point index");
}

void testInvalidInputs()
{
    const Net valid{3, 2, 2, 32, {0x80000000, 0x40000000, 0x80000000, 0xc0000000}};
    const auto reject = [&](auto change, const std::string& name)
    {
        auto net = valid;
        change(net);
        expectThrows<std::invalid_argument>([&] { NewCachedNusPointSet p(net, 0); }, name);
    };
    reject([](Net& n) { n.numPoints = 0; }, "zero points");
    reject([](Net& n) { n.numPoints = std::uint64_t{1} << 31; }, "Java point limit");
    reject([](Net& n) { n.numPoints = std::numeric_limits<std::uint64_t>::max(); }, "large count");
    reject([](Net& n) { n.dimension = 0; }, "zero dimensions");
    reject([](Net& n) { n.dimension = std::numeric_limits<unsigned int>::max(); }, "large dimension");
    for (const unsigned int bits : {0u, 63u, 64u, std::numeric_limits<unsigned int>::max()})
        reject([&](Net& n) { n.outputBits = bits; }, "invalid precision");
    for (const unsigned int columns : {0u, 1u, std::numeric_limits<unsigned int>::max()})
        reject([&](Net& n) { n.numColumns = columns; }, "invalid column count");
    reject([](Net& n) { n.columns.pop_back(); }, "missing matrix value");
    reject([](Net& n) { n.columns.push_back(0); }, "extra matrix value");
    reject([](Net& n) { n.columns[0] = std::uint64_t{1} << 32; }, "wide 32-bit column");
    reject([](Net& n) { n.columns[3] = std::numeric_limits<std::uint64_t>::max(); }, "wide column");
    expectThrows<std::invalid_argument>([] {
        NewCachedNusPointSet p(Net{1, 1, 1, 62, {std::uint64_t{1} << 62}}, 0);
    }, "invalid unused column");
    expectThrows<std::length_error>([] {
        NewCachedNusPointSet p(Net{2147483647, 2147483647u, 31, 62, {}}, 0);
    }, "unsafe flat storage sizes must fail before allocation");

    for (const unsigned int k : {31u, 32u, 63u, 64u, std::numeric_limits<unsigned int>::max()})
        expectThrows<std::invalid_argument>([&] { NewCachedNusPointSet p(1, k, 0); }, "unsafe shift");
    expectThrows<std::invalid_argument>([] { NewCachedNusPointSet p(0, 0, 0); }, "zero Sobol dimension");
    expectThrows<std::out_of_range>([] {
        NewCachedNusPointSet p(static_cast<unsigned int>(pbrt::NSobolDimensions) + 1, 0, 0);
    }, "PBRT dimension limit, including k=0");
    expectThrows<std::out_of_range>([] {
        NewCachedNusPointSet p(std::numeric_limits<unsigned int>::max(), 1, 0);
    }, "very large Sobol dimension");
    checkAgainstTree(Net{1, 1, 1, 62, {(std::uint64_t{1} << 62) - 1}}, 0);
}

void testSobol()
{
    for (const unsigned int k : {0u, 1u, 3u, 7u})
        for (const long seed : {0L, -1L, 1234L, std::numeric_limits<long>::min(),
                                std::numeric_limits<long>::max()})
        {
            Net net{std::uint64_t{1} << k, 4, k, 32, std::vector<std::uint64_t>(4 * k)};
            // Extract the original PBRT table independently of the core calls
            // used by the convenience constructor; do not use double samples.
            for (unsigned int j = 0; j < net.dimension; ++j)
                for (unsigned int c = 0; c < k; ++c)
                    net.columns[j * k + c] = pbrt::SobolMatrices32[j * pbrt::SobolMatrixSize + c];
            NewCachedNusPointSet supplied(net, static_cast<std::uint64_t>(seed));
            NewCachedNusPointSet convenience(4, k, seed);
            require(cacheBits(supplied, 4) == cacheBits(convenience, 4), "initial Sobol zeros");
            require(Access::parents(supplied) == Access::parents(convenience), "Sobol parents");
            require(Access::masks(supplied) == Access::masks(convenience), "Sobol masks");
            for (int rep = 0; rep < 3; ++rep)
            {
                supplied.randomize();
                convenience.randomize();
                require(cacheBits(supplied, 4) == cacheBits(convenience, 4), "PBRT matrix equivalence");
            }
            checkAgainstTree(net, static_cast<std::uint64_t>(seed), 1);
        }
    NewCachedNusPointSet boundary(pbrt::NSobolDimensions, 0, 0);
    boundary.randomize();
    require(boundary.numPoints() == 1, "maximum supported PBRT dimension");
}

void testJavaParity()
{
    // Golden data obtained by running the unchanged Java NewCachedNusPointSet
    // DigitalNetBase2 constructor with these matrices. A RandomStream adapter
    // supplies PCG32(seed=1234, sequence=0) words using the documented assembly.
    // These fixtures test actual Java results, not a C++ copy of its parent code.
    const auto high = std::uint64_t{1} << 61;
    const Net net{7, 3, 3, 62, {3, 4, 1, 0, 2, 2, high | 3, (high >> 1) | 1, (high >> 2) | 5}};
    NewCachedNusPointSet p(net, 1234);
    const auto root = std::numeric_limits<std::uint32_t>::max();
    const std::vector<std::uint32_t> parents{
        root, 0, 1, 2, 3, 2, 1, root, 0, 1, 2, 1, 4, 3, root, 0, 1, 0, 3, 2, 1};
    const std::vector<std::uint64_t> masks{
        3, 7, 3, 1, 1, 1, 0, 3, 0, 0, 0, 0,
        4611686018427387903ULL, 2305843009213693951ULL, 2305843009213693951ULL,
        1152921504606846975ULL, 1152921504606846975ULL, 1152921504606846975ULL};
    require(Access::parents(p) == parents, "Java prepared parent parity");
    require(Access::masks(p) == masks, "Java first-difference mask parity");

    // The first two dimensions have distinct integers sharing long prefixes;
    // Java rounds all seven coordinates to the same double in each dimension.
    const std::array<std::uint64_t, 2> firstDimension{0x3fd6b9e4bebf7ee6ULL, 0x3fd7ab06c94c949aULL};
    const std::array<std::uint64_t, 2> secondDimension{0x3fd5004b80edabe6ULL, 0x3fb6591cf21daf7bULL};
    const std::array<std::uint64_t, 14> lastDimension{
        0x3fe0235de7e4b11eULL, 0x3fb452ffa96e9292ULL, 0x3fd8d78d0f4a8778ULL,
        0x3fea23d6986c9b1fULL, 0x3fedb888f60bf681ULL, 0x3fd39d13b768a052ULL, 0x3fc4ee94dadce851ULL,
        0x3fd127c5a1cd2d26ULL, 0x3fee3deab1ad291fULL, 0x3fe49571e44977bfULL,
        0x3fcdf9ba459bf902ULL, 0x3f972033764183cbULL, 0x3fe1e96c008db929ULL, 0x3fe82b54cc25fc83ULL};
    const std::array<std::uint64_t, 14> integers{
        2325752758694722575ULL, 366128483004621853ULL, 1790054473717479162ULL,
        3767170038632955584ULL, 4283224426317415113ULL, 1413307536519483938ULL, 754153296472844306ULL,
        1236173921459184978ULL, 4358311687104904606ULL, 2966433755175026057ULL,
        1079981152216580200ULL, 104149277820337198ULL, 2581363035703759592ULL, 3483157833392391853ULL};
    for (std::size_t rep = 0; rep < 2; ++rep)
    {
        p.randomize();
        for (std::size_t i = 0; i < 7; ++i)
        {
            require(doubleBits(p.getCoordinate(i, 0)) == firstDimension[rep], "Java dimension 0 parity");
            require(doubleBits(p.getCoordinate(i, 1)) == secondDimension[rep], "Java duplicate parity");
            require(doubleBits(p.getCoordinate(i, 2)) == lastDimension[rep * 7 + i], "Java binary64 parity");
            require(Access::integers(p)[i] == integers[rep * 7 + i], "Java exact integer parity");
        }
    }
}

void testNoAllocations()
{
    NewCachedNusPointSet p(makeNet(65, 4, 62, 0), 1234);
    struct AllocationBlock
    {
        AllocationBlock() { forbidAllocations = true; }
        ~AllocationBlock() { forbidAllocations = false; }
    };
    {
        AllocationBlock block;
        for (int rep = 0; rep < 10; ++rep)
        {
            p.randomize();
            p.next();
            p.reset();
            p.reseed(static_cast<std::uint64_t>(rep));
        }
    }
}

}

int main()
{
    try
    {
        testWords();
        testReference();
        testLifecycleAndSnapshot();
        testInvalidInputs();
        testSobol();
        testJavaParity();
        testNoAllocations();
        std::cout << "NewCachedNusPointSet: all correctness checks passed\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "NewCachedNusPointSet test failed: " << error.what() << '\n';
        return 1;
    }
}
