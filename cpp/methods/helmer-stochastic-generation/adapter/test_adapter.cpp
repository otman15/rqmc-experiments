#include "HelmerStochasticSobol.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <new>
#include <stdexcept>
#include <vector>

#include "rng.h"

namespace
{
bool forbidAllocations = false;
}

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

struct HelmerStochasticSobolTestAccess
{
    static const auto& points(const HelmerStochasticSobol& p) { return p.points_; }
};

namespace
{
using Access = HelmerStochasticSobolTestAccess;
using Init = HelmerStochasticSobol::Initialization;

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

std::uint64_t bits(double value)
{
    std::uint64_t result;
    static_assert(sizeof(result) == sizeof(value), "binary64 required");
    std::memcpy(&result, &value, sizeof(result));
    return result;
}

std::uint64_t cacheHash(const HelmerStochasticSobol& p)
{
    std::uint64_t hash = 14695981039346656037ULL;
    for (const double value : Access::points(p))
    {
        hash ^= bits(value);
        hash *= 1099511628211ULL;
    }
    return hash;
}

void testLegacyFixtures()
{
    // Captured from the unmodified local adapter before the lifecycle refactor.
    // Hash every binary64 value in point-major order (xor then FNV multiply).
    struct Fixture { unsigned int s, k; long seed; std::uint64_t hash; };
    for (const auto& f : std::vector<Fixture>{
        {3, 3, 12345, 0x10c1a0d619442305ULL},
        {64, 5, -7, 0x9ba6e7269f86c325ULL},
        {64, 0, 0, 0x2e1b5272cb3d0825ULL},
        {1, 7, 20171215, 0x36fedcb67f8b6d25ULL}})
    {
        HelmerStochasticSobol p(f.s, f.k, f.seed);
        require(cacheHash(p) == f.hash, "legacy constructor seeded output changed");
        for (std::uint64_t i = 0; i < p.numPoints(); ++i)
        {
            for (unsigned int j = 0; j < f.s; ++j)
                require(bits(p.coordinate(j)) == bits(Access::points(p)[i * f.s + j]),
                        "point/dimension iteration order");
            p.next();
        }
        p.next();
        require(bits(p.coordinate(0)) == bits(Access::points(p)[(p.numPoints() - 1) * f.s]),
                "iteration should stay at the final point");
        p.reset();
        require(bits(p.coordinate(0)) == bits(Access::points(p)[0]), "reset position");
        require(cacheHash(p) == f.hash, "iteration/reset changed the cache");
    }
}

void testRandomization(unsigned int s, unsigned int k, long seed)
{
    HelmerStochasticSobol p(s, k, seed, Init::PrepareOnly);
    HelmerStochasticSobol legacy(s, k, seed);
    sampling::RNG rng(static_cast<std::uint64_t>(seed));
    const auto* storage = Access::points(p).data();
    const auto capacity = Access::points(p).capacity();
    for (const double value : Access::points(p))
        require(bits(value) == 0, "preparation-only constructor generated points");

    std::vector<double> firstCache;
    for (int rep = 0; rep < 4; ++rep)
    {
        // No allocations, even on the first randomization after preparation.
        struct AllocationBlock
        {
            AllocationBlock() { forbidAllocations = true; }
            ~AllocationBlock() { forbidAllocations = false; }
        };
        {
            AllocationBlock block;
            p.randomize();
        }
        if (rep != 0)
            legacy.randomize();
        require(Access::points(p) == Access::points(legacy), "constructor/lifecycle equivalence");
        require(Access::points(p).data() == storage && Access::points(p).capacity() == capacity,
                "randomization replaced cache storage");
        for (unsigned int j = 0; j < s; ++j)
        {
            require(bits(p.coordinate(j)) == bits(rng.GetUniformFloat()),
                    "RNG advancement must consume exactly n*s draws per scrambling");
            require(bits(p.coordinate(j)) == bits(Access::points(p)[j]), "randomize did not rewind");
        }
        for (std::uint64_t i = s; i < p.numPoints() * s; ++i)
            (void) rng.GetUniformInt();
        if (rep == 0)
            firstCache = Access::points(p);
        for (std::uint64_t i = 1; i < p.numPoints(); ++i)
            p.next();
    }
    const auto saved = Access::points(p);
    const auto current = bits(p.coordinate(0));
    p.reseed(static_cast<std::uint64_t>(seed));
    require(Access::points(p) == saved, "reseed changed cached points");
    require(bits(p.coordinate(0)) == current, "reseed changed iteration position");
    p.reset();
    require(Access::points(p) == saved, "reset changed cached points");
    p.randomize();
    require(Access::points(p) == firstCache, "reseed/reset did not replay the first scrambling");
    require(Access::points(p).data() == storage, "reseed replaced storage");

    // Preserve value-copy semantics, including independent future RNG state.
    HelmerStochasticSobol copy(p);
    copy.randomize();
    p.randomize();
    require(Access::points(copy) == Access::points(p), "copy lost independent RNG state");
    HelmerStochasticSobol assigned(1, 0, 0, Init::PrepareOnly);
    assigned = p;
    assigned.randomize();
    p.randomize();
    require(Access::points(assigned) == Access::points(p), "copy assignment lost state");
}

void testInvalidInputs()
{
    for (const auto s : {0u, 65u, std::numeric_limits<unsigned int>::max()})
    {
        bool rejected = false;
        try { HelmerStochasticSobol p(s, 0, 0, Init::PrepareOnly); }
        catch (const std::out_of_range&) { rejected = true; }
        require(rejected, "invalid dimension accepted");
    }
    for (const auto k : {31u, 32u, 64u, std::numeric_limits<unsigned int>::max()})
    {
        bool rejected = false;
        try { HelmerStochasticSobol p(1, k, 0, Init::PrepareOnly); }
        catch (const std::out_of_range&) { rejected = true; }
        require(rejected, "invalid exponent accepted before shift/allocation");
    }
}
}

int main()
{
    try
    {
        testLegacyFixtures();
        for (const auto s : {1u, 3u, 64u})
            for (const auto k : {0u, 1u, 5u})
                for (const long seed : {0L, 12345L, -7L})
                    testRandomization(s, k, seed);
        testInvalidInputs();
        std::cout << "HelmerStochasticSobol: all lifecycle checks passed\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "HelmerStochasticSobol test failed: " << error.what() << '\n';
        return 1;
    }
}
