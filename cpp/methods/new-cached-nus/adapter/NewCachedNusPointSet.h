#ifndef NEW_CACHED_NUS_POINT_SET_H
#define NEW_CACHED_NUS_POINT_SET_H

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

// Full NUS using SSJ's NewCachedNusPointSet preparation and original-point
// traversal. Construction prepares a snapshot only: call randomize() before
// using the initially zero coordinate cache in an experiment.
class NewCachedNusPointSet
{
public:
    struct DigitalNet
    {
        std::uint64_t numPoints;
        unsigned int dimension;
        unsigned int numColumns;
        unsigned int outputBits;
        // Right-aligned outputBits-wide integers: columns[j * numColumns + c].
        // Zero columns are allowed for a one-point net. Extra columns are
        // allowed, but the vector must contain exactly dimension * numColumns.
        std::vector<std::uint64_t> columns;
    };

    // Requires positive Java signed-32-bit sizes, precision 1..62, and enough
    // columns for the point count. Input matrices are not retained.
    NewCachedNusPointSet(const DigitalNet& net, std::uint64_t seed);
    // Prepares 2^k PBRT Sobol points in Gray-code order, with 32 output bits.
    // Requires k < 31 and a supported PBRT dimension. Signed seeds convert
    // modulo 2^64; they are not restricted to the PBRT scramblers' int seeds.
    NewCachedNusPointSet(unsigned int s, unsigned int k, long seed);
    ~NewCachedNusPointSet();

    NewCachedNusPointSet(const NewCachedNusPointSet&) = delete;
    NewCachedNusPointSet& operator=(const NewCachedNusPointSet&) = delete;

    // Allocates nothing, consumes one word per point per dimension (including
    // duplicates), and resets iteration. For 1..32 bits, retain the high
    // outputBits of one PCG32 draw. For 33..62, concatenate the first draw as
    // the high 32 bits and the second as the low 32 bits, then retain the high
    // outputBits. There is no reservoir. Java uses another RNG: equal numeric
    // seeds do not imply equal output; parity requires equal matrices and words.
    void randomize();
    // Changes future randomness only; preserves the cache and iteration position.
    void reseed(std::uint64_t seed);

    // Java's z * 2^(-outputBits) + 2^(-55) conversion can round distinct
    // integers to the same double, including 1.0. Results are not clamped.
    double getCoordinate(std::uint64_t i, unsigned int j) const;
    double coordinate(unsigned int j) const;
    void next(); // Stays at the final point, as in the existing adapters.
    void reset();
    std::uint64_t numPoints() const;

private:
    // Allows correctness tests to inspect exact integers and prepared links
    // without making the upstream RNG or a test RNG part of the public API.
    friend struct NewCachedNusPointSetTestAccess;
    struct RandomState;
    static DigitalNet sobolNet(unsigned int s, unsigned int k);
    void prepare(const DigitalNet& net);
    std::uint64_t randomWord();

    std::size_t numPoints_ = 0;
    unsigned int dimension_ = 0;
    unsigned int outputBits_ = 0;
    double normFactor_ = 0.0;
    std::size_t currentIndex_ = 0;
    std::vector<std::uint32_t> parents_;
    // Dimension major, with numPoints - 1 masks per dimension, as in Java.
    std::vector<std::uint64_t> firstDiffMasks_;
    std::vector<double> coordinates_;
    std::vector<std::uint64_t> scrambledBits_;
    std::unique_ptr<RandomState> random_;
};

#endif
