/*
 * Adapted from sscpp::hups64::CachedNUSPointSet, a port of
 * umontreal.ssj.hups64.CachedNUSPointSet.
 * Copyright (C) 2008-2021 The SSJ Authors, supervised by Pierre L'Ecuyer,
 * Université de Montréal.
 * Licensed under the Apache License, Version 2.0:
 * https://www.apache.org/licenses/LICENSE-2.0
 * Distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied. See the License for permissions and limits.
 */
#ifndef CACHED_NUS_POINT_SET_H
#define CACHED_NUS_POINT_SET_H

#include <cstddef>
#include <cstdint>
#include <vector>

// Full NUS using sscpp's sorted-order traversal. Construction prepares only;
// call randomize() before reading the initially zero coordinate cache.
class CachedNUSPointSet
{
public:
    struct DigitalNet
    {
        std::uint64_t numPoints;
        unsigned int dimension;
        unsigned int numColumns;
        unsigned int outputBits;
        // Right-aligned outputBits-wide columns[j * numColumns + c].
        std::vector<std::uint64_t> columns;
    };

    // Prepares a snapshot of the net, with 1..62 output bits. Matrices are
    // not retained. Point counts and dimensions must fit positive Java ints.
    CachedNUSPointSet(const DigitalNet& net, std::uint64_t seed);
    // Prepares 2^k PBRT Sobol points in Gray-code order with 32 output bits,
    // matching NewCachedNusPointSet's convenience constructor; requires k < 31.
    CachedNUSPointSet(unsigned int s, unsigned int k, long seed);

    CachedNUSPointSet(const CachedNUSPointSet&) = delete;
    CachedNUSPointSet& operator=(const CachedNUSPointSet&) = delete;

    // Reuses all storage, visiting each dimension in sorted coordinate order
    // and writing the scrambled doubles at their original point indices.
    // Consumes one random word per point/dimension, including duplicates.
    // Uses libseq's process-global drand48 stream: 1..32 bits use the high
    // bits of one 32-bit draw; 33..62 concatenate two draws, high draw first.
    void randomize();
    // Seeds the shared stream with the low 32 seed bits, as libseq does.
    // Preserves the current coordinates and iteration position.
    void reseed(std::uint64_t seed);

    // Same z * 2^(-outputBits) + 2^(-55) conversion as sscpp and SSJ.
    double getCoordinate(std::uint64_t i, unsigned int j) const;
    double coordinate(unsigned int j) const;
    void next();
    void reset();
    std::uint64_t numPoints() const;

private:
    static DigitalNet sobolNet(unsigned int s, unsigned int k);
    void prepare(const DigitalNet& net);
    std::uint64_t randomWord();

    std::size_t numPoints_ = 0;
    unsigned int dimension_ = 0;
    unsigned int outputBits_ = 0;
    double normFactor_ = 0.0;
    std::size_t currentIndex_ = 0;
    // originalIndices_[j * n + r] is the original index of sorted rank r.
    std::vector<std::uint32_t> originalIndices_;
    // Dimension-major masks, with numPoints - 1 masks per dimension.
    std::vector<std::uint64_t> firstDiffMasks_;
    std::vector<double> coordinates_;
};

#endif
