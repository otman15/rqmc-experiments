#include "FriedelKellerLibseq.h"

#include <cstdlib>
#include <stdexcept>

#include "RandomizedTSSequence.h"

namespace
{
std::uint64_t checkedPointCount(unsigned int s, unsigned int k)
{
    if (s == 0 || k >= 31)
        throw std::invalid_argument("libseq requires s > 0 and k < 31");
    return std::uint64_t{1} << k;
}
}

FriedelKellerLibseq::FriedelKellerLibseq(
    unsigned int s,
    unsigned int k,
    long seed)
    : FriedelKellerLibseq(s, k, seed,
          static_cast<const std::vector<std::uint32_t>*>(nullptr))
{
}

FriedelKellerLibseq::FriedelKellerLibseq(
    unsigned int s,
    unsigned int k,
    long seed,
    const std::vector<std::uint32_t>& binaryColumns)
    : FriedelKellerLibseq(s, k, seed, &binaryColumns)
{
}

FriedelKellerLibseq::FriedelKellerLibseq(
    unsigned int s,
    unsigned int k,
    long seed,
    const std::vector<std::uint32_t>* binaryColumns)
    : numPoints_(checkedPointCount(s, k)),
      currentIndex_(0)
{
    if (binaryColumns && binaryColumns->size() != std::uint64_t{s} * k)
        throw std::invalid_argument("libseq generating matrix size must equal s * k");
    srand48(seed);

    ring_ = std::make_unique<Ring>(2);

    char generator[] = "Niederreiter";
    matrix_ = std::make_unique<C_matrix>(
        ring_.get(), generator, s);

    if (binaryColumns)
    {
        // C_matrix has no blank constructor. Overwrite every entry through its
        // public API before upstream generates or scrambles any points.
        for (unsigned int j = 0; j < s; ++j)
            for (unsigned int row = 0; row < 32; ++row)
                for (unsigned int c = 0; c < 32; ++c)
                    matrix_->set(j, row, c,
                        c < k ? ((*binaryColumns)[static_cast<std::size_t>(j) * k + c]
                                 >> (31 - row)) & 1u : 0u);
    }

    sequence_ = std::make_unique<RandomizedTSSequence_base_2>(
        matrix_.get(), s, numPoints_);
}

FriedelKellerLibseq::~FriedelKellerLibseq() = default;

void FriedelKellerLibseq::randomize()
{
    sequence_->random_restart();
    currentIndex_ = 0;
}

void FriedelKellerLibseq::reseed(std::uint64_t seed)
{
    srand48(static_cast<long>(static_cast<std::uint32_t>(seed)));
}

void FriedelKellerLibseq::reset()
{
    sequence_->restart();
    currentIndex_ = 0;
}

double FriedelKellerLibseq::coordinate(unsigned int j) const
{
    return (*sequence_)[j];
}

void FriedelKellerLibseq::next()
{
    if (currentIndex_ + 1 < numPoints_)
    {
        ++(*sequence_);
        ++currentIndex_;
    }
}

std::uint64_t FriedelKellerLibseq::numPoints() const
{
    return numPoints_;
}
