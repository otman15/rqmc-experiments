#include "HelmerStochasticSobol.h"

#include <stdexcept>
#include <limits>
#include <vector>

#include "rng.h"
#include "xor_values.h"

namespace
{
std::uint64_t checkedPointCount(unsigned int s, unsigned int k)
{
    // Validate before shifting or allocating, including preparation-only calls.
    if (s == 0 || s > static_cast<unsigned int>(sampling::MAX_SOBOL_DIM))
        throw std::out_of_range("Helmer stochastic Sobol supports dimensions 1 to 64");
    if (k >= 31)
        throw std::out_of_range("Helmer stochastic Sobol requires k < 31");
    const auto n = std::uint64_t{1} << k;
    if (n > std::numeric_limits<std::size_t>::max() / s ||
        n * s > std::vector<double>().max_size())
        throw std::length_error("Helmer stochastic Sobol cache is too large");
    return n;
}
}

struct HelmerStochasticSobol::RandomState
{
    explicit RandomState(std::uint64_t seed) : rng(seed) {}
    sampling::RNG rng;
};

HelmerStochasticSobol::HelmerStochasticSobol(
    unsigned int s,
    unsigned int k,
    long seed)
    : HelmerStochasticSobol(s, k, seed, Initialization::Randomize)
{
}

HelmerStochasticSobol::HelmerStochasticSobol(
    unsigned int s,
    unsigned int k,
    long seed,
    Initialization initialization)
    : dimension_(s),
      numPoints_(checkedPointCount(s, k)),
      currentIndex_(0),
      points_(numPoints_ * s),
      random_(std::make_unique<RandomState>(static_cast<std::uint64_t>(seed)))
{
    if (initialization == Initialization::Randomize)
        randomize();
}

HelmerStochasticSobol::~HelmerStochasticSobol() = default;

HelmerStochasticSobol::HelmerStochasticSobol(const HelmerStochasticSobol& other)
    : dimension_(other.dimension_),
      numPoints_(other.numPoints_),
      currentIndex_(other.currentIndex_),
      points_(other.points_),
      random_(std::make_unique<RandomState>(*other.random_))
{
}

HelmerStochasticSobol& HelmerStochasticSobol::operator=(const HelmerStochasticSobol& other)
{
    if (this != &other)
    {
        points_ = other.points_;
        *random_ = *other.random_;
        dimension_ = other.dimension_;
        numPoints_ = other.numPoints_;
        currentIndex_ = other.currentIndex_;
    }
    return *this;
}

void HelmerStochasticSobol::randomize()
{
    auto& rng = random_->rng;

    for (unsigned int d = 0; d < dimension_; ++d)
        points_[d] = rng.GetUniformFloat();

    for (unsigned int logN = 0;
         (std::uint64_t{1} << logN) < numPoints_;
         ++logN)
    {
        const std::uint64_t previousLength =
            std::uint64_t{1} << logN;

        const std::uint64_t numberOfStrata =
            previousLength * 2;

        const double inverseStrata =
            1.0 /
            static_cast<double>(numberOfStrata);

        for (std::uint64_t i = 0;
             i < previousLength &&
             previousLength + i < numPoints_;
             ++i)
        {
            const std::uint64_t newIndex =
                previousLength + i;

            for (unsigned int d = 0;
                 d < dimension_;
                 ++d)
            {
                const std::uint64_t previousIndex =
                    i ^
                    sampling::sobol_xors[d][logN];

                const int previousStratum =
                    static_cast<int>(
                        points_[
                            previousIndex * dimension_ + d]
                        * numberOfStrata);

                const int stratum =
                    previousStratum ^ 1;

                points_[
                    newIndex * dimension_ + d] =
                    (rng.GetUniformFloat() + stratum) *
                    inverseStrata;
            }
        }
    }
    reset();
}

void HelmerStochasticSobol::reseed(std::uint64_t seed)
{
    random_->rng = sampling::RNG(seed);
}

void HelmerStochasticSobol::reset()
{
    currentIndex_ = 0;
}

double HelmerStochasticSobol::coordinate(
    unsigned int j) const
{
    return points_[
        currentIndex_ * dimension_ + j];
}

void HelmerStochasticSobol::next()
{
    if (currentIndex_ + 1 < numPoints_)
        ++currentIndex_;
}

std::uint64_t HelmerStochasticSobol::numPoints() const
{
    return numPoints_;
}
