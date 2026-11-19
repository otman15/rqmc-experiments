/**
 * @file random_stream.hpp
 * @brief The `RandomStream` interface: a stream of uniform random numbers.
 *
 * C++ counterpart of `umontreal.ssj.rng.RandomStream` (SSJ).
 */
#pragma once

#include <cstdint>
#include <string>

/**
 * @namespace sscpp::rng
 * @brief Uniform random number generators (counterpart of `umontreal.ssj.rng`).
 */
namespace sscpp::rng {

/**
 * @brief Interface for a stream of uniform (pseudo)random numbers.
 *
 * Counterpart of the Java interface `umontreal.ssj.rng.RandomStream`.
 *
 * The period of a generator is cut into adjacent *streams* (segments), and
 * each stream is cut into *substreams*. For a stream @f$g@f$, write
 * @f$I_g@f$ for its initial state, @f$B_g@f$ for the state at the beginning
 * of the current substream, and @f$C_g@f$ for the current state. The reset
 * methods move @f$C_g@f$ (and @f$B_g@f$) between these points.
 *
 * Only the methods needed by the NUS experiments are ported for now. Missing
 * compared with Java: `nextDoubleNonzero`, `nextArrayOfDouble`,
 * `nextArrayOfInt` and `nextArrayOfLong`. As in Java, this interface does not
 * declare `clone()`; each generator class provides its own.
 *
 * Java `int` is mapped to `int` and Java `long` to `std::int64_t`.
 */
class RandomStream {
public:
   virtual ~RandomStream() = default;

   /**
    * @brief Reinitializes the stream to its initial state @f$I_g@f$:
    *        @f$C_g@f$ and @f$B_g@f$ are set to @f$I_g@f$.
    */
   virtual void reset_start_stream() = 0;

   /**
    * @brief Reinitializes the stream to the beginning of its current
    *        substream: @f$C_g@f$ is set to @f$B_g@f$.
    */
   virtual void reset_start_substream() = 0;

   /**
    * @brief Reinitializes the stream to the beginning of its next substream:
    *        @f$N_g@f$ is computed, and @f$C_g@f$ and @f$B_g@f$ are set to
    *        @f$N_g@f$.
    */
   virtual void reset_next_substream() = 0;

   /**
    * @brief Returns a string containing the current state of this stream.
    */
   virtual std::string to_string() const = 0;

   /**
    * @brief Returns a uniform random number over @f$(0,1)@f$ and advances the
    *        state by one step.
    */
   virtual double next_double() = 0;

   /**
    * @brief Returns a random integer uniformly distributed over
    *        @f$\{i, i+1, \dots, j\}@f$.
    * @param i smallest integer that can be generated
    * @param j greatest integer that can be generated
    */
   virtual int next_int(int i, int j) = 0;

   /**
    * @brief Returns a random 64-bit integer uniformly distributed over
    *        @f$\{i, i+1, \dots, j\}@f$.
    * @param i smallest integer that can be generated
    * @param j greatest integer that can be generated
    */
   virtual std::int64_t next_long(std::int64_t i, std::int64_t j) = 0;

   /**
    * @brief Returns a random `b`-bit integer, i.e. uniform over
    *        @f$\{0, \dots, 2^b - 1\}@f$.
    *
    * Default implementation, identical to the Java default method:
    * `nextLong(0, (1L << b) - 1)`. Generators may override it with something
    * faster, such as a right shift of one output word.
    *
    * Java semantics are reproduced exactly: Java takes the shift count modulo
    * 64 and the subtraction wraps around. We compute in `std::uint64_t` to get
    * the same bits without undefined behaviour. (In Java, `b = 63` makes
    * `nextLong` see a range of size @f$2^{63}@f$, and `b = 64` gives the
    * range @f$\{0\}@f$.)
    *
    * @param b number of random bits
    */
   virtual std::int64_t next_bits_long(int b) {
      const std::uint64_t upper = (std::uint64_t{1} << (b & 63)) - 1u;
      return next_long(0, static_cast<std::int64_t>(upper));
   }
};

}  // namespace sscpp::rng
