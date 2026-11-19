/**
 * @file mwc64k3a2.hpp
 * @brief The `MWC64k3a2` generator: multiply-with-carry, base 2^64, order 3.
 *
 * C++ counterpart of `umontreal.ssj.rng.MWC64k3a2` (SSJ).
 */
#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <string>

#include "sscpp/rng/random_stream.hpp"

namespace sscpp::rng {

/**
 * @brief Multiply-with-carry (MWC) generator with base @f$b = 2^{64}@f$ and
 *        order @f$k = 3@f$, with streams and substreams.
 *
 * Bit-exact port of the Java class `umontreal.ssj.rng.MWC64k3a2`: same public
 * methods, same states and same outputs as Java for the same seeds. The code
 * is not a line-by-line copy: Java's emulation of 128-bit products, its
 * BigInteger setup and its hand-unrolled multi-word arithmetic are replaced
 * by simpler C++ that gives the same bits.
 *
 * The recurrence, with @f$a_0 = -1@f$, @f$a_1 = 0@f$,
 * @f$a_2 = 184698970548483715@f$ and @f$a_3 = 6028691832887@f$, is
 * @f[
 *    t = a_2 x_{n-2} + a_3 x_{n-3} + c_{n-1}, \qquad
 *    x_n = t \bmod 2^{64}, \qquad c_n = \lfloor t / 2^{64} \rfloor .
 * @f]
 * Each step returns the previous @f$x_{n-1}@f$ as a 64-bit output word. The
 * state is @f$\{x_{n-3}, x_{n-2}, x_{n-1}, c_{n-1}\}@f$ (this order is used
 * for seeds and for `get_state()`).
 *
 * **Streams and substreams.** The MWC is equivalent to a multiplicative LCG
 * with modulus @f$m = a_3 b^3 + a_2 b^2 - 1@f$, in which one MWC step
 * multiplies the LCG state @f$y@f$ by @f$b^{-1} \bmod m@f$. Streams start
 * @f$2^{169}@f$ steps apart and substreams @f$2^{118}@f$ steps apart; a jump
 * is one multiplication modulo @f$m@f$ by a fixed multiplier (Montgomery
 * multiplication, see `mwc64k3a2.cpp`). The fixed constants are computed once
 * by `tools/mwc64k3a2_constants.py`. The default seed of the first stream is
 * @f$\{1, 3, 4, 5\}@f$.
 *
 * In Java the class extends `RandomStreamBase`; here it implements
 * `RandomStream` directly, and `next_double()` contains what Java splits
 * between `RandomStreamBase.nextDouble()` and `MWC64k3a2.nextValue()` (with
 * the Java defaults: no increased precision, no antithetic output). Not
 * ported: `increasedPrecision`, `setAntithetic` and the `nextArrayOf...`
 * methods of `RandomStreamBase`.
 *
 * Seeds are given as `std::int64_t` (Java `long`); the words are read as
 * unsigned 64-bit integers, as in Java. Copying an `MWC64k3a2` object gives an
 * independent generator with the same state, which is what Java `clone()`
 * does.
 */
class MWC64k3a2 : public RandomStream {
public:
   /// A seed or state @f$\{x_{n-3}, x_{n-2}, x_{n-1}, c_{n-1}\}@f$, as Java `long[4]`.
   using Seed = std::array<std::int64_t, 4>;

   /**
    * @brief Constructs a new stream.
    *
    * Its initial state is the current package seed. The package seed is then
    * advanced by @f$2^{169}@f$ steps, so the next constructed stream starts
    * @f$2^{169}@f$ steps further.
    */
   MWC64k3a2();

   /**
    * @brief Constructs a new stream with the name `name`, printed by
    *        `to_string()` and `to_string_full()`.
    *
    * As in Java, an empty name is treated as no name.
    *
    * @param name name of the stream
    */
   explicit MWC64k3a2(std::string name);

   /**
    * @brief Sets the initial seed of the next created stream (Java
    *        `setPackageSeed`).
    *
    * The carry `seed[3]` must be in @f$[0, a_2 + a_3 - 1]@f$; the all-zero
    * state and the state with all words equal to @f$2^{64}-1@f$ and the
    * maximal carry are not allowed.
    *
    * The package seed is shared by all `MWC64k3a2` objects (static, as in
    * Java). It is not thread-safe.
    *
    * @param seed @f$\{x_{n-3}, x_{n-2}, x_{n-1}, c_{n-1}\}@f$
    * @throws std::invalid_argument if the seed does not satisfy the conditions
    */
   static void set_package_seed(const Seed& seed);

   /**
    * @brief Restarts this stream from the given seed (Java `setSeed`).
    *
    * Only this stream is changed, so streams are no longer spaced
    * @f$2^{169}@f$ steps apart. The seed must satisfy the same conditions as
    * in `set_package_seed()`.
    *
    * @param seed @f$\{x_{n-3}, x_{n-2}, x_{n-1}, c_{n-1}\}@f$
    * @throws std::invalid_argument if the seed does not satisfy the conditions
    */
   void set_seed(const Seed& seed);

   /**
    * @brief Returns the current state
    *        @f$\{x_{n-3}, x_{n-2}, x_{n-1}, c_{n-1}\}@f$ (Java `getState`).
    */
   Seed get_state() const;

   /**
    * @brief Returns an independent copy of this generator (Java `clone`).
    *
    * Same as the copy constructor, but returns a heap object.
    */
   std::unique_ptr<MWC64k3a2> clone() const;

   void reset_start_stream() override;
   void reset_start_substream() override;

   /**
    * @copydoc RandomStream::reset_next_substream
    *
    * The jump ahead of @f$2^{118}@f$ steps is one multiplication modulo
    * @f$m@f$ of the LCG state of @f$B_g@f$.
    */
   void reset_next_substream() override;

   /**
    * @brief Returns "The current state of MWC64k3a2 [name] is: { x3, x2, x1, carry }",
    *        with the words printed as unsigned integers, as in Java.
    */
   std::string to_string() const override;

   /**
    * @brief Returns the stream name and the states @f$I_g@f$, @f$B_g@f$ and
    *        @f$C_g@f$, one per line (Java `toStringFull`).
    *
    * Lines end with `"\n"` (Java uses `System.lineSeparator()`, which is
    * `"\r\n"` on Windows).
    */
   std::string to_string_full() const;

   /**
    * @brief Returns the top 53 bits of the next output word times
    *        @f$2^{-53}@f$, a number in @f$[0, 1)@f$.
    *
    * As in Java, the result can be 0.0 (which the `RandomStream` contract
    * excludes) but never 1.0.
    */
   double next_double() override;

   /**
    * @brief Returns `(int) next_long(i, j)`, as in Java.
    */
   int next_int(int i, int j) override;

   /**
    * @brief Returns a uniform integer in @f$\{i, \dots, j\}@f$.
    *
    * As in Java: if the range size @f$n = j - i + 1@f$ is at most
    * @f$2^{63} - 1@f$, uses the top 63 bits of output words with rejection
    * (no rejection when @f$n@f$ is a power of 2); otherwise draws full output
    * words until one falls in @f$[i, j]@f$.
    *
    * @throws std::invalid_argument if `i > j` (as Java)
    */
   std::int64_t next_long(std::int64_t i, std::int64_t j) override;

   /**
    * @brief Returns the top `b` bits of the next output word, i.e.
    *        `nextNumber() >>> (64 - b)` as in Java (shift count mod 64).
    */
   std::int64_t next_bits_long(int b) override;

   /**
    * @brief Advances the current state by @f$n@f$ steps (Java
    *        `advanceStateByJump`).
    *
    * Only the current state @f$C_g@f$ changes; a later reset restores the
    * saved stream or substream start. The result is the same as @f$n@f$ calls
    * of `next_raw()`, computed with at most about 128 multiplications modulo
    * @f$m@f$.
    *
    * Difference from Java: Java takes a `BigInteger` and throws
    * `IllegalArgumentException` when @f$n < 0@f$. Here @f$n@f$ is an unsigned
    * 64-bit integer, so jumps are limited to @f$n < 2^{64}@f$ and there is no
    * negative case. For these @f$n@f$ the result is bit-identical to Java.
    *
    * @param n number of steps; 0 does nothing
    */
   void advance_state_by_jump(std::uint64_t n);

   /**
    * @brief Returns the next raw 64-bit output word (Java `nextRaw`, for
    *        testing).
    */
   std::int64_t next_raw();

private:
   /// Four 64-bit words: an MWC state {x3, x2, x1, carry} or an LCG state y.
   using Words = std::array<std::uint64_t, 4>;

   /// One step of the recurrence; returns the old x_{n-1}.
   std::uint64_t next_number();

   /// Java `checkSeed`: throws `std::invalid_argument` on a forbidden seed.
   static void check_seed(const Seed& seed);

   std::string name_;  ///< Stream name; empty means no name.

   // Current state C_g: x_{n-1}, x_{n-2}, x_{n-3}, c_{n-1}.
   std::uint64_t x1_ = 0, x2_ = 0, x3_ = 0, carry_ = 0;

   Words stream_{};       ///< Initial state I_g of this stream, {x3, x2, x1, carry}.
   Words stream_y_{};     ///< LCG state of I_g (normal form, not Montgomery form).
   Words substream_{};    ///< Start B_g of the current substream.
   Words substream_y_{};  ///< LCG state of B_g.

   /// Package seed: initial state of the next created stream (Java `nextSeed`).
   static Words next_seed_;
   /// LCG state of the package seed.
   static Words next_seed_y_;
};

}  // namespace sscpp::rng
