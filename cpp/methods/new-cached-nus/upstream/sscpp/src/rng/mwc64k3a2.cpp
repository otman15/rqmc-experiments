/**
 * @file mwc64k3a2.cpp
 * @brief Implementation of `sscpp::rng::MWC64k3a2`.
 *
 * Port of `umontreal.ssj.rng.MWC64k3a2`, with the same states and outputs.
 * Java `long` values are kept in `std::uint64_t`, which gives the same bits
 * and the unsigned reading Java gets with `Long.compareUnsigned` and
 * `Math.unsignedMultiplyHigh`. 128-bit products use `detail::mul_64x64_128`.
 *
 * Jumps. With b = 2^64, the MWC is equivalent to the multiplicative LCG
 *    y_{n+1} = b^-1 * y_n mod m,   m = A3*b^3 + A2*b^2 - 1,
 * through the maps `state_to_lcg` and `lcg_to_state` below (the same maps as
 * Java `stateToLCG` and `stateFromLCGParts4`). Jumping n steps ahead is
 * y -> (b^-1)^n * y mod m. The LCG states are numbers below m < 2^235, stored
 * as 4 little-endian 64-bit words.
 *
 * Multiplications modulo m use Montgomery multiplication with radix
 * R = b^4 = 2^256: mont_mul(x, y) = x * y * R^-1 mod m. Since m = -1 mod b,
 * the Montgomery quotient digit is simply the current low word, and adding a
 * multiple of m only touches two words (see `mont_reduce`).
 *
 * The LCG states are stored in normal form (Java stores them in Montgomery
 * form) and the jump multipliers J in Montgomery form J*R mod m, so a jump is
 * one call mont_mul(J*R, y) = J*y mod m. The value of y is the same as in
 * Java, so the states are identical.
 */
#include "sscpp/rng/mwc64k3a2.hpp"

#include <bit>
#include <stdexcept>
#include <utility>

#include "mul64.hpp"
#include "mwc64k3a2_constants.hpp"

namespace sscpp::rng {

namespace {

using Words = std::array<std::uint64_t, 4>;
using Words8 = std::array<std::uint64_t, 8>;
namespace constants = detail::mwc64k3a2;

/// Coefficients a2 and a3 of the recurrence (Java `A2`, `A3`).
constexpr std::uint64_t kA2 = 184698970548483715ULL;
constexpr std::uint64_t kA3 = 6028691832887ULL;

/// Largest allowed carry in a seed (Java `MAX_CARRY`).
constexpr std::uint64_t kMaxCarry = kA2 + kA3 - 1;

/// 2^-53 (Java `NORM53`).
constexpr double kNorm53 = 0x1.0p-53;

/// m = A3*b^3 + (A2 - 1)*b^2 + (b - 1)*b + (b - 1) (Java `MOD0..MOD3`).
constexpr Words kM = {~0ULL, ~0ULL, kA2 - 1, kA3};

/// b^-1 in Montgomery form: b^-1 * b^4 = b^3 mod m, and b^3 < m.
constexpr Words kBInvMont = {0, 0, 0, 1};

/// Returns t + a*b + carry (low word) and sets carry to the high word.
/// Cannot overflow: (2^64 - 1)^2 + 2 (2^64 - 1) = 2^128 - 1.
inline std::uint64_t mul_add(std::uint64_t a, std::uint64_t b, std::uint64_t t, std::uint64_t& carry) {
   std::uint64_t hi;
   std::uint64_t lo = detail::mul_64x64_128(a, b, hi);
   lo += t;
   hi += (lo < t);
   lo += carry;
   hi += (lo < carry);
   carry = hi;
   return lo;
}

/**
 * Montgomery reduction: returns T * R^-1 mod m, for T < m*R given as 8 words.
 *
 * Round i adds q*m*b^i with q = t[i], which makes word i zero. Since
 * q*m = q*A3*b^3 + q*A2*b^2 - q, this subtracts q from t[i] (giving 0, no
 * borrow) and adds q*A2 at word i+2 and q*A3 at word i+3. After 4 rounds the
 * low 4 words are 0, and the high 4 words are (T + Q*m) / R < 2m (since
 * T < m*R and Q < R); one conditional subtraction of m gives the result.
 * Everything fits in 8 words because T + Q*m < 2*m*R < R^2.
 */
Words mont_reduce(Words8 t) {
   for (int i = 0; i < 4; i++) {
      const std::uint64_t q = t[i];
      t[i] = 0;
      std::uint64_t carry = 0;
      t[i + 2] = mul_add(q, kA2, t[i + 2], carry);
      t[i + 3] = mul_add(q, kA3, t[i + 3], carry);
      for (int j = i + 4; j < 8 && carry != 0; j++) {
         t[j] += carry;
         carry = (t[j] < carry);
      }
   }
   Words r = {t[4], t[5], t[6], t[7]};
   // r >= m ?  Compare from the most significant word.
   bool ge = true;
   for (int i = 3; i >= 0; i--) {
      if (r[i] != kM[i]) {
         ge = r[i] > kM[i];
         break;
      }
   }
   if (ge) {
      std::uint64_t borrow = 0;
      for (int i = 0; i < 4; i++) {
         const std::uint64_t d = r[i] - kM[i] - borrow;
         borrow = (r[i] < kM[i]) || (r[i] == kM[i] && borrow != 0);
         r[i] = d;
      }
   }
   return r;
}

/// Montgomery multiplication: returns a * b * R^-1 mod m, for a, b < m.
Words mont_mul(const Words& a, const Words& b) {
   Words8 t{};
   for (int i = 0; i < 4; i++) {
      std::uint64_t carry = 0;
      for (int j = 0; j < 4; j++)
         t[i + j] = mul_add(a[i], b[j], t[i + j], carry);
      t[i + 4] = carry;
   }
   return mont_reduce(t);
}

/**
 * MWC state {x3, x2, x1, carry} -> LCG state y (normal form), as Java
 * `stateToLCG`: y = (A3*b^3*x3 + b*x2 + b^2*x1 + b^3*carry) mod m.
 *
 * The integer V = b*x2 + b^2*x1 + b^3*(carry + A3*x3) has 5 words and can be
 * larger than m (for example when the carry is large). It is reduced with
 * Montgomery arithmetic: mont_reduce(V) = V*R^-1 mod m, then
 * mont_mul(V*R^-1, R^2 mod m) = V mod m.
 */
Words state_to_lcg(const Words& s) {
   std::uint64_t hi;
   std::uint64_t lo = detail::mul_64x64_128(kA3, s[0], hi);
   lo += s[3];
   hi += (lo < s[3]);
   return mont_mul(mont_reduce(Words8{0, s[1], s[2], lo, hi, 0, 0, 0}), constants::kR2);
}

/**
 * LCG state y (normal form) -> MWC state {x3, x2, x1, carry}, as Java
 * `stateFromLCGParts4`: x3 = y0, x2 = y1, and x1 + b*carry = A2*y0 + y2 + b*y3.
 */
Words lcg_to_state(const Words& y) {
   std::uint64_t hi;
   const std::uint64_t lo = detail::mul_64x64_128(kA2, y[0], hi);
   const std::uint64_t x1 = lo + y[2];
   hi += (x1 < lo);
   return {y[0], y[1], x1, hi + y[3]};
}

/// Seed (Java long bits) -> words.
Words to_words(const MWC64k3a2::Seed& seed) {
   return {static_cast<std::uint64_t>(seed[0]), static_cast<std::uint64_t>(seed[1]),
           static_cast<std::uint64_t>(seed[2]), static_cast<std::uint64_t>(seed[3])};
}

/// "{ a, b, c, d }" with unsigned words (Java `Long.toUnsignedString`).
std::string words_to_string(std::uint64_t a, std::uint64_t b, std::uint64_t c, std::uint64_t d) {
   return "{ " + std::to_string(a) + ", " + std::to_string(b) + ", " + std::to_string(c) + ", "
         + std::to_string(d) + " }";
}

}  // namespace

// Default package seed {1, 3, 4, 5} (Java `nextSeed`) and its LCG state,
// computed by tools/mwc64k3a2_constants.py (constant initialization, no setup
// code runs before main).
MWC64k3a2::Words MWC64k3a2::next_seed_ = {1, 3, 4, 5};
MWC64k3a2::Words MWC64k3a2::next_seed_y_ = constants::kDefaultSeedY;

MWC64k3a2::MWC64k3a2() {
   // This stream starts at the current package seed...
   stream_ = next_seed_;
   stream_y_ = next_seed_y_;
   reset_start_stream();
   // ...and the next created stream will start 2^169 steps further.
   next_seed_y_ = mont_mul(constants::kStreamJump, next_seed_y_);
   next_seed_ = lcg_to_state(next_seed_y_);
}

MWC64k3a2::MWC64k3a2(std::string name) : MWC64k3a2() {
   name_ = std::move(name);
}

void MWC64k3a2::check_seed(const Seed& seed) {
   // Java also checks for null and seed.length != 4; here the array type
   // guarantees both.
   if (seed[3] < 0 || seed[3] > static_cast<std::int64_t>(kMaxCarry))
      throw std::invalid_argument("The carry must be in [0, " + std::to_string(kMaxCarry) + "].");
   if (seed[0] == 0 && seed[1] == 0 && seed[2] == 0 && seed[3] == 0)
      throw std::invalid_argument("The all-zero state is not allowed.");
   if (seed[0] == -1 && seed[1] == -1 && seed[2] == -1 && seed[3] == static_cast<std::int64_t>(kMaxCarry))
      throw std::invalid_argument("The all-ones/max-carry state is not allowed.");
}

void MWC64k3a2::set_package_seed(const Seed& seed) {
   check_seed(seed);
   next_seed_ = to_words(seed);
   next_seed_y_ = state_to_lcg(next_seed_);
}

void MWC64k3a2::set_seed(const Seed& seed) {
   check_seed(seed);
   stream_ = to_words(seed);
   stream_y_ = state_to_lcg(stream_);
   reset_start_stream();
}

MWC64k3a2::Seed MWC64k3a2::get_state() const {
   return {static_cast<std::int64_t>(x3_), static_cast<std::int64_t>(x2_), static_cast<std::int64_t>(x1_),
           static_cast<std::int64_t>(carry_)};
}

std::unique_ptr<MWC64k3a2> MWC64k3a2::clone() const {
   // The default copy constructor copies the state and the saved states,
   // which is the deep copy done by Java clone().
   return std::make_unique<MWC64k3a2>(*this);
}

void MWC64k3a2::reset_start_stream() {
   substream_ = stream_;
   substream_y_ = stream_y_;
   reset_start_substream();
}

void MWC64k3a2::reset_start_substream() {
   x3_ = substream_[0];
   x2_ = substream_[1];
   x1_ = substream_[2];
   carry_ = substream_[3];
}

void MWC64k3a2::reset_next_substream() {
   substream_y_ = mont_mul(constants::kSubstreamJump, substream_y_);
   substream_ = lcg_to_state(substream_y_);
   reset_start_substream();
}

void MWC64k3a2::advance_state_by_jump(std::uint64_t n) {
   if (n == 0)
      return;
   const Words y = state_to_lcg({x3_, x2_, x1_, carry_});
   // p = (b^-1)^n in Montgomery form, by square-and-multiply from the top bit of n.
   Words p = kBInvMont;
   for (int bit = static_cast<int>(std::bit_width(n)) - 2; bit >= 0; bit--) {
      p = mont_mul(p, p);
      if ((n >> bit) & 1u)
         p = mont_mul(p, kBInvMont);
   }
   const Words s = lcg_to_state(mont_mul(p, y));  // (b^-1)^n * y mod m
   x3_ = s[0];
   x2_ = s[1];
   x1_ = s[2];
   carry_ = s[3];
}

std::string MWC64k3a2::to_string() const {
   std::string s = "The current state of MWC64k3a2";
   if (!name_.empty())
      s += " " + name_;
   return s + " is: " + words_to_string(x3_, x2_, x1_, carry_);
}

std::string MWC64k3a2::to_string_full() const {
   std::string s = "MWC64k3a2 stream";
   if (!name_.empty())
      s += " " + name_;
   s += ":\n";
   s += " Ig = " + words_to_string(stream_[0], stream_[1], stream_[2], stream_[3]) + "\n";
   s += " Bg = " + words_to_string(substream_[0], substream_[1], substream_[2], substream_[3]) + "\n";
   s += " Cg = " + words_to_string(x3_, x2_, x1_, carry_) + "\n";
   return s;
}

std::uint64_t MWC64k3a2::next_number() {
   // t = A2*x2 + A3*x3 + carry < 2^122: the new x1 is its low word and the
   // new carry its high word.
   const std::uint64_t out = x1_;
   std::uint64_t hi2, hi3;
   const std::uint64_t lo2 = detail::mul_64x64_128(kA2, x2_, hi2);
   const std::uint64_t lo3 = detail::mul_64x64_128(kA3, x3_, hi3);
   std::uint64_t lo = lo2 + lo3;
   std::uint64_t hi = hi2 + hi3 + (lo < lo2);
   lo += carry_;
   hi += (lo < carry_);
   x3_ = x2_;
   x2_ = x1_;
   x1_ = lo;
   carry_ = hi;
   return out;
}

double MWC64k3a2::next_double() {
   // Java RandomStreamBase.nextDouble() with its defaults simply returns
   // MWC64k3a2.nextValue(), which is the code below. The value is below 2^53,
   // so its conversion to double is exact.
   return static_cast<double>(next_number() >> 11) * kNorm53;
}

int MWC64k3a2::next_int(int i, int j) {
   return static_cast<int>(next_long(static_cast<std::int64_t>(i), static_cast<std::int64_t>(j)));
}

std::int64_t MWC64k3a2::next_long(std::int64_t i, std::int64_t j) {
   if (i > j)
      throw std::invalid_argument(std::to_string(i) + " is larger than " + std::to_string(j) + ".");
   // n = j - i + 1, wrapping around modulo 2^64 like Java long arithmetic.
   // The sums below are also done in unsigned arithmetic (same bits as Java,
   // no C++ signed-overflow undefined behaviour).
   const std::uint64_t n = static_cast<std::uint64_t>(j) - static_cast<std::uint64_t>(i) + 1u;
   if (static_cast<std::int64_t>(n) > 0) {
      std::uint64_t r = next_number() >> 1;  // 63 random bits.
      const std::uint64_t m = n - 1u;
      if ((n & m) == 0u)  // n is a power of 2: no rejection.
         return static_cast<std::int64_t>(static_cast<std::uint64_t>(i) + (r & m));
      // Java: while (u + m - (r = u % n) < 0L). The test is negative exactly
      // when u lies in the last, incomplete block of n values below 2^63
      // (the signed sum wraps around).
      std::uint64_t u = r;
      while (static_cast<std::int64_t>(u + m - (r = u % n)) < 0)
         u = next_number() >> 1;
      return static_cast<std::int64_t>(static_cast<std::uint64_t>(i) + r);
   }
   // Range size larger than 2^63 - 1: draw full words until one is in [i, j].
   std::int64_t r;
   do {
      r = static_cast<std::int64_t>(next_number());
   } while (r < i || r > j);
   return r;
}

std::int64_t MWC64k3a2::next_bits_long(int b) {
   // Java: nextNumber() >>> (64 - b), where Java takes the shift count mod 64.
   // Unsigned subtraction gives the same count without signed overflow (b = INT_MIN).
   return static_cast<std::int64_t>(next_number() >> ((64u - static_cast<unsigned>(b)) & 63u));
}

std::int64_t MWC64k3a2::next_raw() {
   return static_cast<std::int64_t>(next_number());
}

}  // namespace sscpp::rng
