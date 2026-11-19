/**
 * @file mul64.hpp
 * @brief Full 64 x 64 -> 128-bit unsigned multiplication (internal header).
 *
 * Java has no 128-bit type and uses `Math.unsignedMultiplyHigh`; standard
 * C++20 has none either. `mul_64x64_128` uses the fastest instruction
 * available and falls back to standard C++ elsewhere:
 * - GCC and Clang (64-bit targets): `unsigned __int128`;
 * - MSVC on x64: `_umul128`; MSVC on ARM64: `__umulh`;
 * - any other compiler: `mul_64x64_128_portable` (32-bit halves).
 * All versions give exactly the same bits. The portable version is always
 * compiled, so the tests can compare it with the fast one.
 */
#pragma once

#include <cstdint>

#if defined(_MSC_VER) && !defined(__SIZEOF_INT128__) && (defined(_M_X64) || defined(_M_ARM64))
#include <intrin.h>
#endif

namespace sscpp::rng::detail {

/**
 * @brief Standard C++ version: returns the low 64 bits of @f$a b@f$ and
 *        stores the high 64 bits in `hi`.
 *
 * Splits each factor into 32-bit halves: four 32 x 32 -> 64-bit products,
 * none of which can overflow.
 */
constexpr std::uint64_t mul_64x64_128_portable(std::uint64_t a, std::uint64_t b, std::uint64_t& hi) noexcept {
   constexpr std::uint64_t kLow32 = 0xFFFFFFFFULL;
   const std::uint64_t a_lo = a & kLow32, a_hi = a >> 32;
   const std::uint64_t b_lo = b & kLow32, b_hi = b >> 32;
   const std::uint64_t p0 = a_lo * b_lo;
   const std::uint64_t p1 = a_lo * b_hi;
   const std::uint64_t p2 = a_hi * b_lo;
   const std::uint64_t p3 = a_hi * b_hi;
   // Bits 32..95 of the product, before carrying into the high word: < 3 * 2^32.
   const std::uint64_t mid = (p0 >> 32) + (p1 & kLow32) + (p2 & kLow32);
   hi = p3 + (p1 >> 32) + (p2 >> 32) + (mid >> 32);
   return (mid << 32) | (p0 & kLow32);
}

/**
 * @brief Returns the low 64 bits of @f$a b@f$ and stores the high 64 bits
 *        in `hi`, with the fastest method of the compiler.
 */
inline std::uint64_t mul_64x64_128(std::uint64_t a, std::uint64_t b, std::uint64_t& hi) noexcept {
#if defined(__SIZEOF_INT128__)
   // `__extension__` silences -Wpedantic: __int128 is a GCC/Clang extension.
   __extension__ typedef unsigned __int128 U128;
   const U128 p = static_cast<U128>(a) * b;
   hi = static_cast<std::uint64_t>(p >> 64);
   return static_cast<std::uint64_t>(p);
#elif defined(_MSC_VER) && defined(_M_X64)
   return _umul128(a, b, &hi);
#elif defined(_MSC_VER) && defined(_M_ARM64)
   hi = __umulh(a, b);
   return a * b;
#else
   return mul_64x64_128_portable(a, b, hi);
#endif
}

}  // namespace sscpp::rng::detail
