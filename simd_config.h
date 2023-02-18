#pragma once

#include <cstdint>
#include <cstddef>

// define FORCE_INLINE
#if defined(__GNUC__)
#   define FORCE_INLINE __attribute__((always_inline)) inline
#elif defined(_MSC_VER) || defined(__INTEL_COMPILER)
#   define FORCE_INLINE __forceinline
#else
#   define FORCE_INLINE inline
#endif

#if defined(_MSC_VER) && (_M_IX86_FP==2 || defined(_M_X64))
#  define __SSE2__
#  define __SSE4_1__
#  define __SSE4_2__
#endif

#if defined(__AVX512F__)
#   define SIMD_N_BITS 512
#elif defined(__AVX2__)
#   define SIMD_N_BITS 256
#elif defined(__AVX__)
#   error AVX2 is needed
#elif defined(__SSE4_1__)
#   define SIMD_N_BITS 128
#endif

#ifdef SIMD_N_BITS
#   include <immintrin.h>
#   ifdef _MSC_VER
#       include <intrin.h>
#   endif
#else
#error "SIMD_N_BITS not defined"
#endif

inline int64_t popcnt(int64_t x)
{
#ifdef _MSC_VER
    return __popcnt64(x);
#else
    return __builtin_popcountll(x);
#endif
}