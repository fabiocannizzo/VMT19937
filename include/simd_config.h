#pragma once

#include <cstdint>
#include <cstddef>

#include "macros.h"

namespace xvmt {

enum class ISA {
    Scalar,
    SSE2,
    SSE42,
    AVX2,
    AVX512,
    NEON,
    SVE256
};

template <ISA isa> struct IsaTraits;
template <> struct IsaTraits<ISA::Scalar> { static constexpr size_t s_hwBitLen = 32; };
template <> struct IsaTraits<ISA::SSE2>   { static constexpr size_t s_hwBitLen = 128; };
template <> struct IsaTraits<ISA::SSE42>  { static constexpr size_t s_hwBitLen = 128; };
template <> struct IsaTraits<ISA::AVX2>   { static constexpr size_t s_hwBitLen = 256; };
template <> struct IsaTraits<ISA::AVX512> { static constexpr size_t s_hwBitLen = 512; };
template <> struct IsaTraits<ISA::NEON>   { static constexpr size_t s_hwBitLen = 128; };
template <> struct IsaTraits<ISA::SVE256> { static constexpr size_t s_hwBitLen = 256; };

namespace details {

template <size_t Bits> struct BitLenToIsa;
template <> struct BitLenToIsa<32>  { static constexpr ISA s_isa = ISA::Scalar; };
template <> struct BitLenToIsa<128> {
#if defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__aarch64__) || defined(_M_ARM64) || defined(__arm__)
    static constexpr ISA s_isa = ISA::NEON;
#else
    static constexpr ISA s_isa = ISA::SSE42;
#endif
};
template <> struct BitLenToIsa<64>  { static constexpr ISA s_isa = ISA::Scalar; };
template <> struct BitLenToIsa<256> {
#if defined(__ARM_FEATURE_SVE) && (defined(__aarch64__) || defined(_M_ARM64))
    static constexpr ISA s_isa = ISA::SVE256;
#else
    static constexpr ISA s_isa = ISA::AVX2;
#endif
};
template <> struct BitLenToIsa<512> { static constexpr ISA s_isa = ISA::AVX512; };

} // namespace details

} // namespace xvmt

#if defined(_MSC_VER) && (_M_IX86_FP==2 || defined(_M_X64))
#  define __SSE2__
#  define __SSE4_1__
#  define __SSE4_2__
#endif

namespace xvmt {
namespace details {
    template <size_t Bits>
    struct BestIsa {
#if defined(__AVX512F__)
        static constexpr ISA s_isa = (Bits >= 512) ? ISA::AVX512 : (Bits >= 256 ? ISA::AVX2 : ISA::SSE42);
#elif defined(__AVX2__)
        static constexpr ISA s_isa = (Bits >= 256) ? ISA::AVX2 : ISA::SSE42;
#elif defined(__SSE4_2__)
        static constexpr ISA s_isa = ISA::SSE42;
#elif defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__aarch64__) || defined(_M_ARM64) || defined(__arm__)
#  if defined(__ARM_FEATURE_SVE) && (defined(__aarch64__) || defined(_M_ARM64))
        static constexpr ISA s_isa = (Bits >= 256) ? ISA::SVE256 : ISA::NEON;
#  else
        static constexpr ISA s_isa = ISA::NEON;
#  endif
#else
        static constexpr ISA s_isa = ISA::Scalar;
#endif
    };
} // namespace details
} // namespace xvmt

#ifndef SIMD_N_BITS
#if defined(__AVX512F__)
#   define SIMD_N_BITS 512
#   define SIMD_ISA ISA::AVX512
#elif defined(__AVX2__)
#   define SIMD_N_BITS 256
#   define SIMD_ISA ISA::AVX2
#elif defined(__AVX__)
#   error AVX2 is needed
#elif defined(__SSE4_2__)
#   define SIMD_N_BITS 128
#   define SIMD_ISA ISA::SSE42
#elif defined(__SSE2__)
#   error SSE4.2 is needed
#elif defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__aarch64__) || defined(_M_ARM64) || defined(__arm__)
#   if defined(__ARM_FEATURE_SVE) && (defined(__aarch64__) || defined(_M_ARM64))
#       define SIMD_N_BITS 256
#       define SIMD_ISA ISA::SVE256
#   else
#       define SIMD_N_BITS 128
#       define SIMD_ISA ISA::NEON
#   endif
#endif
#endif /* SIMD_N_BITS */

#ifndef SIMD_ISA
#  if SIMD_N_BITS == 512
#    define SIMD_ISA ISA::AVX512
#  elif SIMD_N_BITS == 256
#    if defined(__ARM_FEATURE_SVE) && (defined(__aarch64__) || defined(_M_ARM64))
#      define SIMD_ISA ISA::SVE256
#    else
#      define SIMD_ISA ISA::AVX2
#    endif
#  elif SIMD_N_BITS == 128
#    if defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__aarch64__) || defined(_M_ARM64) || defined(__arm__)
#      define SIMD_ISA ISA::NEON
#    else
#      define SIMD_ISA ISA::SSE42
#    endif
#  else
#    define SIMD_ISA ISA::Scalar
#  endif
#endif /* SIMD_ISA */

#ifdef SIMD_N_BITS
#   if defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__aarch64__) || defined(_M_ARM64) || defined(__arm__)
#       include <arm_neon.h>
#       if defined(__ARM_FEATURE_SVE)
#           include <arm_sve.h>
#       endif
#   else
#       include <immintrin.h>
#       ifdef _MSC_VER
#           include <intrin.h>
#       endif
#   endif
#else
#error "SIMD_N_BITS not defined"
#endif

