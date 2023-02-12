#pragma once

#include "simd_config.h"

template <size_t NumBits>
struct SimdRegister;

template <>
struct SimdRegister<128>
{
    __m128i m_v;

    typedef SimdRegister<128> XV;

    SimdRegister() {}
    SimdRegister(const void* p) : m_v(_mm_loadu_si128((const __m128i*)p)) {}
    SimdRegister(__m128i v) : m_v(v) {}

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return _mm_and_si128(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return _mm_xor_si128(a.m_v, b.m_v); }

    static FORCE_INLINE XV zero() { return _mm_setzero_si128(); }

    union Konst
    {
        Konst() : m_simd(zero().m_v) {}
        Konst(__m128i v) : m_simd(v) {}

        void operator^=(XV v) { m_simd = (XV(m_simd) ^ v).m_v; }
        uint8_t parity() const
        {
            return popcnt(m_u64[0] ^ m_u64[1]) & 0x1;
        }
        __m128i m_simd;
        uint64_t m_u64[2];
    };
};

#if SIMD_N_BITS>=256
template <>
struct SimdRegister<256>
{
    __m256i m_v;

    typedef SimdRegister<256> XV;

    SimdRegister() {}
    SimdRegister(const void* p) : m_v(_mm256_load_si256((const __m256i*)p)) {}
    SimdRegister(__m256i v) : m_v(v) {}

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return _mm256_and_si256(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return _mm256_xor_si256(a.m_v, b.m_v); }

    static FORCE_INLINE XV zero() { return _mm256_setzero_si256(); }

    union Konst
    {
        Konst() : m_simd(zero().m_v) {}
        Konst(__m256i v) : m_simd(v) {}

        void operator^=(XV v) { m_simd = (XV(m_simd) ^ v).m_v; }

        uint8_t parity() const
        {
            typedef SimdRegister<128> v128_t;
            return v128_t::Konst(v128_t(_mm_xor_si128(m_u128[0], m_u128[1])).m_v).parity();
        }
        __m256i m_simd;
        __m128i m_u128[2];
    };
};
#endif

#if SIMD_N_BITS>=512
template <>
struct SimdRegister<512>
{
    __m512i m_v;

    typedef SimdRegister<512> XV;

    SimdRegister() {}
    SimdRegister(const void* p) : m_v(_mm512_load_si512((const __m512i*)p)) {}
    SimdRegister(__m512i v) : m_v(v) {}

    //    static FORCE_INLINE XV load(const void* p) { return _mm512_loadu_si512((const __m512i*) p); }

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return _mm512_and_si512(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return _mm512_xor_si512(a.m_v, b.m_v); }

    static FORCE_INLINE XV zero() { return _mm512_setzero_si512(); }

    union Konst
    {
        Konst() : m_simd(zero().m_v) {}

        void operator^=(XV v) { m_simd = (XV(m_simd) ^ v).m_v; }

        uint8_t parity() const
        {
            typedef SimdRegister<256> v256_t;
            return v256_t::Konst(v256_t(_mm256_xor_si256(m_u256[0], m_u256[1])).m_v).parity();
        }
        __m512i m_simd;
        __m256i m_u256[2];
    };
};
#endif




