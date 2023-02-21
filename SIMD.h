#pragma once

#include "simd_config.h"

template <size_t NumBits>
struct SimdRegister;

template <>
struct SimdRegister<32>
{
    uint32_t m_v;

    typedef SimdRegister<32> XV;

    SimdRegister() {}
    SimdRegister(const void* p) : m_v(*(const uint32_t*)p) {}
    SimdRegister(uint32_t v) : m_v(v) {}

    friend FORCE_INLINE XV operator&(const XV a, const XV b) { return a.m_v & b.m_v; }
    friend FORCE_INLINE XV operator^(const XV a, const XV b) { return a.m_v ^ b.m_v; }
    friend FORCE_INLINE XV operator|(const XV a, const XV b) { return a.m_v | b.m_v; }
    //friend FORCE_INLINE XV operator>(const XV& a, const XV& b) { return a.m_v > b.m_v ? uint32_t(-1) : uint32_t(0); }
    friend FORCE_INLINE XV operator>>(const XV a, int n) { return uint32_t(a.m_v >> n); }
    friend FORCE_INLINE XV operator<<(const XV a, int n) { return uint32_t(a.m_v << n); }

    FORCE_INLINE XV ifOddValueElseZero(const XV value) const
    {
#if defined(__GNUC__) && defined(__x86_64__)
        // force the use of cmov with gcc
        uint32_t z;
        __asm__(
            "mov %[a], %[z]\n"
            "and $0x1, %[z]\n"
            "cmovne %[b], %[z]\n"
            : [z] "=r"(z)
            : [a] "r"(m_v), [b] "r"(value.m_v)
            : "cc"
        );
        return z;
#else
        const uint32_t lowestBit = m_v & 0x1;
        return lowestBit ? value.m_v : 0;
#endif
    }

    static FORCE_INLINE XV zero() { return uint32_t(0); }

    struct Konst
    {
        Konst() : m_simd(0) {}
        Konst(uint32_t v) : m_simd(v) {}

        void operator^=(XV v) { m_simd ^= v.m_v; }

        uint8_t parity() const
        {
            return popcnt(m_simd) % 2;
        }
        uint32_t m_simd;
    };
};

template <>
struct SimdRegister<64>
{
    uint64_t m_v;

    typedef SimdRegister<64> XV;

    SimdRegister() {}
    SimdRegister(const void* p) : m_v(*(const uint64_t*)p) {}
    SimdRegister(uint32_t v) : m_v(v | (uint64_t(v) << 32)) {}
    SimdRegister(uint64_t v) : m_v(v) {}

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return a.m_v & b.m_v; }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return a.m_v ^ b.m_v; }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { return a.m_v | b.m_v; }
    friend FORCE_INLINE XV operator>>(const XV& a, int n)
    {
        uint64_t mask = (uint64_t(0xFFFFFFFF) << 32) | ((uint32_t(1) << (32 - n)) - 1);
        return uint64_t(a.m_v >> n) & mask;
    }
    friend FORCE_INLINE XV operator<<(const XV& a, int n)
    {
        uint64_t mask = uint64_t(-1) & ~(((uint64_t(1) << n) - 1) << 32);
        return uint64_t(a.m_v << n);
    }

    FORCE_INLINE XV ifOddValueElseZero(const XV& value) const
    {
        const uint64_t maskHi = m_v & (uint64_t(1) << 32) ? (uint64_t(0xFFFFFFFF) << 32) : 0;
        const uint64_t maskLo = m_v & 0x1 ? uint64_t(0xFFFFFFFF) : 0;
        return value.m_v & (maskHi | maskLo);
    }

    static FORCE_INLINE XV zero() { return uint64_t(0); }

    struct Konst
    {
        Konst() : m_simd(0) {}
        Konst(uint64_t v) : m_simd(v) {}

        void operator^=(XV v) { m_simd ^= v.m_v; }

        uint8_t parity() const
        {
            return popcnt(m_simd) % 2;
        }
        uint64_t m_simd;
    };
};

#if SIMD_N_BITS>=128
template <>
struct SimdRegister<128>
{
    __m128i m_v;

    typedef SimdRegister<128> XV;

    SimdRegister() {}
    SimdRegister(uint32_t v) : m_v(_mm_set1_epi32(v)){}
    SimdRegister(const void* p) : m_v(_mm_loadu_si128((const __m128i*)p)) {}
    SimdRegister(__m128i v) : m_v(v) {}

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return _mm_and_si128(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return _mm_xor_si128(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { return _mm_or_si128(a.m_v, b.m_v); }
    //friend FORCE_INLINE XV operator>(const XV& a, const XV& b) { return _mm_cmpgt_epi32(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n) { return _mm_slli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n) { return _mm_srli_epi32(a.m_v, n); }

    static FORCE_INLINE XV zero() { return _mm_setzero_si128(); }

    FORCE_INLINE XV ifOddValueElseZero(const XV& value) const
    {
        const __m128i z = zero().m_v;
        const __m128i lowestBit = _mm_slli_epi32(m_v, 31);
        const __m128i isOdd = _mm_cmpgt_epi32(z, lowestBit);
        return _mm_and_si128(isOdd, value.m_v);
    }

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
#endif

#if SIMD_N_BITS>=256
template <>
struct SimdRegister<256>
{
    __m256i m_v;

    typedef SimdRegister<256> XV;

    SimdRegister() {}
    SimdRegister(const void* p) : m_v(_mm256_load_si256((const __m256i*)p)) {}
    SimdRegister(uint32_t v) : m_v(_mm256_set1_epi32(v)) {}
    SimdRegister(__m256i v) : m_v(v) {}

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return _mm256_and_si256(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return _mm256_xor_si256(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { return _mm256_or_si256(a.m_v, b.m_v); }
    //friend FORCE_INLINE XV operator>(const XV& a, const XV& b) { return _mm256_cmpgt_epi32(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n) { return _mm256_slli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n) { return _mm256_srli_epi32(a.m_v, n); }

    static FORCE_INLINE XV zero() { return _mm256_setzero_si256(); }

    FORCE_INLINE XV ifOddValueElseZero(const XV& value) const
    {
        const __m256i z = zero().m_v;
        const __m256i lowestBit = _mm256_slli_epi32(m_v, 31);
        const __m256i isOdd = _mm256_cmpgt_epi32(z, lowestBit);
        return _mm256_and_si256(isOdd, value.m_v);
    }

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
    SimdRegister(uint32_t v) : m_v(_mm512_set1_epi32(v)) {}
    SimdRegister(__m512i v) : m_v(v) {}

    //    static FORCE_INLINE XV load(const void* p) { return _mm512_loadu_si512((const __m512i*) p); }

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return _mm512_and_si512(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return _mm512_xor_si512(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { return _mm512_or_si512(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n) { return _mm512_slli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n) { return _mm512_srli_epi32(a.m_v, n); }

    FORCE_INLINE XV ifOddValueElseZero(const XV& value) const
    {
        const __m512i z = zero().m_v;
        const __m512i lowestBit = _mm512_slli_epi32(m_v, 31);
        const __mmask16 isOdd = _mm512_cmpneq_epi32_mask(lowestBit, z);
        return _mm512_mask_blend_epi32(isOdd, z, value.m_v);
    }

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




