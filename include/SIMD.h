#pragma once

#include "simd_config.h"
#include "portable.h"

#include <algorithm>

template <size_t NumBits>
struct SimdRegister;

template <size_t NumBits>
struct SimdRegisterEmulator;

template <>
struct SimdRegister<32>
{
    uint32_t m_v;

    typedef SimdRegister<32> XV;

    SimdRegister() : m_v(0) {}
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
#if 0 && defined(__GNUC__) && defined(__x86_64__)
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
#elif defined(_MSC_VER)
        const uint32_t lowestBit = m_v & 0x1;
        return lowestBit ? value.m_v : 0;
#else
        const uint32_t x[2] = {0, value.m_v};
        const uint32_t lowestBit = m_v & 0x1;
        return x[lowestBit];
#endif
    }

    static FORCE_INLINE XV zero() { return uint32_t(0); }

    uint8_t parity() const { return popcnt(m_v) % 2; }
};

// This class is for debugging and testing only
// It emulates a regsiter with size Mx32 nbits, in case that is not natively available on the hardware
template <size_t M>
struct SimdRegisterEmulator
{
    struct A
    {
        A() : a{ {} } {}
        A(uint32_t v) { std::fill_n(a, M, v); }
        uint32_t a[M];
    };
    A m_v;

    typedef SimdRegisterEmulator<M> XV;

    SimdRegisterEmulator() {}
    SimdRegisterEmulator(uint32_t v) : m_v(v) {}
    SimdRegisterEmulator(const uint32_t* p) : m_v(*p) {}
    SimdRegisterEmulator(const A& v) : m_v(v) {}

    template <bool A>
    void store(uint32_t* dst) { std::copy_n(m_v.a, M, dst); }

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { XV r; for (size_t i = 0; i < M; ++i) r.m_v.a[i] = a.m_v.a[i] & b.m_v.a[i]; return r; }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { XV r; for (size_t i = 0; i < M; ++i) r.m_v.a[i] = a.m_v.a[i] ^ b.m_v.a[i]; return r; }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { XV r; for (size_t i = 0; i < M; ++i) r.m_v.a[i] = a.m_v.a[i] | b.m_v.a[i]; return r; }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n) { XV r; for (size_t i = 0; i < M; ++i) r.m_v.a[i] = a.m_v.a[i] << n; return r; }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n) { XV r; for (size_t i = 0; i < M; ++i) r.m_v.a[i] = a.m_v.a[i] >> n; return r; }

    static FORCE_INLINE XV zero() { return XV(uint32_t(0)); }

    FORCE_INLINE XV ifOddValueElseZero(const XV& value) const { XV r; for (size_t i = 0; i < M; ++i) r.m_v.a[i] = (m_v.a[i] % 2) ? value.m_v.a[i] : 0; return r; }

    uint8_t parity() const
    {
        size_t n;
        for (auto v : m_v.a)
            n += v;
        return popcnt(n) % 2;
    }
};

template <>
struct SimdRegister<64>
{
    uint64_t m_v;

    typedef SimdRegister<64> XV;

    SimdRegister() : m_v(0) {}
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

    uint8_t parity() const { return popcnt(m_v) % 2; }
};

#if SIMD_N_BITS>=128
template <>
struct SimdRegister<128>
{
    __m128i m_v;

    typedef SimdRegister<128> XV;

    SimdRegister() : m_v(_mm_undefined_si128()) {}
    SimdRegister(uint32_t v) : m_v(_mm_set1_epi32(v)){}
    SimdRegister(const void* p) : m_v(_mm_loadu_si128((const __m128i*)p)) {}
    SimdRegister(__m128i v) : m_v(v) {}

    template <bool A>
    void store(uint32_t* dst) { if (A) _mm_store_si128((__m128i*)dst, m_v); else _mm_storeu_si128((__m128i*)dst, m_v); }

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return _mm_and_si128(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return _mm_xor_si128(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { return _mm_or_si128(a.m_v, b.m_v); }
    //friend FORCE_INLINE XV operator>(const XV& a, const XV& b) { return _mm_cmpgt_epi32(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n) { return _mm_slli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n) { return _mm_srli_epi32(a.m_v, n); }

    static FORCE_INLINE XV zero() { return _mm_setzero_si128(); }

    FORCE_INLINE XV ifOddValueElseZero(const XV& value) const
    {
#if 0
        const __m128 z = _mm_setzero_ps();
        const __m128 lowestBit = _mm_castsi128_ps(_mm_slli_epi32(m_v, 31));
        return _mm_castps_si128(_mm_blendv_ps(z, _mm_castsi128_ps(value.m_v), lowestBit));
#else
        const __m128i z = zero().m_v;
        const __m128i lowestBit = _mm_slli_epi32(m_v, 31);
        const __m128i isOdd = _mm_cmpgt_epi32(z, lowestBit);
        return _mm_and_si128(isOdd, value.m_v);
#endif
    }

    uint8_t parity() const
    {
        __m128i hi(_mm_shuffle_epi32(m_v, 2 | (3 << 2)));
        __m128i mix = _mm_xor_si128(hi, m_v);
        double d = _mm_cvtsd_f64(_mm_castsi128_pd(mix));
        return popcnt(*(uint64_t*)(&d)) & 1;
    }
};
#endif

#if SIMD_N_BITS>=256
template <>
struct SimdRegister<256>
{
    __m256i m_v;

    typedef SimdRegister<256> XV;

    SimdRegister() : m_v(_mm256_undefined_si256()) {}
    SimdRegister(const void* p) : m_v(_mm256_load_si256((const __m256i*)p)) {}
    SimdRegister(uint32_t v) : m_v(_mm256_set1_epi32(v)) {}
    SimdRegister(__m256i v) : m_v(v) {}

    template <bool A>
    void store(uint32_t* dst) { if (A) _mm256_store_si256((__m256i*)dst, m_v); else _mm256_storeu_si256((__m256i*)dst, m_v); }

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
        const __m256i lowestBit = _mm256_slli_epi32(m_v, 31); // move least significant bit to most significant bit
#if 0
        const __m256i isOdd = _mm256_cmpgt_epi32(z, lowestBit);
        return _mm256_and_si256(value.m_v, isOdd);
#else
        const __m256 mask = _mm256_castsi256_ps(lowestBit);
        return _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(z), _mm256_castsi256_ps(value.m_v), mask));
#endif
    }

    uint8_t parity() const
    {
        __m128i hi(_mm256_extracti128_si256(m_v, 1));
        __m128i lo(_mm256_castsi256_si128(m_v));
        return SimdRegister<128>(_mm_xor_si128(lo, hi)).parity();
    }
};
#elif defined(SIMD_EMULATION)
template <>
struct SimdRegister<256> : public SimdRegisterEmulator<8>
{
    using SimdRegisterEmulator<8>::SimdRegisterEmulator;
    SimdRegister<256>(const SimdRegisterEmulator<8>& v) : SimdRegisterEmulator<8>(v) {}
};
#endif

#if SIMD_N_BITS>=512
template <>
struct SimdRegister<512>
{
    __m512i m_v;

    typedef SimdRegister<512> XV;

    SimdRegister() : m_v(_mm512_undefined_si512()) {}
    SimdRegister(const void* p) : m_v(_mm512_load_si512((const __m512i*)p)) {}
    SimdRegister(uint32_t v) : m_v(_mm512_set1_epi32(v)) {}
    SimdRegister(__m512i v) : m_v(v) {}

    template <bool A>
    void store(uint32_t* dst) { if (A) _mm512_store_si512((__m512i*)dst, m_v); else _mm512_storeu_si512((__m512i*)dst, m_v); }

    //    static FORCE_INLINE XV load(const void* p) { return _mm512_loadu_si512((const __m512i*) p); }

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return _mm512_and_si512(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return _mm512_xor_si512(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { return _mm512_or_si512(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n) { return _mm512_slli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n) { return _mm512_srli_epi32(a.m_v, n); }

    FORCE_INLINE XV ifOddValueElseZero(const XV& value) const
    {
        const __m512i lowestBit = _mm512_slli_epi32(m_v, 31); // move least significant bit to most significant bit
        const __mmask16 isOdd = _mm512_movepi32_mask(lowestBit);
        return _mm512_maskz_mov_epi32(isOdd, value.m_v);
    }

    static FORCE_INLINE XV zero() { return _mm512_setzero_si512(); }

    uint8_t parity() const
    {
        __m256i hi(_mm512_extracti64x4_epi64(m_v, 1));
        __m256i lo(_mm512_castsi512_si256(m_v));
        return SimdRegister<256>(_mm256_xor_si256(lo, hi)).parity();
    }
};
#elif defined(SIMD_EMULATION)
template <>
struct SimdRegister<512> : public SimdRegisterEmulator<16>
{
    using SimdRegisterEmulator<16>::SimdRegisterEmulator;
    SimdRegister<512>(const SimdRegisterEmulator<16>& v) : SimdRegisterEmulator<16>(v) {}
};
#endif




