#pragma once

/* 
   Converted to C++ from the original code downloaded from:
   mt19937ar.c
   http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/MT2002/emt19937ar.html
   and vectorized using SIMD
*/

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

// define MT19937_SIMD_VEC_LEN
#if defined(__GNUC__) || defined(__INTEL_COMPILER)
#   if defined(__AVX__)
#       define MT19937_SIMD_VEC_LEN 8
#   elif defined(__SSE4_1__)
#       define MT19937_SIMD_VEC_LEN 4
#   endif
#elif defined(_MSC_VER)
#   if _M_IX86_FP==2 || defined(_M_X64) || defined(__SSE2__)
#       if defined(__AVX__)
#           define MT19937_SIMD_VEC_LEN 8
#       else
#           define MT19937_SIMD_VEC_LEN 4
#       endif
#   endif
#else
#    define MT19937_SIMD_VEC_LEN 1
#endif

//#undef MT19937_SIMD_VEC_LEN
//#define MT19937_SIMD_VEC_LEN 1
//#undef FORCE_INLINE
//#define FORCE_INLINE __declspec(noinline)

#if MT19937_SIMD_VEC_LEN>0
#   include <immintrin.h>
#else
#error d
#endif

template <size_t L>
struct V;

#if MT19937_SIMD_VEC_LEN>4
template <>
struct V<8>
{
    __m256i m_v;

    typedef V<8> XV;

    V() {}
    V(uint32_t v) : m_v(_mm256_set1_epi32(v)) {}
    V(int32_t v) : m_v(_mm256_set1_epi32(v)) {}
    V(__m256i v) : m_v(v) {}

    template <bool Aligned>
    static FORCE_INLINE XV load(const uint32_t* p) { return Aligned ? _mm256_load_si256((const __m256i*) p) : _mm256_loadu_si256((const __m256i*) p); }

    template <bool Aligned>
    FORCE_INLINE void store(uint32_t* p) { if (Aligned) _mm256_store_si256((__m256i*) p, m_v); else _mm256_storeu_si256((__m256i*) p, m_v); }

    template <size_t N>
    FORCE_INLINE void storeN(uint32_t* p) {
        _mm256_maskstore_epi32((int32_t*)p, _mm256_set_epi32(0, N > 6, N > 5, N > 4, N > 3, N > 2, N > 1, 1), m_v);
    }

    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { return _mm256_or_si256(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return _mm256_and_si256(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return _mm256_xor_si256(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator==(const XV& a, const XV& b) { return _mm256_cmpeq_epi32(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator>(const XV& a, const XV& b) { return _mm256_cmpgt_epi32(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator<(const XV& a, const XV& b) { return _mm256_cmpgt_epi32(b.m_v, a.m_v); }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n) { return _mm256_slli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n) { return _mm256_srli_epi32(a.m_v, n); }

    // shift left the first 32-bit element of b int a: {a1, a2, a3, a4, a5, a6, a7, b0}
    static FORCE_INLINE XV shiftLeft(const XV& a, const XV& b) { return _mm256_permutevar8x32_epi32(_mm256_blend_epi32(a.m_v, b.m_v, 0x1), _mm256_set_epi32(0, 7, 6, 5, 4, 3, 2, 1)); }

    // returns value if v is odd, zero otherwise
    FORCE_INLINE XV ifOddValueElseZero(const XV& value) const { return ((*this & XV(1)) == XV(1)) & value; }
};
#endif


template <>
struct V<4>
{
    __m128i m_v;

    typedef V<4> XV;

    V() {}
    V(uint32_t v) : m_v(_mm_set1_epi32(v)){}
    V(int32_t v) : m_v(_mm_set1_epi32(v)) {}
    V(__m128i v) : m_v(v) {}
#if MT19937_SIMD_VEC_LEN>4
    V(const V<8>& v) : m_v(_mm256_castsi256_si128(v.m_v)) {}
#endif

    template <bool Aligned>
    static FORCE_INLINE XV load(const uint32_t* p) { return Aligned ? _mm_load_si128((const __m128i*) p) : _mm_loadu_si128((const __m128i*) p); }

    template <bool Aligned>
    FORCE_INLINE void store(uint32_t* p) { if (Aligned) _mm_store_si128((__m128i*) p, m_v); else _mm_storeu_si128((__m128i*) p, m_v); }

    template <size_t N>
    FORCE_INLINE void storeN(uint32_t* p) {
        const char y1 = N > 1 ? -1 : 0;  // store element 1 ?
        const char y2 = N > 2 ? -1 : 0;  // store element 2 ?
        _mm_maskmoveu_si128(m_v, _mm_set_epi8(0, 0, 0, 0, y2, y2, y2, y2, y1, y1, y1, y1, -1, -1, -1, -1), (char*)p);
    }

    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { return _mm_or_si128(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return _mm_and_si128(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return _mm_xor_si128(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator==(const XV& a, const XV& b) { return _mm_cmpeq_epi32(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator>(const XV& a, const XV& b) { return _mm_cmpgt_epi32(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator<(const XV& a, const XV& b) { return _mm_cmplt_epi32(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n) { return _mm_slli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n) { return _mm_srli_epi32(a.m_v, n); }

    // shift left the first 32-bit element of b int a: {a1, a2, a3, b0}
    static FORCE_INLINE XV shiftLeft(const XV& a, const XV& b) { return _mm_shuffle_epi32(_mm_blend_epi16(a.m_v, b.m_v, 0x3), 1 + (2 << 2) + (3 << 4)); }

    // returns value if v is odd, zero otherwise
    FORCE_INLINE XV ifOddValueElseZero(const XV& value) const { return ((*this & XV(1)) == XV(1)) & value; }
};

template <>
struct V<1>
{
    uint32_t m_v;

    typedef V<1> XV;

    V() {}
    V(uint32_t v) : m_v(v) {}
    V(int32_t v) : m_v(v) {}

    template <bool Aligned>
    static FORCE_INLINE XV load(const uint32_t* p) { return *p; }

    template <bool Aligned>
    FORCE_INLINE void store(uint32_t* p) { *p = m_v; }

    template <size_t N>
    FORCE_INLINE void storeN(uint32_t* p) { throw; /* we should never get here */ }

    friend FORCE_INLINE XV operator|(XV a, XV b) { return a.m_v | b.m_v; }
    friend FORCE_INLINE XV operator&(XV a, XV b) { return a.m_v & b.m_v; }
    friend FORCE_INLINE XV operator^(XV a, XV b) { return a.m_v ^ b.m_v; }
    friend FORCE_INLINE XV operator<<(XV a, uint32_t n) { return a.m_v << n; }
    friend FORCE_INLINE XV operator>>(XV a, uint32_t n) { return a.m_v >> n; }

    // return b
    static FORCE_INLINE XV shiftLeft(XV a, XV b) { return b; }

    // returns value if v is odd, zero otherwise
    FORCE_INLINE XV ifOddValueElseZero(XV value) const { return m_v & 0x1 ? value : 0; }
};

template <size_t _VecLen = MT19937_SIMD_VEC_LEN>
class MT19937
{
    const static size_t VecLen = _VecLen;
    typedef V<VecLen> XV;

    template <typename VecTy, size_t N, uint8_t ALIGN>
    struct AlignedArray
    {
        AlignedArray() : m_data(calcDataPtr()) {}

        FORCE_INLINE const VecTy& operator[](size_t i) const { return m_data[i]; }
        FORCE_INLINE VecTy& operator[](size_t i) { return m_data[i]; }

        FORCE_INLINE const VecTy* begin() const { return m_data; }
        FORCE_INLINE VecTy* begin() { return m_data; }

    private:
        FORCE_INLINE VecTy* calcDataPtr()
        {
            uint8_t offset = static_cast<uint8_t>(ALIGN - reinterpret_cast<size_t>(&m_mem[0]) % ALIGN) % ALIGN;
            return  reinterpret_cast<VecTy*>(m_mem + offset);
        }

    private:
        VecTy* m_data;
        unsigned char m_mem[N * sizeof(VecTy) + ALIGN - 1];
    };

    // Period parameters
    static const uint32_t N = 624;
    static const uint32_t M = 397;
    static const uint32_t s_matrixA = 0x9908b0dfUL;   // constant vector a
    static const int32_t s_upperMask = 0x80000000UL; // most significant w-r bits
    static const int32_t s_lowerMask = 0x7fffffffUL; // least significant r bits
    static const uint32_t s_temperMask1 = 0x9d2c5680UL;
    static const uint32_t s_temperMask2 = 0xefc60000UL;

    // both arrays are larger than necessary, because when we refill we read beyond the N-th element
    AlignedArray<uint32_t, N + VecLen, 64> mt;  // the array for the state vector
    AlignedArray<uint32_t, N + VecLen, 64> u32; // a cache of uniform discrete random numbers in the range [0,0xffffffff]
    uint32_t mti = N + 1;    // mti==N+1 means mt[N] is not initialized

    template <typename VecTy>
    FORCE_INLINE VecTy temper(VecTy y)
    {
        // Tempering
        y = y ^ (y >> 11);
        y = y ^ (y << 7) & VecTy(s_temperMask1);
        y = y ^ (y << 15) & VecTy(s_temperMask2);
        y = y ^ (y >> 18);
        return y;
    }

    template <size_t N_ELEM, bool Align, size_t L>
    FORCE_INLINE V<L> body(const V<L>& curState, uint32_t* pmtCur, const uint32_t* pmtFar, uint32_t* pu32)
    {
        typedef V<L> VecTy;

        VecTy nextState(VecTy::template load<Align>(pmtCur + L));
        VecTy farState(VecTy::template load<!Align>(pmtFar));

        VecTy cusStateP(VecTy::shiftLeft(curState, nextState));
            
        VecTy y = (curState & s_upperMask) | (cusStateP & s_lowerMask);
        VecTy mag = y.ifOddValueElseZero(VecTy(s_matrixA));
        y = farState ^ (y >> 1) ^ mag;

        VecTy u32 = temper(y);
            
        if (N_ELEM == L) {
            y.template store<Align>(pmtCur);
            u32.template store<Align>(pu32);
        }
        else {
            y.template storeN<N_ELEM>(pmtCur);
            u32.template storeN<N_ELEM>(pu32);
        }
        return nextState;
    }

    void refill()
    {
        if (mti == N + 1)   // if init_genrand() has not been called, 
            reinit(uint32_t(5489)); // a default initial seed is used 

        uint32_t* pmt;
        uint32_t* pu32;
        const uint32_t* pmt_end;

        // ************************************
        // process first N-M elements
        //

        const size_t nFull1 = (N - M) / VecLen;
        const size_t nLeft1 = (N - M) % VecLen;

        pmt = mt.begin();
        pu32 = u32.begin();
        pmt_end = pmt + nFull1 * VecLen;

        XV curState = XV::template load<true>(pmt);
        do {
            curState = body<VecLen, true>(curState, pmt, pmt + M, pu32);
            pmt += VecLen;
            pu32 += VecLen;
        } while (pmt != pmt_end);

        // in this iteration we read beyond the end of the state buffer
        // which is why we dimensioned the state buffer a little bit larger then necessary
        if (nLeft1)
            body<nLeft1, true, (nLeft1 <= 1 ? 1 : nLeft1 <= 4 ? 4 :8)>(curState, pmt, pmt + M, pu32);

        // process remaining M-1 elements

        const size_t nFull2 = (M - 1) / VecLen;
        const size_t nLeft2 = (M - 1) % VecLen;

        pmt = mt.begin() + (N - M);
        pu32 = u32.begin() + (N - M);
        pmt_end = pmt + nFull2 * VecLen;

        curState = XV::template load<false>(pmt);
        do {
            curState = body<VecLen, false>(curState, pmt, pmt - (N-M), pu32);
            pmt += VecLen;
            pu32 += VecLen;
        } while (pmt != pmt_end);

        if (nLeft2)
            body<nLeft2, true, (nLeft2 <= 1 ? 1 : nLeft2 <= 4 ? 4 : 8)>(curState, pmt, pmt - (N - M), pu32);

        // ****************************
        // process last element
        //

        //static uint32_t mag01[2] = { 0, s_matrixA };
        uint32_t y = (mt[N - 1] & s_upperMask) | (mt[0] & s_lowerMask);
        y = mt[M - 1] ^ (y >> 1) ^ (y & 0x1UL ? s_matrixA : 0);
        mt[N - 1] = y;
        u32[N - 1] = temper(y);

        mti = 0;
    }

public:
    // constructors
    MT19937() : mti(N + 1) {}

    MT19937(uint32_t seed)
    {
        reinit(seed);
    }

    MT19937(uint32_t seeds[], uint32_t n_seeds)
    {
        reinit(seeds, n_seeds);
    }

    // initializes mt[N] with a seed
    void reinit(uint32_t s)
    {
        mt[0] = s & 0xffffffffUL;
        for (mti = 1; mti < N; mti++) {
            mt[mti] =
                (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
            // See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
            // In the previous versions, MSBs of the seed affect
            // only MSBs of the array mt[].
            // 2002/01/09 modified by Makoto Matsumoto
            // mt[mti] &= 0xffffffffUL;  // for >32 bit machines
        }
    }

    // initialize by an array with array-length
    // init_key is the array for initializing keys
    // key_length is its length
    // slight change for C++, 2004/2/26
    void reinit(uint32_t seeds[], uint32_t n_seeds)
    {
        uint32_t i, j, k;
        reinit(uint32_t(19650218));
        i = 1; j = 0;
        k = (N > n_seeds ? N : n_seeds);
        for (; k; k--) {
            mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525UL))
                + seeds[j] + j; // non linear
            //mt[i] &= 0xffffffffUL; // for WORDSIZE > 32 machines
            i++; j++;
            if (i >= N) { mt[0] = mt[N - 1]; i = 1; }
            if (j >= n_seeds) j = 0;
        }
        for (k = N - 1; k; k--) {
            mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941UL))
                - i; // non linear 
            //mt[i] &= 0xffffffffUL; // for WORDSIZE > 32 machines 
            i++;
            if (i >= N) { mt[0] = mt[N - 1]; i = 1; }
        }

        mt[0] = 0x80000000UL; // MSB is 1; assuring non-zero initial array 
    }

    // generates a random number on [0,0xffffffff] interval 
    uint32_t genrand_uint32()
    {
        if (mti >= N) // generate N words at one time 
            refill();
        uint32_t y = u32[mti++];
        return y;
    }

    // generates a random number on [0,0x7fffffff]-int32_terval 
    uint32_t genrand_uint31(void)
    {
        return (long)(genrand_uint32() >> 1);
    }

    // generates a random number on [0,1]-real-int32_terval 
    double genrand_real1(void)
    {
        return genrand_uint32() * (1.0 / 4294967295.0);
        // divided by 2^32-1 
    }

    // generates a random number on [0,1)-real-int32_terval 
    double genrand_real2(void)
    {
        return genrand_uint32() * (1.0 / 4294967296.0);
        // divided by 2^32 
    }

    // generates a random number on (0,1)-real-int32_terval 
    double genrand_real3(void)
    {
        return (((double)genrand_uint32()) + 0.5) * (1.0 / 4294967296.0);
        // divided by 2^32 
    }

    // generates a random number on [0,1) with 53-bit resolution
    double genrand_res53(void)
    {
        uint32_t a = genrand_uint32() >> 5, b = genrand_uint32() >> 6;
        return(a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
    }
    // These real versions are due to Isaku Wada, 2002/01/09 added 
};

