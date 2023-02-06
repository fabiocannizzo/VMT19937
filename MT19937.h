#pragma once

/* 
   Converted to C++ from the original code downloaded from:
   mt19937ar.c
   http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/MT2002/emt19937ar.html
*/

#include <cstdint>
#include <cstddef>

// define FORCE_INLINE
#if defined(__GNUC__)

#elif defined(_MSC_VER) || defined(__INTEL_COMPILER)
#       define FORCE_INLINE __forceinline
#else
#       define FORCE_INLINE inline
#endif


#if defined(__GNUC__)
#   define FORCE_INLINE __attribute__((always_inline)) inline
#   if defined(__AVX__)
#       define VECLEN 8
#   elif defined(__SSE4_1__)
#       define VECLEN 4
#   endif
#elif defined(_MSC_VER)
#   define FORCE_INLINE __forceinline
#   if _M_IX86_FP==2 || defined(_M_X64) || defined(__SSE2__)
#       if defined(__AVX__)
#           define VECLEN 8
#       else
#           define VECLEN 4
#       endif
#   endif
#endif

#if VECLEN>0
#   include <immintrin.h>
#endif

template <size_t L>
struct V;

template <>
struct V<4>
{
    __m128i v;

    typedef V<4> XV;

    V() {}
    V(uint32_t v) : v(_mm_set1_epi32(v)){}
    V(int32_t v) : v(_mm_set1_epi32(v)) {}
    V(__m128i v) : v(v) {}

    template <bool Aligned>
    static FORCE_INLINE XV load(const uint32_t* p) { return Aligned ? _mm_load_si128((const __m128i*) p) : _mm_loadu_si128((const __m128i*) p); }

    template <bool Aligned>
    FORCE_INLINE void store(uint32_t* p) { if (Aligned) _mm_store_si128((__m128i*) p, v); else _mm_storeu_si128((__m128i*) p, v); }

    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { return _mm_or_si128(a.v, b.v); }
    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return _mm_and_si128(a.v, b.v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return _mm_xor_si128(a.v, b.v); }
    friend FORCE_INLINE XV operator==(const XV& a, const XV& b) { return _mm_cmpeq_epi32(a.v, b.v); }
    friend FORCE_INLINE XV operator>(const XV& a, const XV& b) { return _mm_cmpgt_epi32(a.v, b.v); }
    friend FORCE_INLINE XV operator<(const XV& a, const XV& b) { return _mm_cmplt_epi32(a.v, b.v); }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n) { return _mm_slli_epi32(a.v, n); }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n) { return _mm_srli_epi32(a.v, n); }

    // shift left the first 32-bit element of b int a: {a1, a2, a3, b0}
    static FORCE_INLINE XV shiftLeft(const XV& a, const XV& b) { return _mm_shuffle_epi32(_mm_blend_epi16(a.v, b.v, 0x3), 1 + (2 << 2) + (3 << 4)); }

    static FORCE_INLINE XV zero() { return _mm_setzero_si128(); }

    // returns 0xFFFFFFFF bitmask is the vale is odd
    FORCE_INLINE V<4> isOdd() const { return (*this & XV(1)) == XV(1); }
};

class MT19937
{
    template <typename T, size_t N, uint8_t ALIGN>
    struct AlignedArray
    {
        AlignedArray() : m_data(calcDataPtr()) {}

        FORCE_INLINE const T& operator[](size_t i) const { return m_data[i]; }
        FORCE_INLINE T& operator[](size_t i) { return m_data[i]; }

        FORCE_INLINE const T* begin() const { return m_data; }
        FORCE_INLINE T* begin() { return m_data; }

    private:
        FORCE_INLINE T* calcDataPtr()
        {
            uint8_t offset = static_cast<uint8_t>(ALIGN - reinterpret_cast<size_t>(&m_mem[0]) % ALIGN) % ALIGN;
            return  reinterpret_cast<T*>(m_mem + offset);
        }

    private:
        T* m_data;
        unsigned char m_mem[N * sizeof(T) + ALIGN - 1];
    };

    // Period parameters
    static const uint32_t N = 624;
    static const uint32_t M = 397;
    static const uint32_t MATRIX_A = 0x9908b0dfUL;   // constant vector a
    static const int32_t UPPER_MASK = 0x80000000UL; // most significant w-r bits
    static const int32_t LOWER_MASK = 0x7fffffffUL; // least significant r bits
    static const uint32_t c1 = 0x9d2c5680UL;
    static const uint32_t c2 = 0xefc60000UL;

    // both arrays are larger than necessary, because when we refill we read beyond the N-th element
    AlignedArray<uint32_t, N + VECLEN, 64> mt;  // the array for the state vector
    AlignedArray<uint32_t, N + VECLEN, 64> u32; // a cache of uniform discrete random numbers in the range [0,0xffffffff]
    uint32_t mti = N + 1;    // mti==N+1 means mt[N] is not initialized

    FORCE_INLINE uint32_t temper(uint32_t y)
    {
        // Tempering 
        y ^= (y >> 11);
        y ^= (y << 7) & c1;
        y ^= (y << 15) & c2;
        y ^= (y >> 18);
        return y;
    }

#if VECLEN>1
    struct Looper
    {
        const __m128i upper_mask = _mm_set1_epi32(UPPER_MASK);
        const __m128i lower_mask = _mm_set1_epi32(LOWER_MASK);
        const __m128i matrixa = _mm_set1_epi32(MATRIX_A);
        //const __m128i one = _mm_set1_epi32(1);
        const __m128i xc1 = _mm_set1_epi32(c1);
        const __m128i xc2 = _mm_set1_epi32(c2);

        template <size_t L>
        FORCE_INLINE V<L> temper(V<L> y)
        {
            // Tempering
            y = y ^ (y >> 11);
            y = y ^ (y << 7) & c1;
            y = y ^ (y << 15) & c2;
            y = y ^ (y >> 18);
            return y;
        }

        template <size_t N_ELEM, bool Align, size_t L>
        FORCE_INLINE V<L> body(const V<L>& curState, uint32_t* pmtCur, const uint32_t* pmtFar, uint32_t* pu32)
        {
            typedef V<L> XV;
            
            XV nextState = XV::template load<Align>(pmtCur + L);
            XV farState  = XV::template load<!Align>(pmtFar);

            XV cusStateP = XV::shiftLeft(curState, nextState);
            
            XV y = (curState & upper_mask) | (cusStateP & lower_mask);
            XV mag = y.isOdd() & matrixa;
            y = farState ^ (y >> 1) ^ mag;

            XV u32 = temper(y);
            
            if (N_ELEM == 4) {
                y.template store<Align>(pmtCur);
                u32.template store<Align>(pu32);
            }
            else {
                union { __m128i u128; uint32_t u32[4]; } yAux, u32Aux;
                yAux.u128 = y.v;
                u32Aux.u128 = u32.v;
                for (size_t i = 0; i < N_ELEM; ++i) {
                    pmtCur[i] = yAux.u32[i];
                    pu32[i] = u32Aux.u32[i];
                }
            }
            return nextState;
        }
        
    };
#endif

    void refill()
    {
        static uint32_t mag01[2] = { 0, MATRIX_A };

        if (mti == N + 1)   // if init_genrand() has not been called, 
            reinit(uint32_t(5489)); // a default initial seed is used 

#if VECLEN==4

        const size_t vecLen = 4;
        typedef V<vecLen> XV;

        Looper looper;

        {
            // process N-M elements

            const size_t nFull = (N - M) / VECLEN;

            uint32_t* pmt = mt.begin();
            uint32_t* pu32 = u32.begin();
            const uint32_t* pmt_end = pmt + nFull * VECLEN;

            XV curState = XV::load<true>(pmt);
            do {
                curState = looper.body<VECLEN, true>(curState, pmt, pmt + M, pu32);
                pmt += VECLEN;
                pu32 += VECLEN;
            } while (pmt != pmt_end);

            // in this iteration we read beyond the end of the state buffer
            // which is why we dimensioned the state buffer a little bit larger then necessary
            looper.body<(N - M) - nFull * VECLEN, true>(curState, pmt, pmt + M, pu32);
        }

        {
            // process M-1 elements

            const size_t nFull = (M - 1) / VECLEN;

            uint32_t* pmt = mt.begin() + (N - M);
            uint32_t* pu32 = u32.begin() + (N - M);
            const uint32_t* pmt_end = pmt + nFull * VECLEN;

            XV curState = XV::load<false>(pmt);
            do {
                curState = looper.body<VECLEN, false>(curState, pmt, pmt - (N-M), pu32);
                pmt += VECLEN;
                pu32 += VECLEN;
            } while (pmt != pmt_end);
        }

#else
        
        for (size_t kk = 0; kk < N - M; kk++) {
            uint32_t y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            y = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
            mt[kk] = y;
            u32[kk] = temper(y);
        }
        for (size_t kk = N - M; kk < N - 1; kk++) {
            uint32_t y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            y = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
            mt[kk] = y;
            u32[kk] = temper(y);
        }
#endif

        // process last element
        uint32_t y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        y = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];
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
            mt[mti] &= 0xffffffffUL;
            // for >32 bit machines
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
            mt[i] &= 0xffffffffUL; // for WORDSIZE > 32 machines
            i++; j++;
            if (i >= N) { mt[0] = mt[N - 1]; i = 1; }
            if (j >= n_seeds) j = 0;
        }
        for (k = N - 1; k; k--) {
            mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941UL))
                - i; // non linear 
            mt[i] &= 0xffffffffUL; // for WORDSIZE > 32 machines 
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

