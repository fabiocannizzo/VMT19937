#pragma once

/* 
   Converted to C++ from the original code downloaded from:
   mt19937ar.c
   http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/MT2002/emt19937ar.html
*/

#include <cstdint>
#include <cstddef>

class MT19937
{
    template <typename T, size_t N, uint8_t ALIGN>
    struct AlignedArray
    {
        AlignedArray() : m_data(calcDataPtr()) {}

        const T& operator[](size_t i) const { return m_data[i]; }
        T& operator[](size_t i) { return m_data[i]; }

        const T* begin() const { return m_data; }
        T* begin() { return m_data; }

        uint32_t size() const { return N; }

    private:
        T* calcDataPtr()
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
    static const uint32_t UPPER_MASK = 0x80000000UL; // most significant w-r bits
    static const uint32_t LOWER_MASK = 0x7fffffffUL; // least significant r bits

    AlignedArray<uint32_t, N, 64> mt;  // the array for the state vector
    AlignedArray<uint32_t, N, 64> u32; // a cache of uniform discrete random numbers in the range [0,0xffffffff]
    uint32_t mti = N + 1;    // mti==N+1 means mt[N] is not initialized

    uint32_t temper(uint32_t y)
    {
        const uint32_t c1 = 0x9d2c5680UL;
        const uint32_t c2 = 0xefc60000UL;
        // Tempering 
        y ^= (y >> 11);
        y ^= (y << 7) & c1;
        y ^= (y << 15) & c2;
        y ^= (y >> 18);
        return y;
    }

    void refill()
    {
        static uint32_t mag01[2] = { 0x0UL, MATRIX_A };

        if (mti == N + 1)   // if init_genrand() has not been called, 
            reinit(uint32_t(5489)); // a default initial seed is used 

        size_t kk;

        for (kk = 0; kk < N - M; kk++) {
            uint32_t y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            y = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
            mt[kk] = y;
            u32[kk] = temper(y);
        }
        for (; kk < N - 1; kk++) {
            uint32_t y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            y = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
            mt[kk] = y;
            u32[kk] = temper(y);
        }
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

