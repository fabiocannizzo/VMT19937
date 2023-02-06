#define _CRT_SECURE_NO_WARNINGS

#include "MT19937.h"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <ctime>

template <size_t VecLen>
double testPerformance()
{
    const uint64_t nRandom = 2000000000;
    uint32_t init[4] = { 0x123, 0x234, 0x345, 0x456 }, length = 4;

    std::cout << "Generating " << nRandom << " discrete uniform 32-bit random numbers in [0,2^32) with SIMD length " << VecLen << "\n";

    MT19937<VecLen> mt(init, length);

    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < 3000000000; ++i)
        mt.genrand_uint32();
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << "Completed in: " << elapsed_seconds.count() << "s" << std::endl;

    return 0;
}

template <size_t VecLen1, size_t VecLen2>
void testEquivalence()
{
    std::cout << "Testing equivalence of generators with SIMD length " << VecLen1 << " and " << VecLen2 << " ... ";

    const uint64_t nRandom = 500000;
    uint32_t init[4] = { 0x123, 0x234, 0x345, 0x456 }, length = 4;

    MT19937<VecLen1> mt1(init, length);
    MT19937<VecLen2> mt2(init, length);

    for (size_t i = 0; i < nRandom; ++i) {
        uint32_t r1 = mt1.genrand_uint32();
        uint32_t r2 = mt2.genrand_uint32();
        if (r1 != r2) {
            std::cout << "FAILED!\n";
            throw;
        }
    }

    std::cout << "SUCCESS!\n";
}

int main()
{
    testEquivalence<1, 4>();
#if MT19937_SIMD_VEC_LEN > 4
    testEquivalence<1, 8>();
#endif

    testPerformance<1>();
    testPerformance<4>();
#if MT19937_SIMD_VEC_LEN > 4
    testPerformance<8>();
#endif

    return 0;
}
