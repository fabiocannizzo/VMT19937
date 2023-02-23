#define _CRT_SECURE_NO_WARNINGS

#include "MT19937-SIMD.h"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>

const uint32_t seedlength = 4;
const uint32_t seedinit[seedlength] = { 0x123, 0x234, 0x345, 0x456 };

const uint64_t nRandomPerf = uint64_t(624) * 16 * 500000;

extern "C" unsigned long genrand_int32();
extern "C" void init_by_array(unsigned long init_key[], int key_length);

template <size_t VecLen>
double testPerformance()
{
    std::cout << "Generate " << nRandomPerf << " random numbers with SIMD length " << VecLen << " ... ";

    MT19937SIMD<VecLen> mt(seedinit, seedlength, NULL);

    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < nRandomPerf; ++i)
        mt.genrand_uint32();
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << "done in: " << std::fixed << std::setprecision(2) << elapsed_seconds.count() << "s" << std::endl;

    return 0;
}

template <size_t VecLen>
double testVecPerformance()
{
    const size_t BlkSize = 64 / 4;// 624 * VecLen / 32;
    std::cout << "Generate " << nRandomPerf << " random numbers with SIMD length " << VecLen << " in blocks of " << BlkSize << " ... ";

    MT19937SIMD<VecLen> mt(seedinit, seedlength, NULL);

    //std::vector<uint32_t> dst(BlkSize);
    alignas(64) uint32_t dst[BlkSize];
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < nRandomPerf / BlkSize; ++i)
        mt.genrand_uint32_blk64(dst);
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << "done in: " << std::fixed << std::setprecision(2) << elapsed_seconds.count() << "s" << std::endl;

    return 0;
}

void originalPerformance()
{
    unsigned long init[seedlength];
    for (size_t i = 0; i < seedlength; ++i)
        init[i] = seedinit[i];

    init_by_array(init, seedlength);
    std::cout << "Generate " << nRandomPerf << " random numbers with original code ... ";
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < nRandomPerf; ++i)
        genrand_int32();
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();
    std::cout << "done in: " << std::fixed << std::setprecision(2) << nSeconds << "s\n";
}

int main()
{
    originalPerformance();
    testPerformance<32>();
    //testPerformance<64>();
    //testPerformance<128>();
    testVecPerformance<128>();
#if SIMD_N_BITS >= 256
    testPerformance<256>();
#endif
#if SIMD_N_BITS >= 512
    testPerformance<512>();
#endif

    return 0;
}
