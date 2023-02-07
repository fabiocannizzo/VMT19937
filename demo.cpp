#define _CRT_SECURE_NO_WARNINGS

#include "MT19937.h"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>

const uint32_t seedlength = 4;
const uint32_t seedinit[seedlength] = { 0x123, 0x234, 0x345, 0x456 };

const uint64_t nRandomTest = 500000;
const uint64_t nRandomPerf = 2000000000;

std::vector<uint32_t> benchmark(nRandomTest);

void printSome(const std::vector<uint32_t>& v)
{
    std::cout << "\n";
    for (size_t i = 0; i < 16; ++i)
        std::cout << std::setw(10) << v[i] << ((i+1) % 8 == 0 ? "\n" : " ");
    std::cout << "...\n";
    for (size_t i = 240; i < 240+16; ++i)
        std::cout << std::setw(10) << v[i] << ((i + 1) % 8 == 0 ? "\n" : " ");
    std::cout << "...\n";
    for (size_t i = 624 - 16; i < 624; ++i)
        std::cout << std::setw(10) << v[i] << ((i + 1) % 8 == 0 ? "\n" : " ");
    std::cout << "\n";
}

template <size_t VecLen>
double testPerformance()
{
    std::cout << "Generate " << nRandomPerf << " random numbers with SIMD length " << VecLen << " ... ";

    MT19937<VecLen> mt(seedinit, seedlength);

    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < nRandomPerf; ++i)
        mt.genrand_uint32();
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << "done in: " << std::fixed << std::setprecision(2) << elapsed_seconds.count() << "s" << std::endl;

    return 0;
}

template <size_t VecLen>
void testEquivalence()
{
    std::cout << "Testing equivalence of generators with SIMD length " << VecLen << " ... ";

    MT19937<VecLen> mt(seedinit, seedlength);

    for (size_t i = 0; i < nRandomTest; ++i) {
        uint32_t r2 = mt.genrand_uint32();
        if (benchmark[i] != r2) {
            std::cout << "FAILED!\n"
                      << "Difference found at index " << i << ": expected " << benchmark[i] << ", but got " << r2 << "\n";
            throw;
        }
    }

    std::cout << "SUCCESS!\n";
}

extern "C" unsigned long genrand_int32();
extern "C" void init_by_array(unsigned long init_key[], int key_length);


void generateBenchmark()
{
    unsigned long init[seedlength];
    for (size_t i = 0; i < seedlength; ++i)
        init[i] = seedinit[i];

    std::cout << "Generate random numbers with the original source code ... ";
    init_by_array(init, seedlength);
    for (size_t i = 0; i < nRandomTest; ++i)
        benchmark[i] = (uint32_t) genrand_int32();
    std::cout << "done!\n";
    printSome(benchmark);
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
    generateBenchmark();

    testEquivalence<1>();
    testEquivalence<4>();
#if MT19937_SIMD_VEC_LEN > 4
    //testEquivalence<8>();
#endif

    std::cout << "\n";

    originalPerformance();
    testPerformance<1>();
    testPerformance<4>();
#if MT19937_SIMD_VEC_LEN > 4
    //testPerformance<8>();
#endif

    return 0;
}
