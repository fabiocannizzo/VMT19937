#define _CRT_SECURE_NO_WARNINGS

#include "MT19937-SIMD.h"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <tuple>
#include <set>

const uint32_t seedlength = 4;
const uint32_t seedinit[seedlength] = { 0x123, 0x234, 0x345, 0x456 };

const uint64_t nRandomPerf = uint64_t(624) * 16 * 500000;

extern "C" unsigned long genrand_int32();
extern "C" void init_by_array(unsigned long init_key[], int key_length);

struct Result
{
    Result(bool orig, size_t nb, size_t blk, size_t id, double dt)
        : original(orig), nBits(nb), blkSize(blk), runId(id), time(dt) {}
    bool original;
    size_t nBits;
    size_t blkSize;
    mutable size_t runId;
    mutable double time;
    bool operator<(const Result& rhs) const {
        return std::tuple(!original, nBits, blkSize) < std::tuple(!rhs.original, rhs.nBits, rhs.blkSize);
    }
};

template <size_t VecLen, size_t BlkSize = 1>
Result testPerformance(size_t runId)
{
    std::cout << "Generate " << nRandomPerf << " random numbers with SIMD length " << VecLen;
    if (BlkSize > 1)
        std::cout << " in blocks of " << BlkSize;
    std::cout << " ... ";

    std::vector<uint32_t> dst(BlkSize + 64 / sizeof(uint32_t));
    uint32_t* aligneddst = (uint32_t*)((intptr_t)dst.data() + (64 - ((intptr_t)dst.data() % 64)));

    MT19937SIMD<VecLen> mt(seedinit, seedlength, NULL);

    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < nRandomPerf / BlkSize; ++i)
        switch (BlkSize) {
            case 1: aligneddst[0] = mt.genrand_uint32(); break;
            case 16: mt.genrand_uint32_blk64(aligneddst); break;
            case (624 * (VecLen / 32)):  mt.genrand_uint32_stateBlk(aligneddst); break;
            default: throw std::invalid_argument("not implemented");
        };
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();

    std::cout << "done in: " << std::fixed << std::setprecision(2) << nSeconds << "s" << std::endl;

    return Result(runId, false, VecLen, BlkSize, nSeconds);
}


Result originalPerformance(size_t runId)
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

    return Result(runId, true, 32, 1, nSeconds);
}

void usage()
{
    std::cerr
        << "Invalid command line arguments\n"
        << "Example:\n"
        << "jump [-n nRepeats]\n"
        << "  nRepeats defaults to 1\n"
        ;
    std::exit(-1);
}

int main(int argc, const char** argv)
{
    size_t nRepeat = 1;
    // parse command line arguments
    for (int i = 1; i < argc; i += 2) {
        string key(argv[i]);
        const char* value = argv[i + 1];
        if (key == "-n") {
            nRepeat = atoi(value);
        }
        else
            usage();
    }
    std::cout << "nRepeat = " << nRepeat << "\n";

    std::vector<Result> res;
    for (size_t i = 0; i < nRepeat; ++i) {

        res.push_back(originalPerformance(i));
        res.push_back(testPerformance<32>(i));
        //testPerformance<64>();
        res.push_back(testPerformance<128>(i));
        res.push_back(testPerformance<128, 16>(i));
        res.push_back(testPerformance<128, 624 * (128 / 32)>(i));
#if SIMD_N_BITS >= 256
        res.push_back(testPerformance<256>(i));
        res.push_back(testPerformance<256, 16>(i));
        res.push_back(testPerformance<256, 624 * (256 / 32)>(i));
#endif
#if SIMD_N_BITS >= 512
        res.push_back(testPerformance<512>(i));
        res.push_back(testPerformance<512, 16>(i));
        res.push_back(testPerformance<512, 624 * (512 / 32)>(i));
#endif
    }

    std::set<Result> avgresult;
    for (const auto& r : res) {
        auto [iter, ins] = avgresult.insert(r);
        if (ins)
            iter->runId = 1;
        else {
            ++iter->runId;
            iter->time += r.time;
        }
    }
    std::cout << std::setw(8) << std::right << "method"
        << std::setw(8) << std::right << "n-bits"
        << std::setw(8) << std::right << "blksize"
        << std::setw(8) << std::right << "time"
        << "\n";
    for (auto& r : avgresult) {
        std::cout << std::setw(8) << std::right << (r.original ? "orig" : "simd")
            << std::setw(8) << std::right << r.nBits
            << std::setw(8) << std::right << r.blkSize
            << std::setw(8) << std::right << std::fixed << std::setprecision(2) << r.time / nRepeat
            << "\n";
    }

    return 0;
}