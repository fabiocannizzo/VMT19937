#include "TestUtils.h"
#include "cli_args.h"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <tuple>
#include <set>
#include <numeric>
#include <algorithm>
#include <string>
#include <cmath>
#include <random>

using namespace std;
using namespace xvmt;
using namespace xvmt::details;

#ifdef NO_MKL
constexpr bool c_skip_mkl = true;
#else
#   if __has_include(<mkl.h>)
#       define MKL_AVAIL 1
        constexpr bool c_skip_mkl = false;
#   else
#       if !defined(__arm__) && !defined(__aarch64__)
#           pragma message("MKL not found, disabling MKL tests")
#       endif
        constexpr bool c_skip_mkl = true;
#   endif
#endif

#ifdef MKL_AVAIL
#   include <mkl.h>
#endif

#ifdef NO_ORIG
constexpr bool c_skip_original = true;
#else
constexpr bool c_skip_original = false;
#endif

#if !c_skip_original
#   define SFMT_MEXP 19937
#   include "../SFMT-src-1.5.1/SFMT.h"
#endif

#ifdef NO_STL
constexpr bool c_skip_stl = true;
#else
constexpr bool c_skip_stl = false;
#endif

#ifdef NO_XMT
constexpr bool c_skip_xmt = true;
#else
constexpr bool c_skip_xmt = false;
#endif

#ifdef NO_SFMT
constexpr bool c_skip_sfmt = true;
#else
constexpr bool c_skip_sfmt = false;
#endif

#ifdef NO_VMT
constexpr bool c_skip_vmt = true;
#else
constexpr bool c_skip_vmt = false;
#endif

#ifdef NO_VMT_SMALL
constexpr bool c_skip_vmt_small = true;
#else
constexpr bool c_skip_vmt_small = false;
#endif

#ifdef NO_GEN_32
constexpr bool c_skip_gen_32 = true;
#else
constexpr bool c_skip_gen_32 = false;
#endif

#ifdef NO_GEN_64
constexpr bool c_skip_gen_64 = true;
#else
constexpr bool c_skip_gen_64 = false;
#endif

#ifdef NO_GEN_SFMT
constexpr bool c_skip_gen_sfmt = true;
#else
constexpr bool c_skip_gen_sfmt = false;
#endif

// global variables
bool g_skip_mkl = c_skip_mkl;
bool g_skip_original = c_skip_original;
bool g_skip_stl = c_skip_stl;
bool g_skip_xmt = c_skip_xmt;
bool g_skip_sfmt = c_skip_sfmt;
bool g_skip_vmt = c_skip_vmt;
bool g_skip_vmt_small = c_skip_vmt_small;
bool g_skip_gen_32 = c_skip_gen_32;
bool g_skip_gen_64 = c_skip_gen_64;
bool g_skip_gen_sfmt = c_skip_gen_sfmt;

bool g_skip_qry1 = false;
bool g_skip_qry16 = false;
bool g_skip_qryN = false;
string g_dir = "dat";
size_t g_nRepeat = 1;
double g_stDev = 0.0;

// this might be changed via cli arguments
size_t g_nRandom = size_t(624) * 32 * 800;

constexpr uint32_t s_seedlength = 4;
constexpr uint32_t s_seedinit[s_seedlength] = { 0x123, 0x234, 0x345, 0x456 };


extern "C" unsigned long genrand_int32();
extern "C" void init_by_array(unsigned long init_key[], int key_length);
extern "C" unsigned long long genrand64_int64();
extern "C" void init_by_array64(unsigned long long init_key[], unsigned long long key_length);

enum GenMode {orig, sfmt, mkl_mt, mkl_sfmt, xmt32, vmt, vsfmt, xsfmt, stl_mt, xmt64, stl_mt64, orig64, vmt64};

const char* modename[] = {"ORIG-MT19937", "ORIG-SFMT19937", "MKL-MT19937", "MKL-SFMT19937", "X-MT19937", "V-MT19937", "V-SFMT19937", "X-SFMT19937", "STL-MT19937", "X-MT19937-64", "STL-MT19937-64", "ORIG-MT19937-64", "V-MT19937-64"};

const size_t anySize[] = {/* 1, 4, 16, 64, 256, 624, 1024, 4096,*/ 16384 };

template <GenMode G>
struct GenTraits;

// for maximum period, we should select the file based on the number of states
// but these periods are so large anyway that who do not care!
const auto pmt = std::make_unique<MT19937Matrix<32>>(g_dir + "/mt32/F19933.bits");
// for maximum period, we should select the file based on the number of states
// but these periods are so large anyway that who do not care!
const auto psfmt = std::make_unique<SFMT19937Matrix>(g_dir + "/sfmt/F19935.bits");
const auto pvmt64 = std::make_unique<MT19937Matrix<64>>(g_dir + "/mt64/F19933.bits");

// use the same destination memory in all tests to avoid spurious difference in test results due to memory layout
// sized for the largest 64-bit anySize block (each element = 8 bytes)
AlignedVector<uint32_t, 64> aligneddst(anySize[sizeof(anySize)/sizeof(anySize[0])-1] * 2);

template <>
struct GenTraits<vmt>
{
    static const GenMode mode = vmt;
    static const MT19937Matrix<32>* jumpMatrix() { return pmt.get(); }

    template <size_t RegBitLen, QryMode QM, size_t HwRegBitLen>
    using gen_t = VMT19937<RegBitLen, QM == QM_Block16, BitLenToIsa<HwRegBitLen>::isa>;
};

template <>
struct GenTraits<xmt32>
{
    static const GenMode mode = xmt32;
    static const MT19937Matrix<32>* jumpMatrix() { return nullptr; }

    template <size_t RegBitLen, QryMode QM, size_t HwRegBitLen, std::enable_if_t<RegBitLen == HwRegBitLen, int> = 0>
    using gen_t = XMT19937<BitLenToIsa<HwRegBitLen>::isa, QM == QM_Block16>;
};

template <>
struct GenTraits<vsfmt>
{
    static const GenMode mode = vsfmt;
    static const SFMT19937Matrix* jumpMatrix() { return psfmt.get(); }

    template <size_t RegBitLen, QryMode QM, size_t HwRegBitLen>
    using gen_t = VSFMT19937<RegBitLen, QM == QM_Block16, BitLenToIsa<HwRegBitLen>::isa>;
};

// X-SFMT19937: single-state SFMT (VRegBitLen == 128 == s_stateWordBits, so s_nStates == 1)
template <>
struct GenTraits<xsfmt>
{
    static const GenMode mode = xsfmt;
    static const SFMT19937Matrix* jumpMatrix() { return nullptr; }

    template <size_t RegBitLen, QryMode QM, size_t HwRegBitLen, std::enable_if_t<RegBitLen == HwRegBitLen, int> = 0>
    using gen_t = VSFMT19937<RegBitLen, QM == QM_Block16, BitLenToIsa<HwRegBitLen>::isa>;
};

template <>
struct GenTraits<xmt64>
{
    static const GenMode mode = xmt64;
    static const MT19937Matrix<32>* jumpMatrix() { return nullptr; }

    template <size_t RegBitLen, QryMode QM, size_t HwRegBitLen, std::enable_if_t<RegBitLen == HwRegBitLen, int> = 0>
    using gen_t = XMT19937_64<BitLenToIsa<HwRegBitLen>::isa, QM == QM_Block16>;
};

template <>
struct GenTraits<vmt64>
{
    static const GenMode mode = vmt64;
    static const MT19937Matrix<64>* jumpMatrix() { return pvmt64.get(); }

    template <size_t RegBitLen, QryMode QM, size_t HwRegBitLen>
    using gen_t = VMT19937_64<RegBitLen, QM == QM_Block16, BitLenToIsa<HwRegBitLen>::isa>;
};

const size_t s_messageSpacing[] = { 15, 9, 8, 8, 8, 12 };

struct Results
{
    Results(GenMode _mode, size_t _nb, size_t _ib, size_t _blk, QryMode _qryMode)
        : mode(_mode), nBits(_nb), nBitsHw(_ib), blkSize(_blk), qryMode(_qryMode), nStates(1), mi(0), ma(0), avg(0), stdev(0)
    {}
    GenMode mode;
    size_t nBits, nBitsHw;
    size_t blkSize;
    QryMode qryMode;
    size_t nStates;

    mutable std::vector<double> singleRuns;
    mutable double mi, ma, avg, stdev;

    void print() const
    {
        size_t i = 0;
        std::cout
            << std::setw(s_messageSpacing[i++]) << modename[mode]
            << std::setw(s_messageSpacing[i++]) << nBits
            << std::setw(s_messageSpacing[i++]) << nStates
            << std::setw(s_messageSpacing[i++]) << nBitsHw
            << std::setw(s_messageSpacing[i++]) << blkSize
            << std::setw(s_messageSpacing[i++]) << queryModeName(qryMode)
            << " ... ";
    }

    bool operator<(const Results& b) const
    {
        return std::tuple(mode, nBits, nBitsHw, blkSize, (int)qryMode) < std::tuple(b.mode, b.nBits, b.nBitsHw, b.blkSize, (int)b.qryMode);
    }
};

struct TableCompare
{
    bool operator()(const Results& a, const Results& b) const
    {
        return std::tuple(a.nBitsHw, a.blkSize, a.avg) < std::tuple(b.nBitsHw, b.blkSize, b.avg);
    }
};


std::set<Results> results;

void done(double nSeconds)
{
    std::cout << "done in: " << std::setw(8) << std::fixed << std::setprecision(2) << nSeconds << "s\n";
}

struct Results;
bool alreadyHaveEnoughIter(const Results& key);

// add results and update statistics
void addResult(const Results& key, double seconds)
{
    const Results& r = *results.insert(key).first;
    auto& v = r.singleRuns;
    v.push_back(seconds);
    double f = v.front();
    double n = (double) v.size();
    r.mi = f;
    r.ma = f;
    double s = f, s2 = f*f;
    for (size_t i = 1; i < n; ++i) {
        f = v[i];
        r.mi = std::min(r.mi, f);
        r.ma = std::max(r.ma, f);
        s += f;
        s2 += f * f;
    }
    r.avg = s / n;
    if (n > 1)
        r.stdev = std::sqrt((s2 - s * r.avg) / (n - 1));
}

NO_INLINE
bool alreadyHaveEnoughIter(const Results& key)
{
    auto iter = results.find(key);
    bool notEnough = (iter == results.end())
        || (iter->singleRuns.size() < g_nRepeat)
        || (std::abs(iter->stdev / iter->avg) > g_stDev / 100.0);
    //if (notEnough)
    //    std::cout << "running test...\n";
    return !notEnough;
}

#ifndef NO_ORIG
void mtOrigPerformance()
{
    Results key(orig, 32, 32, 1, QM_Scalar);

    key.print();

    if (alreadyHaveEnoughIter(key)) {
        std::cout << "skip\n";
        return;
    }

    std::vector<uint32_t> dst(1);

    unsigned long init[s_seedlength];
    for (size_t i = 0; i < s_seedlength; ++i)
        init[i] = s_seedinit[i];

    init_by_array(init, s_seedlength);

    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < g_nRandom; ++i) {
        volatile uint32_t val = genrand_int32();
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();
    done(nSeconds);

    addResult(key, nSeconds);
}

void mtOrig64Performance()
{
    Results key(orig64, 64, 64, 1, QM_Scalar);

    key.print();

    if (alreadyHaveEnoughIter(key)) {
        std::cout << "skip\n";
        return;
    }

    unsigned long long init[s_seedlength];
    for (size_t i = 0; i < s_seedlength; ++i)
        init[i] = s_seedinit[i];

    init_by_array64(init, s_seedlength);

    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < g_nRandom; ++i) {
        volatile uint64_t val = genrand64_int64();
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();
    done(nSeconds);

    addResult(key, nSeconds);
}
#endif

#ifndef NO_STL
void stlMtPerformance()
{
    Results key(stl_mt, 32, 32, 1, QM_Scalar);
    key.print();
    if (alreadyHaveEnoughIter(key)) {
        std::cout << "skip\n";
        return;
    }

    std::mt19937 gen(5489);

    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < g_nRandom; ++i) {
        volatile uint32_t val = gen();
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();
    done(nSeconds);

    addResult(key, nSeconds);
}

void stlMtPerformanceVectorial()
{
    Results key(stl_mt, 32, 32, 1, QM_Any);
    key.print();
    if (alreadyHaveEnoughIter(key)) {
        std::cout << "skip\n";
        return;
    }

    std::mt19937 gen(5489);

    auto start = std::chrono::system_clock::now();
    uint32_t sum = 0;
    for (size_t i = 0; i < g_nRandom; ++i)
        sum += gen();
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();
    done(nSeconds);

    addResult(key, nSeconds);
}

void stlMt64Performance()
{
    Results key(stl_mt64, 64, 64, 1, QM_Scalar);
    key.print();
    if (alreadyHaveEnoughIter(key)) {
        std::cout << "skip\n";
        return;
    }

    std::mt19937_64 gen(5489ULL);

    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < g_nRandom; ++i) {
        volatile uint64_t val = gen();
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();
    done(nSeconds);

    addResult(key, nSeconds);
}

void stlMt64PerformanceVectorial()
{
    Results key(stl_mt64, 64, 64, 1, QM_Any);
    key.print();
    if (alreadyHaveEnoughIter(key)) {
        std::cout << "skip\n";
        return;
    }

    std::mt19937_64 gen(5489ULL);

    auto start = std::chrono::system_clock::now();
    uint64_t sum = 0;
    for (size_t i = 0; i < g_nRandom; ++i)
        sum += gen();
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();
    done(nSeconds);

    addResult(key, nSeconds);
}
#endif

#ifndef NO_ORIG
template <bool ScalarQry>
void sfmtOrigPerformance(size_t BlkSize)
{
    Results key(sfmt, 128, 128, BlkSize, ScalarQry ? QM_Scalar : QM_Any);
    key.print();
    if (alreadyHaveEnoughIter(key)) {
        std::cout << "skip\n";
        return;
    }

    MYASSERT((BlkSize == 1) || (BlkSize % 4 == 0 && BlkSize >= SFMT_N32), "BlkSize must be a multiple of 4 and >=156*128");
    MYASSERT((g_nRandom % BlkSize) == 0, "nRandom must be a multiple of BlkSize");

    sfmt_t sfmtgen;
    sfmt_init_gen_rand(&sfmtgen,12345);


    auto start = std::chrono::system_clock::now();
    for (size_t i = 0, n = g_nRandom / BlkSize; i < n; ++i) {
        if constexpr (ScalarQry) {
            volatile uint32_t val = sfmt_genrand_uint32(&sfmtgen);
        }
        else
            sfmt_fill_array32(&sfmtgen, aligneddst.data(), (int) BlkSize);
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();
    done(nSeconds);

    addResult(key, nSeconds);
}
#endif

#ifdef MKL_AVAIL

template <int N>
struct MKLTraits
    ;
template <>
struct MKLTraits<VSL_BRNG_MT19937>
{
    static const GenMode s_mode = mkl_mt;
    static const size_t s_wordSize = 32;
};

template <>
struct MKLTraits<VSL_BRNG_SFMT19937>
{
    static const GenMode s_mode = mkl_sfmt;
    static const size_t s_wordSize = 128;
};

void mklPerformance(MKL_INT GenCode, MKL_INT BlkSize)
{
    GenMode mode;
    size_t wordSize;
    switch (GenCode) {
        case VSL_BRNG_MT19937:
            mode = mkl_mt;
            wordSize = 32;
            break;
        case VSL_BRNG_SFMT19937:
            mode = mkl_sfmt;
            wordSize = 128;
            break;
        default: THROW("how did we get here?");
    }

    MYASSERT(g_nRandom % BlkSize == 0, "incorrect count");

    Results key(mode, wordSize, SIMD_N_BITS, BlkSize, QM_Any);
    key.print();

    if (alreadyHaveEnoughIter(key)) {
        std::cout << "skip\n";
        return;
    }

    VSLStreamStatePtr stream;
    vslNewStream(&stream, GenCode, 5489);

    auto start = std::chrono::system_clock::now();
    for (size_t i = 0, n = g_nRandom / BlkSize; i < n; ++i)
        viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, (int)BlkSize, aligneddst.data());
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();
    done(nSeconds);

    // Deleting the stream
    vslDeleteStream(&stream);

    addResult(key, nSeconds);
}
#endif

template <GenMode Mode, size_t L, size_t I, QryMode QM>
void vRandGenPerformance5(size_t blkSize)
{
    MYASSERT(((blkSize > 0) && ((g_nRandom % blkSize) == 0)), "invalid blkSize " << blkSize);

    using Gen = typename GenTraits<Mode>::template gen_t<L, QM, I>;

    Results key(Mode, L, I, blkSize, QM);
    key.nStates = Gen::s_nStates;

    key.print();

    if (alreadyHaveEnoughIter(key)) {
        std::cout << "skip\n";
        return;
    }

    using output_word_t = typename Gen::output_word_t;

    const typename Gen::matrix_t *jumpMatrixPtr = nullptr;
    if constexpr (Gen::s_nStates > 1)
        jumpMatrixPtr = GenTraits<Mode>::jumpMatrix();

    auto makeGen = [&]() -> Gen {
        if constexpr (sizeof(output_word_t) == 8)
            return Gen(uint64_t(5489), 0, nullptr, jumpMatrixPtr);
        else
            return Gen(s_seedinit, s_seedlength, 0, nullptr, jumpMatrixPtr);
    };
    Gen mt = makeGen();

    auto start = std::chrono::system_clock::now();

    if constexpr (QM == QM_Scalar) {
        if constexpr (sizeof(output_word_t) == 8) {
            for (size_t i = 0, n = g_nRandom / blkSize; i < n; ++i) {
                volatile uint64_t val = mt.genrand_uint64();
            }
        }
        else {
            for (size_t i = 0, n = g_nRandom / blkSize; i < n; ++i) {
                volatile uint32_t val = mt.genrand_uint32();
            }
        }
    }
    else {
        for (size_t i = 0, n = g_nRandom / blkSize; i < n; ++i) {
            if constexpr (QM == QM_Block16) {
                if constexpr (sizeof(output_word_t) == 8)
                    mt.genrand_word_blk(reinterpret_cast<output_word_t*>(aligneddst.data()));
                else
                    mt.genrand_uint32_blk16(aligneddst.data());
            }
            else if constexpr (QM == QM_Any) {
                if constexpr (sizeof(output_word_t) == 8)
                    mt.genrand_word_anySize(reinterpret_cast<output_word_t*>(aligneddst.data()), blkSize);
                else
                    mt.genrand_uint32_anySize(aligneddst.data(), blkSize);
            }
            else
                NOT_IMPLEMENTED;
        }
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();

    done(nSeconds);

    addResult(key, nSeconds);
}


template <GenMode Mode, size_t L, size_t I, QryMode QM>
void vRandGenPerformance4()
{
    using output_word_t = typename GenTraits<Mode>::template gen_t<L, QM, I>::output_word_t;
    if constexpr (QM == QM_Any) {
        if (!g_skip_qryN)
            for (auto sz : anySize)
                vRandGenPerformance5<Mode, L, I, QM>(sz);
    }
    else if constexpr (QM == QM_Block16) {
        if (!g_skip_qry16) {
            constexpr size_t blkSz = (sizeof(output_word_t) == 8) ? 8 : 16;
            vRandGenPerformance5<Mode, L, I, QM>(blkSz);
        }
    }
    else if constexpr (QM == QM_Scalar) {
        if (!g_skip_qry1)
            vRandGenPerformance5<Mode, L, I, QM>(1);
    }
}

template <GenMode Mode, size_t L, size_t I, QryMode...QMs>
void vRandGenPerformance2()
{
    constexpr size_t M = std::min<size_t>(L, SIMD_N_BITS);
    if constexpr (I <= M && (Mode != xmt32 || I == L) && (Mode != xmt64 || I == L) && (Mode != xsfmt || I == L)) {
        if constexpr (!c_skip_vmt_small || I == L) {
            if (!g_skip_vmt_small || I == L) {
                (vRandGenPerformance4<Mode, L, I, QMs>(), ...);
            }
        }
    }
}

template <GenMode Mode, size_t L, size_t...Is>
void vRandGenPerformance1()
{
    (vRandGenPerformance2<Mode, L, Is, QM_Scalar, QM_Block16, QM_Any>(), ...);
}

template <GenMode Mode, size_t...Ls>
void vRandGenPerformance0()
{
    (vRandGenPerformance1<Mode, Ls, /*32,*/ 128, 256, 512>(), ...);
}

void syntax()
{
    std::cout
        << "Usage: perf [OPTIONS]\n"
        << "Options:\n"
        << "  -h, --help            Show this help message\n"
        << "  -n=<nRepeats>         Number of performance test iterations (default 1)\n"
        << "  -s=<nRndScaler>       Multiply nRandom by this scaler\n"
        << "  --slow                Increase the number of random numbers generated by 1000 times\n"
        << "  --dir=<datpath>       Folder where to find jump matrix files (default: dat)\n"
        << "  --stdev=<val>         Stop iterating when stdev/avg < val% (nRepeats is treated as minimum)\n"
        << "  -wait                 Wait for debugger attachment\n"
        << "Exclusion flags:\n"
        << "  --no-original         Skip original MT19937/SFMT19937 tests\n"
        << "  --no-stl              Skip std::mt19937 tests\n"
        << "  --no-mkl              Skip MKL tests\n"
        << "  --no-xmt              Skip XMT tests\n"
        << "  --no-sfmt             Skip SFMT-related tests\n"
        << "  --no-vmt              Skip VMT tests\n"
        << "  --no-vmt-small        Skip VMT/VSFMT tests where HwRegBitLen < VRegBitLen\n"
        << "  --no-gen-32           Skip all 32-bit generators (MT19937)\n"
        << "  --no-gen-64           Skip all 64-bit generators (MT19937-64)\n"
        << "  --no-gen-sfmt         Skip all SFMT generators\n"
        << "  --no-qry1             Skip scalar query tests\n"
        << "  --no-qry16            Skip block-16 query tests\n"
        << "  --no-qryN             Skip vectorial query tests\n";
}

// CPU initialization
void initCpu()
{
    try {
        detectCpuInfo().print();
        setCpuAffinity(1);
        setPriorityHigh();
    }
    catch (const std::exception& ex) {
        std::cerr << "Error during CPU initialization: " << ex.what() << "\n";
        std::exit(-1);
    }
    catch (...) {
        std::cerr << "Unknown error during CPU initialization\n";
        std::exit(-1);
    }
}

#ifdef MKL_AVAIL
void mklPerformance(MKL_INT GenCode, MKL_INT BlkSize);
#endif

// init MKL dispatching
void initMKL()
{
#ifdef MKL_AVAIL
    if (!g_skip_mkl) {
        try {
#  if (SIMD_N_BITS==128)
            std::cout << "Force MKL dispatching to SSE2\n";
            MYASSERT(mkl_enable_instructions(MKL_ENABLE_SSE4_2), "SSE2 not supported");
#  elif (SIMD_N_BITS==256)
            std::cout << "set MKL dispacthing to AVX2\n";
            MYASSERT(mkl_enable_instructions(MKL_ENABLE_AVX2), "AVX2 not supported");
#  elif (SIMD_N_BITS==512)
            std::cout << "set MKL dispacthing to AVX512\n";
            MYASSERT(mkl_enable_instructions(MKL_ENABLE_AVX512), "AVX512 not supported");
#  else
            NOT_IMPLEMENTED;
#  endif
        }
        catch (const std::exception& ex) {
            std::cerr << "Error during MKL initialization: " << ex.what() << "\n";
            std::exit(-1);
        }
        catch (...) {
            std::cerr << "Unknown error during MKL initialization\n";
            std::exit(-1);
        }
    }
#endif
}

void parseCliArgs(ArgMap& args)
{
    // parse command line arguments
    try {
        if (consumeArg(args, "h") || consumeArg(args, "help")) {
            syntax();
            std::exit(0);
        }

        consumeArg(args, "n", false, g_nRepeat);

        if (consumeArg(args, "no-mkl")) g_skip_mkl = true;
        if (consumeArg(args, "no-original")) g_skip_original = true;
        if (consumeArg(args, "no-stl")) g_skip_stl = true;
        if (consumeArg(args, "no-vmt")) g_skip_vmt = true;
        if (consumeArg(args, "no-vmt-small")) g_skip_vmt_small = true;
        if (consumeArg(args, "no-sfmt")) g_skip_sfmt = true;
        if (consumeArg(args, "no-xmt")) g_skip_xmt = true;
        if (consumeArg(args, "no-gen-32")) g_skip_gen_32 = true;
        if (consumeArg(args, "no-gen-64")) g_skip_gen_64 = true;
        if (consumeArg(args, "no-gen-sfmt")) g_skip_gen_sfmt = true;

        if (consumeArg(args, "no-qry1")) g_skip_qry1 = true;
        if (consumeArg(args, "no-qry16")) g_skip_qry16 = true;
        if (consumeArg(args, "no-qryN")) g_skip_qryN = true;

        size_t rndScaler = 1;
        if (consumeArg(args, "s", false, rndScaler)) g_nRandom *= rndScaler;

        consumeArg(args, "stdev", false, g_stDev);
        consumeArg(args, "dir", false, g_dir);
        if (consumeArg(args, "slow")) g_nRandom *= 1000;

        if (!args.empty()) {
            syntax();
            std::exit(-1);
        }

        if (g_stDev > 0.0)
            g_nRepeat = std::max<size_t>(g_nRepeat, 2);
        MYASSERT(g_nRepeat > 0, "nRepeat must be positive");
    }
    catch (const std::exception& ex) {
        std::cerr << "Error parsing CLI arguments: " << ex.what() << "\n";
        syntax();
        std::exit(-1);
    }
}

int main(int argc, const char** argv)
{
    ArgMap args = parseArgs(argc, argv);
    waitForDebugger(args);

    // detect CPU and set affinity and priority
    initCpu();

    // init MKL dispatching
    initMKL();

    // parse command line arguments
    parseCliArgs(args);

    // print some test information
    std::cout << "Target hardware SIMD register size (bits): " << SIMD_N_BITS << "\n";
    std::cout << "nRepeat = " << g_nRepeat << "\n";
    std::cout << "stDev = " << g_stDev << "\n";
    std::cout << "nRandom = " << g_nRandom << "\n";
    std::cout << (!g_skip_mkl ? "including" : "skipping") << " MKL tests\n";
    std::cout << (!g_skip_original ? "including" : "skipping") << " original implementation tests\n";
    std::cout << (!g_skip_stl ? "including" : "skipping") << " STL tests\n";
    std::cout << (!g_skip_vmt ? "including" : "skipping") << " VMT tests\n";
    std::cout << (!g_skip_xmt ? "including" : "skipping") << " XMT tests\n";
    std::cout << (!g_skip_sfmt ? "including" : "skipping") << " SFMT tests\n";

    // run all tests
    try {

        auto nResults = []() {
            return std::accumulate(results.begin(), results.end(), size_t(0),
                [](size_t s, const Results& r) { return s + r.singleRuns.size(); });
        };

        for (size_t i = 0; i < g_nRepeat || g_stDev > 0; ++i) {

            size_t nResBefore = nResults();

            std::cout << "Iteration: " << std::setw(2) << i + 1 << ": "
                      << "generating " << g_nRandom << " 32-bits random numbers\n";
            {
                size_t m = 0;
                std::cout
                    << std::setw(s_messageSpacing[m++]) << "Generator"
                    << std::setw(s_messageSpacing[m++]) << "VReg"
                    << std::setw(s_messageSpacing[m++]) << "nStates"
                    << std::setw(s_messageSpacing[m++]) << "HwReg"
                    << std::setw(s_messageSpacing[m++]) << "BlkSize"
                    << std::setw(s_messageSpacing[m++]) << "QueryMode"
                    << "\n";
            }

            // Original & STL MT32
            if constexpr (!c_skip_gen_32) {
#ifndef NO_ORIG
                if (!g_skip_original && !g_skip_gen_32) {
                    mtOrigPerformance();
                }
#endif
                if constexpr (!c_skip_stl) {
                    if (!g_skip_stl && !g_skip_gen_32) {
                        if (!g_skip_qry1) stlMtPerformance();
                        if (!g_skip_qryN) stlMtPerformanceVectorial();
                    }
                }
            }

            // Original SFMT
            if constexpr (!c_skip_original && !c_skip_sfmt && !c_skip_gen_sfmt) {
#ifndef NO_ORIG
                if (!g_skip_original && !g_skip_sfmt && !g_skip_gen_sfmt) {
                    if (!g_skip_qry1)
                        sfmtOrigPerformance<true>(1);
                    if (!g_skip_qryN) {
                        for (auto sz : anySize)
                            if (sz >= 624)
                                sfmtOrigPerformance<false>(sz);
                    }
                }
#endif
            }

            // MKL MT32 & SFMT
            if constexpr (!c_skip_mkl) {
#ifdef MKL_AVAIL
                if (!g_skip_mkl) {
                    if constexpr (!c_skip_gen_32) {
                        if (!g_skip_gen_32) {
                            if (!g_skip_qry1)
                                mklPerformance(VSL_BRNG_MT19937, 1);
                            if (!g_skip_qryN) {
                                for (auto sz : anySize)
                                    mklPerformance(VSL_BRNG_MT19937, (MKL_INT)sz);
                            }
                        }
                    }
                    if constexpr (!c_skip_sfmt && !c_skip_gen_sfmt) {
                        if (!g_skip_sfmt && !g_skip_gen_sfmt) {
                            if (!g_skip_qry1)
                                mklPerformance(VSL_BRNG_SFMT19937, 1);
                            if (!g_skip_qryN) {
                                for (auto sz : anySize)
                                    mklPerformance(VSL_BRNG_SFMT19937, (MKL_INT)sz);
                            }
                        }
                    }
                }
#endif
            }

            // VMT MT32
            if constexpr (!c_skip_vmt && !c_skip_gen_32) {
                if (!g_skip_vmt && !g_skip_gen_32) {
                    vRandGenPerformance0<vmt, SIMD_N_BITS>();
                }
            }

            // XMT MT32
            if constexpr (!c_skip_xmt && !c_skip_gen_32) {
                if (!g_skip_xmt && !g_skip_gen_32) {
                    vRandGenPerformance0<xmt32, SIMD_N_BITS>();
                }
            }

            // SFMT (VMT & XMT)
            if constexpr (!c_skip_sfmt && !c_skip_gen_sfmt) {
                if (!g_skip_sfmt && !g_skip_gen_sfmt) {
                    // xsfmt is single-state VSFMT, treated as "small" or just SFMT
                    vRandGenPerformance0<xsfmt, 128>();
                    if constexpr (!c_skip_vmt) {
                        if (!g_skip_vmt) {
#if SIMD_N_BITS >= 256
                            vRandGenPerformance0<vsfmt, SIMD_N_BITS>();
#endif
                        }
                    }
                }
            }

            // 64-bit MT
            if constexpr (!c_skip_gen_64) {
                if (!g_skip_gen_64) {
#ifndef NO_ORIG
                    if (!g_skip_original) mtOrig64Performance();
#endif
                    if constexpr (!c_skip_stl) {
                        if (!g_skip_stl) {
                            if (!g_skip_qry1) stlMt64Performance();
                            if (!g_skip_qryN) stlMt64PerformanceVectorial();
                        }
                    }
                    if constexpr (!c_skip_xmt) {
                        if (!g_skip_xmt) vRandGenPerformance0<xmt64, SIMD_N_BITS>();
                    }
                    if constexpr (!c_skip_vmt) {
                        if (!g_skip_vmt) vRandGenPerformance0<vmt64, SIMD_N_BITS>();
                    }
                }
            }

            size_t nResAfter = nResults();
            if (nResAfter == nResBefore)
                break; // no new results added
        }

        std::set<Results, TableCompare> sortedResults(results.begin(), results.end());

        const size_t spacing[] = { 20, 8, 8, 8, 8, 10, 6, 8, 8, 8, 8, 11, 12 };
        size_t s = 0;
        std::cout << "\n"
            << std::setw(spacing[s++]) << std::right << "prng"
            << std::setw(spacing[s++]) << std::right << "VReg"
            << std::setw(spacing[s++]) << std::right << "nStates"
            << std::setw(spacing[s++]) << std::right << "HwReg"
            << std::setw(spacing[s++]) << std::right << "blksize"
            << std::setw(spacing[s++]) << std::right << "qrymode"
            << std::setw(spacing[s++]) << std::right << "nruns"
            << std::setw(spacing[s++]) << std::right << "tmin"
            << std::setw(spacing[s++]) << std::right << "tmax"
            << std::setw(spacing[s++]) << std::right << "tavg"
            << std::setw(spacing[s++]) << std::right << "tdev"
            << std::setw(1 + spacing[s++]) << std::right << "tdev/tavg"
            << std::setw(spacing[s++]) << std::right << "throughput"
            << "\n";
        for (auto& r : sortedResults) {
            s = 0;
            std::cout << std::setw(spacing[s++]) << std::right << modename[r.mode]
                << std::setw(spacing[s++]) << std::right << r.nBits
                << std::setw(spacing[s++]) << std::right << r.nStates
                << std::setw(spacing[s++]) << std::right << r.nBitsHw
                << std::setw(spacing[s++]) << std::right << r.blkSize
                << std::setw(spacing[s++]) << std::right << queryModeName(r.qryMode)
                << std::setw(spacing[s++]) << std::right << r.singleRuns.size()
                << std::setw(spacing[s++]) << std::right << std::fixed << std::setprecision(3) << r.mi
                << std::setw(spacing[s++]) << std::right << std::fixed << std::setprecision(3) << r.ma
                << std::setw(spacing[s++]) << std::right << std::fixed << std::setprecision(3) << r.avg
                << std::setw(spacing[s++]) << std::right << std::fixed << std::setprecision(3) << r.stdev
                << std::setw(spacing[s++]) << std::right << std::fixed << std::setprecision(3) << r.stdev / r.avg * 100 << "%"
                << std::setw(spacing[s++]) << std::right << std::fixed << std::setprecision(1) << g_nRandom / 1.0e6 / r.avg
                << "\n";
        }
    }
    catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return -1;
    }
    catch (...) {
        std::cerr << "Unknown error\n";
        return -1;
    }

    return 0;
}
