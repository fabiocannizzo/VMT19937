#include "TestUtils.h"

#include "../SFMT-src-1.5.1/SFMT.h"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <tuple>
#include <map>
#include <numeric>
#include <algorithm>
#include <cmath>

using namespace std;

#define TEST_MKL 1
#define TEST_VMT 1
#define TEST_ORIG 1

#if TEST_MKL==1
#   if __has_include(<mkl.h>)
#       include <mkl.h>
#   else
#       pragma message "MKL not found, disabling MKL tests"
#       undef TEST_MKL
#       define TEST_MKL 0
#   endif
#endif

#if TEST_ORIG==1
#   include "../SFMT-src-1.5.1/SFMT.h"
#endif

// global variables
bool g_testMkl = TEST_MKL;
bool g_testOriginal = true;
bool g_testVMT = true;
size_t g_nRepeat = 1;
std::string dir = "dat";

// this might be changed via cli arguments
size_t g_nRandom = size_t(624) * 32 * 800;

constexpr uint32_t s_seedlength = 4;
constexpr uint32_t s_seedinit[s_seedlength] = { 0x123, 0x234, 0x345, 0x456 };


extern "C" unsigned long genrand_int32();
extern "C" void init_by_array(unsigned long init_key[], int key_length);

enum Mode {orig, sfmt, mkl_mt, mkl_sfmt, vmt, vsfmt};

const char* modename[] = {"ORIG-MT19937", "ORIG-SFMT19937", "MKL-MT19937", "MKL-SFMT19937", "V-MT19937", "V-SFMT19937" };

const size_t anySize[] = { 1, 4, 16, 64, 256, 624, 1024, 4096, 16384 };

template <typename G>
struct GenTraits;

// for maximum period, we should select the file based on the number of states
// but these periods are so large anyway that who do not care!
std::unique_ptr<Details::VMT19937Base<32, 32>::matrix_t> pmt(new Details::VMT19937Base<32, 32>::matrix_t(dir + "/mt/F19933.bits"));
// for maximum period, we should select the file based on the number of states
// but these periods are so large anyway that who do not care!
std::unique_ptr<Details::VSFMT19937Base<128, 32>::matrix_t> psfmt(new Details::VSFMT19937Base<128, 32>::matrix_t(dir + "/sfmt/F19935.bits"));


template <size_t RegBitLen, VRandGenQueryMode QueryMode, size_t RegBitLenHw>
struct GenTraits<VMT19937<RegBitLen, QueryMode, RegBitLenHw>>
{
    static const Mode mode = vmt;
    static const auto* matrix() { return pmt.get(); }
};

template <size_t RegBitLen, VRandGenQueryMode QueryMode, size_t RegBitLenHw>
struct GenTraits<VSFMT19937<RegBitLen, QueryMode, RegBitLenHw>>
{
    static const Mode mode = vsfmt;
    static const auto* matrix() { return psfmt.get(); }
};

const size_t s_messageSpacing[] = { 15, 9, 8, 8, 12 };

struct ResultKey
{
    ResultKey(Mode _mode, size_t _nb, size_t _ib, size_t _blk, VRandGenQueryMode _qryMode)
        : mode(_mode), nBits(_nb), nImplBits(_ib), blkSize(_blk), qryMode(_qryMode) {}
    Mode mode;
    size_t nBits, nImplBits;
    size_t blkSize;
    VRandGenQueryMode qryMode;
    void print() const
    {
        size_t i = 0;
        std::cout
            << std::setw(s_messageSpacing[i++]) << modename[mode]
            << std::setw(s_messageSpacing[i++]) << nBits
            << std::setw(s_messageSpacing[i++]) << nImplBits
            << std::setw(s_messageSpacing[i++]) << blkSize
            << std::setw(s_messageSpacing[i++]) << queryModeName(qryMode)
            << " ... ";
    }
    bool operator<(const ResultKey& rhs) const
    {
        return std::tuple(mode, nBits, nImplBits, blkSize, (int) qryMode) < std::tuple(rhs.mode, rhs.nBits, rhs.nImplBits, rhs.blkSize, (int) rhs.qryMode);
    }
};

struct ResultValues
{
    ResultValues() : mi(0), ma(0), avg(0), stdev(0) {}
    std::vector<double> singleRuns;
    double mi, ma, avg, stdev;
};

std::map<ResultKey, ResultValues> results;

void done(double nSeconds)
{
    std::cout << "done in: " << std::setw(8) << std::fixed << std::setprecision(2) << nSeconds << "s\n";
}

void addResult(const ResultKey& key, double seconds)
{
    auto& res = results.insert({ key, ResultValues{} }).first->second;
    auto& v = res.singleRuns;
    v.push_back(seconds);
    double f = v.front();
    double n = (double) v.size();
    res.mi = f;
    res.ma = f;
    double s = f, s2 = f*f;
    for (size_t i = 1; i < n; ++i) {
        f = v[i];
        res.mi = std::min(res.mi, f);
        res.ma = std::max(res.ma, f);
        s += f;
        s2 += f * f;
    }
    res.avg = s / n;
    if (n > 1)
        res.stdev = std::sqrt((s2 - s * res.avg) / (n - 1));
}

bool alreadyHaveEnoughIter(const ResultKey& key)
{
    auto iter = results.find(key);
    bool notEnough = (iter == results.end())
        || (iter->second.singleRuns.size() <= 1)
        || (std::abs(iter->second.stdev / iter->second.avg) > 0.01);
    return !notEnough;
}

#if TEST_ORIG==1
void mtOrigPerformance()
{
    ResultKey key(orig, 32, 32, 1, QM_Scalar);

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
    for (size_t i = 0; i < g_nRandom; ++i)
        dst[0] = genrand_int32();
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();
    done(nSeconds);

    addResult(key, nSeconds);
}
#endif

template <bool ScalarQry>
void sfmtOrigPerformance(size_t BlkSize)
{
    ResultKey key(sfmt, 128, 128, BlkSize, ScalarQry ? QM_Scalar : QM_Any);
    key.print();
    if (alreadyHaveEnoughIter(key)) {
        std::cout << "skip\n";
        return;
    }

    MYASSERT((BlkSize == 1) || (BlkSize % 4 == 0 && BlkSize >= SFMT_N32), "BlkSize must be a multiple of 4 and >=156*128");
    MYASSERT((g_nRandom % BlkSize) == 0, "nRandom must be a multiple of BlkSize");
    AlignedVector<uint32_t, 64> aligneddst(BlkSize);

    sfmt_t sfmtgen;
    sfmt_init_gen_rand(&sfmtgen,12345);


    auto start = std::chrono::system_clock::now();
    for (size_t i = 0, n = g_nRandom / BlkSize; i < n; ++i) {
        if constexpr (ScalarQry)
            aligneddst[0] = sfmt_genrand_uint32(&sfmtgen);
        else
            sfmt_fill_array32(&sfmtgen, aligneddst.data(), (int) BlkSize);
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();
    done(nSeconds);

    addResult(key, nSeconds);
}

#if TEST_MKL==1

template <int N>
struct MKLTraits
    ;
template <>
struct MKLTraits<VSL_BRNG_MT19937>
{
    static const Mode s_mode = mkl_mt;
    static const size_t s_wordSize = 32;
};

template <>
struct MKLTraits<VSL_BRNG_SFMT19937>
{
    static const Mode s_mode = mkl_sfmt;
    static const size_t s_wordSize = 128;
};

void mklPerformance(MKL_INT GenCode, MKL_INT BlkSize)
{
    Mode mode;
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

    ResultKey key(mode, wordSize, SIMD_N_BITS, BlkSize, QM_Any);
    key.print();

    if (alreadyHaveEnoughIter(key)) {
        std::cout << "skip\n";
        return;
    }

    AlignedVector<uint32_t, 64> aligneddst(BlkSize);

    VSLStreamStatePtr stream;
    vslNewStream(&stream, GenCode, 5489);

    auto start = std::chrono::system_clock::now();
    for (size_t i = 0, n = g_nRandom / BlkSize; i < n; ++i)
        viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, BlkSize, aligneddst.data());
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();
    done(nSeconds);

    // Deleting the stream
    vslDeleteStream(&stream);

    addResult(key, nSeconds);
}
#endif

template <typename Gen>
void vRandGenPerformance4(size_t blkSize)
{
    const Mode mode = GenTraits<Gen>::mode;

    MYASSERT(((blkSize > 0) && ((g_nRandom % blkSize) == 0)), "invalid blkSize " << blkSize);

    ResultKey key(mode, Gen::s_regLenBits, Gen::s_regLenBitsHw, blkSize, Gen::s_queryMode);

    key.print();

    if (alreadyHaveEnoughIter(key)) {
        std::cout << "skip\n";
        return;
    }

    AlignedVector<uint32_t, 64> aligneddst(blkSize);

    const Gen::matrix_t *jumpMatrixPtr = nullptr;
    if constexpr (Gen::s_nStates > 1)
        jumpMatrixPtr = GenTraits<Gen>::matrix();
    Gen mt(s_seedinit, s_seedlength, 0, nullptr, jumpMatrixPtr);

    auto start = std::chrono::system_clock::now();

    for (size_t i = 0, n = g_nRandom / blkSize; i < n; ++i)
        if constexpr (Gen::s_queryMode == QM_Scalar)
            aligneddst[0] = mt.genrand_uint32();
        else if constexpr (Gen::s_queryMode == QM_Block16)
            mt.genrand_uint32_blk16(aligneddst.data());
        else if constexpr (Gen::s_queryMode == QM_StateSize)
            mt.genrand_uint32_stateBlk(aligneddst.data());
        else if constexpr (Gen::s_queryMode == QM_Any)
            mt.genrand_uint32_anySize(aligneddst.data(), blkSize);
        else
            NOT_IMPLEMENTED;

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();

    done(nSeconds);

    addResult(key, nSeconds);
}


template <typename Gen>
void vRandGenPerformance3()
{
    size_t blkSize;
    switch (Gen::s_queryMode) {
        case QM_Scalar: blkSize = 1; break;
        case QM_Block16: blkSize = 16; break;
        case QM_Any: blkSize = 0; break;
        case QM_StateSize: blkSize = Gen::s_n32InFullState; break;
        default: THROW("how did we get here?");
    }
    bool fst = true;
    for (auto sz : anySize) {
        if (Gen::s_queryMode != QM_Any && !fst)
            break;
        fst = false;
        vRandGenPerformance4<Gen>(Gen::s_queryMode != QM_Any?  blkSize: sz);
    }
}

template <template <size_t, VRandGenQueryMode, size_t> class Gen, size_t L, size_t I, VRandGenQueryMode...QMs>
void vRandGenPerformance2()
{
    if constexpr (I <= std::min<size_t>(L, SIMD_N_BITS))
        (vRandGenPerformance3<Gen<L, QMs, I>>(), ...);
}

template <template <size_t, VRandGenQueryMode, size_t> class Gen, size_t L, size_t...Is>
void vRandGenPerformance1()
{
    (vRandGenPerformance2<Gen, L, Is, QM_Scalar, QM_Block16, QM_StateSize, QM_Any>(), ...);
}

template <template <size_t, VRandGenQueryMode, size_t> class Gen, size_t...Ls>
void vRandGenPerformance0()
{
    (vRandGenPerformance1<Gen, Ls, /*32,*/ 128, 256, 512>(), ...);
}

void usage()
{
    std::cerr
        << "Invalid command line arguments\n"
        << "Example:\n"
        << "perf [-n nRepeats] [-s slow] [-m minRepeats]\n"
        << "  nRepeats defaults to 1\n"
        << "  slow must be 0 or 1, defaults to 0\n"
        << "  minRepeats defaults to nRepeats\n"
        ;
}

void syntax()
{
    std::cout
        << "perf [-n nRepeats] [--slow] [--no-mkl] [--no-original] [--no-vmt]\n"
        << "  nRepeats: number of performance test iterations (default 1)\n"
        << "  --no-mkl: skip MKL tests\n"
        << "  --no-original: skip original implementations tests\n"
        << "  --slow: increase the number of random numbers generated by 1000 times\n"
        ;
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

// init MKL dispatching
void initMKL()
{
#if TEST_MKL==1
    if (g_testMkl) {
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

void parseCliArgs(int argc, const char** argv)
{
    // parse command line arguments
    g_nRepeat = 1;
    try {
        for (int i = 1; i < argc; ++i) {
            string key(argv[i]);
            if (key == "-n") {
            	MYASSERT(++i < argc, "-n must be followed by a number");
                const char* value = argv[i];
                g_nRepeat = stoul(string(value));
            }
            else if (key == "--no-mkl")
                g_testMkl = false;
            else if (key == "--no-original")
                g_testOriginal = false;
            else if (key == "--no-vmt")
                g_testVMT = false;
            else if (key == "--slow")
                g_nRandom *= 1000;
            else if (key == "--dir") {
            	MYASSERT(++i < argc, "--dir must be followed by a path");
                dir = argv[i];
            }
            else {
                usage();
                std::exit(-1);
            }
        }
        MYASSERT(g_nRepeat > 0, "nRepeat must be positive");
    }
    catch (const std::exception& ex) {
        std::cerr << "Error parsing command line arguments: " << ex.what() << "\n";
        usage();
        std::exit(-1);
    }
    catch (...) {
        usage();
        std::exit(-1);
    }
}

int main(int argc, const char** argv)
{
    // detect CPU and set affinity and priority
    initCpu();

    // init MKL dispatching
    initMKL();

    // parse command line arguments
    parseCliArgs(argc, argv);

    // print some test information
    std::cout << "nRepeat = " << g_nRepeat << "\n";
    std::cout << "nRandom = " << g_nRandom << "\n";
    std::cout << (g_testMkl ? "including" : "skipping") << " MKL tests\n";
    std::cout << (g_testOriginal ? "including" : "skipping") << " original implementation tests\n";

    // run all tests
    try {

        for (size_t i = 0; i < g_nRepeat; ++i) {
            std::cout << "Iteration: " << std::setw(2) << i + 1 << ": "
                << "generating " << g_nRandom << " 32-bits random numbers\n";
            {
                size_t m = 0;
                std::cout
                    << std::setw(s_messageSpacing[m++]) << "Generator"
                    << std::setw(s_messageSpacing[m++]) << "WordSize"
                    << std::setw(s_messageSpacing[m++]) << "RegSize"
                    << std::setw(s_messageSpacing[m++]) << "BlkSize"
                    << std::setw(s_messageSpacing[m++]) << "QueryMode"
                    << "\n";
            }

#if TEST_ORIG==1
            if (g_testOriginal) {
                // original Matsumoto - Tsukamoto MT19937 implementation
                mtOrigPerformance();
                sfmtOrigPerformance<true>(1);
                for (auto sz : anySize)
                    if (sz >= 624)
                        sfmtOrigPerformance<false>(sz);
            }
#endif

#if TEST_MKL==1
            if (g_testMkl) {
                for (auto sz : anySize)
                    mklPerformance(VSL_BRNG_MT19937, (MKL_INT)sz);
                for (auto sz : anySize)
                    mklPerformance(VSL_BRNG_SFMT19937, (MKL_INT)sz);
            }
#endif
#if TEST_VMT==1
            if (g_testVMT) {
                vRandGenPerformance0<VMT19937, /*32, */128/*, 256, 512*/>();
            }
            //vRandGenPerformance0<VSFMT19937, 128, 256, 512>();
#endif
        }

        const size_t spacing[] = { 20, 8, 8, 8, 10, 6, 8, 8, 8, 8, 11, 12 };
        size_t s = 0;
        std::cout << "\n"
            << std::setw(spacing[s++]) << std::right << "prng"
            << std::setw(spacing[s++]) << std::right << "g-bits"
            << std::setw(spacing[s++]) << std::right << "r-bits"
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
        for (auto& [k, v] : results) {
            s = 0;
            std::cout << std::setw(spacing[s++]) << std::right << modename[k.mode]
                << std::setw(spacing[s++]) << std::right << k.nBits
                << std::setw(spacing[s++]) << std::right << k.nImplBits
                << std::setw(spacing[s++]) << std::right << k.blkSize
                << std::setw(spacing[s++]) << std::right << queryModeName(k.qryMode)
                << std::setw(spacing[s++]) << std::right << v.singleRuns.size()
                << std::setw(spacing[s++]) << std::right << std::fixed << std::setprecision(3) << v.mi
                << std::setw(spacing[s++]) << std::right << std::fixed << std::setprecision(3) << v.ma
                << std::setw(spacing[s++]) << std::right << std::fixed << std::setprecision(3) << v.avg
                << std::setw(spacing[s++]) << std::right << std::fixed << std::setprecision(3) << v.stdev
                << std::setw(spacing[s++]) << std::right << std::fixed << std::setprecision(3) << v.stdev / v.avg * 100 << "%"
                << std::setw(spacing[s++]) << std::right << std::fixed << std::setprecision(1) << g_nRandom / 1.0e6 / v.avg
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
