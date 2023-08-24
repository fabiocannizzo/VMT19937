#include "TestUtils.h"

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

using namespace std;

#define TEST_MKL 1
#define TEST_VMT 1
#define TEST_VSFMT 1
#define TEST_XMT 1
#define TEST_ORIG 1

#if TEST_MKL==1
#   if __has_include(<mkl.h>)
#       include <mkl.h>
#   else
#       pragma message("MKL not found, disabling MKL tests")
#       undef TEST_MKL
#       define TEST_MKL 0
#   endif
#endif

#if TEST_ORIG==1
#   define HAVE_SSE2
#   define SFMT_MEXP 19937
#   include "../SFMT-src-1.5.1/SFMT.h"
#endif

// global variables
bool g_testMkl = TEST_MKL;
bool g_testOriginal = true;
bool g_testVMT = true;
bool g_testXMT = true;
bool g_testSFMT = true;
bool g_testQry1 = true;
bool g_testQryN = true;
bool g_testQry16 = true;
size_t g_nRepeat = 1;
double g_stDev = 0.0;
std::string dir = "dat";

// this might be changed via cli arguments
size_t g_nRandom = size_t(624) * 32 * 800;

constexpr uint32_t s_seedlength = 4;
constexpr uint32_t s_seedinit[s_seedlength] = { 0x123, 0x234, 0x345, 0x456 };


extern "C" unsigned long genrand_int32();
extern "C" void init_by_array(unsigned long init_key[], int key_length);

enum GenMode {orig, sfmt, mkl_mt, mkl_sfmt, xmt, vmt, vsfmt};

const char* modename[] = {"ORIG-MT19937", "ORIG-SFMT19937", "MKL-MT19937", "MKL-SFMT19937", "X-MT19937", "V-MT19937", "V-SFMT19937" };

const size_t anySize[] = {/* 1, 4, 16, 64, 256, 624, 1024, 4096,*/ 16384 };

template <GenMode G>
struct GenTraits;

// for maximum period, we should select the file based on the number of states
// but these periods are so large anyway that who do not care!
const auto pmt = std::make_unique<MT19937Matrix>(dir + "/mt/F19933.bits");
// for maximum period, we should select the file based on the number of states
// but these periods are so large anyway that who do not care!
const auto psfmt = std::make_unique<SFMT19937Matrix>(dir + "/sfmt/F19935.bits");

// use the same destination memory in all tests to avoid spurious difference in test results due to memory layout
AlignedVector<uint32_t, 64> aligneddst(anySize[sizeof(anySize)/sizeof(anySize[0])-1]);

template <>
struct GenTraits<vmt>
{
    static const GenMode mode = vmt;
    static const MT19937Matrix* jumpMatrix() { return pmt.get(); }

    template <size_t RegBitLen, QryMode QM, size_t RegBitLenHw>
    using gen_t = VMT19937<RegBitLen, QM == QM_Block16, RegBitLenHw>;
};

template <>
struct GenTraits<xmt>
{
    static const GenMode mode = xmt;
    static const MT19937Matrix* jumpMatrix() { return nullptr; }

    template <size_t RegBitLen, QryMode QM, size_t RegBitLenHw, std::enable_if_t<RegBitLen == RegBitLenHw, int> = 0>
    using gen_t = XMT19937<RegBitLen, QM == QM_Block16>;
};

template <>
struct GenTraits<vsfmt>
{
    static const GenMode mode = vsfmt;
    static const SFMT19937Matrix* jumpMatrix() { return psfmt.get(); }

    template <size_t RegBitLen, QryMode QM, size_t RegBitLenHw>
    using gen_t = VSFMT19937<RegBitLen, QM == QM_Block16, RegBitLenHw>;
};

const size_t s_messageSpacing[] = { 15, 9, 8, 8, 12 };

struct Results
{
    Results(GenMode _mode, size_t _nb, size_t _ib, size_t _blk, QryMode _qryMode)
        : mode(_mode), nBits(_nb), nBitsHw(_ib), blkSize(_blk), qryMode(_qryMode), mi(0), ma(0), avg(0), stdev(0)
    {}
    GenMode mode;
    size_t nBits, nBitsHw;
    size_t blkSize;
    QryMode qryMode;

    mutable std::vector<double> singleRuns;
    mutable double mi, ma, avg, stdev;

    void print() const
    {
        size_t i = 0;
        std::cout
            << std::setw(s_messageSpacing[i++]) << modename[mode]
            << std::setw(s_messageSpacing[i++]) << nBits
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

__declspec(noinline)
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

#if TEST_ORIG==1
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

template <GenMode Mode, size_t L, size_t I, QryMode QM>
void vRandGenPerformance5(size_t blkSize)
{
    MYASSERT(((blkSize > 0) && ((g_nRandom % blkSize) == 0)), "invalid blkSize " << blkSize);

    Results key(Mode, L, I, blkSize, QM);

    key.print();

    if (alreadyHaveEnoughIter(key)) {
        std::cout << "skip\n";
        return;
    }

    using Gen = typename GenTraits<Mode>::template gen_t<L, QM, I>;

    const typename Gen::matrix_t *jumpMatrixPtr = nullptr;
    if constexpr (Gen::s_nStates > 1)
        jumpMatrixPtr = GenTraits<Mode>::jumpMatrix();
    Gen mt(s_seedinit, s_seedlength, 0, nullptr, jumpMatrixPtr);

    auto start = std::chrono::system_clock::now();

    for (size_t i = 0, n = g_nRandom / blkSize; i < n; ++i) {
        if constexpr (QM == QM_Scalar)
            aligneddst[0] = mt.genrand_uint32();
        else if constexpr (QM == QM_Block16)
            mt.genrand_uint32_blk16(aligneddst.data());
        else if constexpr (QM == QM_Any)
            mt.genrand_uint32_anySize(aligneddst.data(), blkSize);
        else
            NOT_IMPLEMENTED;
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
    if constexpr (QM == QM_Any) {
        if (g_testQryN)
            for (auto sz : anySize)
                vRandGenPerformance5<Mode, L, I, QM>(sz);
    }
    else if constexpr (QM == QM_Block16) {
        if (g_testQry16)
            vRandGenPerformance5<Mode, L, I, QM>(16);
    }
    else if constexpr (QM == QM_Scalar) {
        if (g_testQry1)
            vRandGenPerformance5<Mode, L, I, QM>(1);
    }
}

template <GenMode Mode, size_t L, size_t I, QryMode...QMs>
void vRandGenPerformance2()
{
    constexpr size_t M = std::min<size_t>(L, SIMD_N_BITS);
    if constexpr (I <= M && (Mode != xmt || I == L))
        (vRandGenPerformance4<Mode, L, I, QMs>(), ...);
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
    std::cerr
        << "Invalid command line arguments\n"
        << "perf [-n nRepeats] [-s nRndScaler] [--no-mkl] [--no-original] [--no-vmt] [--no-sfmt]\n"
        << "     [--no - xmt][--dir datpath][--no - qry1][--no - qry16][--no - qryN][--stdev 1.0]\n"
        << "  -n nRepeats: number of performance test iterations (default 1)\n"
        << "  --no-mkl: skip MKL tests\n"
        << "  --no-vmt: skip VMT tests\n"
        << "  --no-xmt: skip XMT tests\n"
        << "  --no-sfmt: skip SFMT tests\n"
        << "  --no-qry1: skip scalar query tests\n"
        << "  --no-qry16: skip block-16 query tests\n"
        << "  --no-qryN: skip vectorial query tests\n"
        << "  --no-original: skip original implementations tests\n"
        << "  --slow: increase the number of random numbers generated by 1000 times\n"
        << "  --dir datpath: folder where to find jump matrix files"
        << "  --stdev 1.0: stop iteratng when stdev/avg<1.0% (nRepeats is treated like a minimum)"
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
            else if (key == "--no-sfmt")
                g_testSFMT = false;
            else if (key == "--no-xmt")
                g_testVMT = false;
            else if (key == "--no-qry1")
                g_testQry1 = false;
            else if (key == "--no-qry16")
                g_testQry16 = false;
            else if (key == "--no-qryN")
                g_testQryN = false;
            else if (key == "-s") {
                MYASSERT(++i < argc, "-s must be followed by a number");
                const char* value = argv[i];
                g_nRandom *= stoul(string(value));
            }
            else if (key == "--stdev") {
                MYASSERT(++i < argc, "--stdev must be a number");
                const char* value = argv[i];
                g_stDev = stod(string(value));
            }
            else if (key == "--dir") {
            	MYASSERT(++i < argc, "--dir must be followed by a path");
                dir = argv[i];
            }
            else {
                THROW("Invalid command line argument: " << argv[i]);
            }
        }
        if (g_stDev > 0.0)
            g_nRepeat = std::max<size_t>(g_nRepeat, 2);
        MYASSERT(g_nRepeat > 0, "nRepeat must be positive");
    }
    catch (const std::exception& ex) {
        std::cerr << "Error parsing command line arguments: " << ex.what() << "\n";
        syntax();
        std::exit(-1);
    }
    catch (...) {
        syntax();
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
    std::cout << "Target hardware SIMD register size (bits): " << SIMD_N_BITS << "\n";
    std::cout << "nRepeat = " << g_nRepeat << "\n";
    std::cout << "stDev = " << g_stDev << "\n";
    std::cout << "nRandom = " << g_nRandom << "\n";
    std::cout << (g_testMkl ? "including" : "skipping") << " MKL tests\n";
    std::cout << (g_testOriginal ? "including" : "skipping") << " original implementation tests\n";

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
                if (g_testSFMT) {
                    if (g_testQry1)
                        sfmtOrigPerformance<true>(1);
                    if (g_testQryN) {
                        for (auto sz : anySize)
                            if (sz >= 624)
                                sfmtOrigPerformance<false>(sz);
                    }
                }
            }
#endif

#if TEST_MKL==1
            if (g_testMkl) {
                for (auto sz : anySize)
                    mklPerformance(VSL_BRNG_MT19937, (MKL_INT)sz);
                if (g_testSFMT) {
                    for (auto sz : anySize)
                        mklPerformance(VSL_BRNG_SFMT19937, (MKL_INT)sz);
                }
            }
#endif
#if TEST_VMT==1
            if (g_testVMT) {
                vRandGenPerformance0<vmt, /*32, 128, 256, 512*/ SIMD_N_BITS>();
            }
#endif
#if TEST_VSFMT==1
            if (g_testSFMT) {
                vRandGenPerformance0<vsfmt, /*32, 128, 256, 512*/ SIMD_N_BITS>();
            }
#endif
#if TEST_XMT==1
            if (g_testXMT) {
                vRandGenPerformance0<xmt, /*32, 128, 256, 512*/ SIMD_N_BITS>();
            }
#endif
            size_t nResAfter = nResults();
            if (nResAfter == nResBefore)
                break; // no new results added
        }

        std::set<Results, TableCompare> sortedResults(results.begin(), results.end());

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
        for (auto& r : sortedResults) {
            s = 0;
            std::cout << std::setw(spacing[s++]) << std::right << modename[r.mode]
                << std::setw(spacing[s++]) << std::right << r.nBits
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
