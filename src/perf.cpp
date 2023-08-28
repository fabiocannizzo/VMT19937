#include "VRandGen.h"

#include "../SFMT-src-1.5.1/SFMT.h"

#ifdef TEST_MKL
#include "mkl.h"
#endif

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

const uint32_t seedlength = 4;
const uint32_t seedinit[seedlength] = { 0x123, 0x234, 0x345, 0x456 };

size_t nRandomPerf = size_t(624) * 16 * 800;
const size_t anySizeQuery = 1000;

extern "C" unsigned long genrand_int32();
extern "C" void init_by_array(unsigned long init_key[], int key_length);

enum Mode {orig, sfmt, mkl_mt, mkl_sfmt, vmt, vsfmt};

const char* modename[] = {"ORIG-MT19937", "ORIG-SFMT19937", "MKL-MT19937", "MKL-SFMT19937", "V-MT19937", "V-SFMT19937" };

const char* queryModeName(VRandGenQueryMode qm)
{
    switch (qm) {
        case QM_Any: return "AnySize";
        case QM_Scalar: return "Scalar";
        case QM_Block16: return "Block16";
        case QM_StateSize: return "State";
        default: THROW("how did we get here?")
    }
}

template <template <size_t, size_t> class Gen>
struct GenTraits;

// for maximum period, we should select the file based on the number of states
// but these periods are so large anyway that who do not care!
std::unique_ptr<Details::VMT19937Base<32, 32>::matrix_t> pmt(new Details::VMT19937Base<32, 32>::matrix_t("dat/mt/F19933.bits"));
// for maximum period, we should select the file based on the number of states
// but these periods are so large anyway that who do not care!
std::unique_ptr<Details::VSFMT19937Base<128, 32>::matrix_t> psfmt(new Details::VSFMT19937Base<128, 32>::matrix_t("dat/sfmt/F19935.bits"));


template <>
struct GenTraits<Details::VMT19937Base>
{
    static const Mode mode = vmt;
    static const char* name() { return "VMT19937"; }
    static const auto* matrix() { return pmt.get(); }
};

template <>
struct GenTraits<Details::VSFMT19937Base>
{
    static const Mode mode = vsfmt;
    static const char* name() { return "VSFMT19937"; }
    static const auto* matrix() { return psfmt.get(); }
};

const size_t messageSpacing[] = { 15, 9, 8, 8, 12 };

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
            << std::setw(messageSpacing[i++]) << modename[mode]
            << std::setw(messageSpacing[i++]) << nBits
            << std::setw(messageSpacing[i++]) << nImplBits
            << std::setw(messageSpacing[i++]) << blkSize
            << std::setw(messageSpacing[i++]) << queryModeName(qryMode)
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

template <template <size_t, size_t> class Gen, size_t L, size_t I, VRandGenQueryMode QryMode>
void vRandGenPerformance3()
{
    if constexpr (I <= std::min<size_t>(L, SIMD_N_BITS)) {

        typedef Details::VRandGen<Gen<L, I>, QryMode> gen_t;
        size_t blkSize;
        switch (QryMode) {
            case QM_Scalar: blkSize = 1; break;
            case QM_Block16: blkSize = 16; break;
            case QM_Any: blkSize = anySizeQuery; break;
            case QM_StateSize: blkSize = gen_t::s_n32InFullState; break;
            default: THROW("how did we get here?");
        }
        const Mode mode = GenTraits<Gen>::mode;

        ResultKey key(mode, L, I, blkSize, QryMode);

        key.print();

        if (alreadyHaveEnoughIter(key)) {
            std::cout << "skip\n";
            return;
        }

        AlignedVector<uint32_t, 64> aligneddst(blkSize);

        // we provide a jump matrix, although it is redundant for the purpose of just measuring performace
        gen_t mt(seedinit, seedlength, 0, nullptr, GenTraits<Gen>::matrix());

        auto start = std::chrono::system_clock::now();

        for (size_t i = 0, n = nRandomPerf / blkSize; i < n; ++i)
            if constexpr (QryMode == QM_Scalar)
                aligneddst[0] = mt.genrand_uint32();
            else if constexpr (QryMode == QM_Block16)
                mt.genrand_uint32_blk16(aligneddst.data());
            else if constexpr (QryMode == QM_StateSize)
                mt.genrand_uint32_stateBlk(aligneddst.data());
            else
                NOT_IMPLEMENTED;

        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        double nSeconds = elapsed_seconds.count();

        std::cout << "done in: " << std::fixed << std::setprecision(2) << nSeconds << "s" << std::endl;

        addResult(key, nSeconds);
    }
}


void originalPerformance()
{
    ResultKey key(orig, 32, 32, 1, QM_Scalar);

    key.print();

    if (alreadyHaveEnoughIter(key)) {
        std::cout << "skip\n";
        return;
    }

    std::vector<uint32_t> dst(1);

    unsigned long init[seedlength];
    for (size_t i = 0; i < seedlength; ++i)
        init[i] = seedinit[i];

    init_by_array(init, seedlength);

    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < nRandomPerf; ++i)
        dst[0] = genrand_int32();
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();
    std::cout << "done in: " << std::fixed << std::setprecision(2) << nSeconds << "s\n";

    addResult(key, nSeconds);
}

template <bool ScalarQry>
void sfmtPerformance(size_t BlkSize)
{
    ResultKey key(sfmt, 128, 128, BlkSize, QM_Any);
    key.print();
    if (alreadyHaveEnoughIter(key)) {
        std::cout << "skip\n";
        return;
    }

    MYASSERT((BlkSize == 1) || (BlkSize % 4 == 0 && BlkSize >= SFMT_N32), "BlkSize must be a multiple of 4 and >=156*128");
    MYASSERT((nRandomPerf % BlkSize) == 0, "nRandomPerf must be a multiple of BlkSize");
    AlignedVector<uint32_t, 64> aligneddst(BlkSize);

    sfmt_t sfmtgen;
    sfmt_init_gen_rand(&sfmtgen,12345);

    
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0, n = nRandomPerf / BlkSize; i < n; ++i) {
        if constexpr (ScalarQry)
            aligneddst[0] = sfmt_genrand_uint32(&sfmtgen);
        else
            sfmt_fill_array32(&sfmtgen, aligneddst.data(), (int) BlkSize);
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();
    std::cout << "done in: " << std::fixed << std::setprecision(2) << nSeconds << "s\n";

    addResult(key, nSeconds);
}

#ifdef TEST_MKL

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

    MYASSERT(nRandomPerf % BlkSize == 0, "incorrect count");

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
    for (size_t i = 0, n = nRandomPerf / BlkSize; i < n; ++i)
        viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, BlkSize, aligneddst.data());
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();
    std::cout << "done in: " << std::fixed << std::setprecision(2) << nSeconds << "s\n";

    // Deleting the stream
    vslDeleteStream(&stream);

    addResult(key, nSeconds);
}
#endif

template <template <size_t, size_t> class Gen, size_t L, size_t I, VRandGenQueryMode...QMs>
void vRandGenPerformance2()
{
    (vRandGenPerformance3<Gen, L, I, QMs>(), ...);
}

template <template <size_t, size_t> class Gen, size_t L, size_t...Is>
void vRandGenPerformance1()
{
    (vRandGenPerformance2<Gen, L, Is, QM_Scalar, QM_Block16, QM_StateSize>(), ...);
}

template <template <size_t, size_t> class Gen, size_t...Ls>
void vRandGenPerformance0()
{
    (vRandGenPerformance1<Gen, Ls, 32, 128, 256, 512>(), ...);
}

void usage()
{
    std::cerr
        << "Invalid command line arguments\n"
        << "Example:\n"
        << "perf [-n nRepeats] [-s slow]\n"
        << "  nRepeats defaults to 1\n"
        << "  slow must be 0 or 1, defaults to 0\n"
        ;
}

int main(int argc, const char** argv)
{
#ifdef TEST_MKL
# if (SIMD_N_BITS==128)
    std::cout << "Force MKL dispatching to SSE2\n";
    MYASSERT(mkl_enable_instructions(MKL_ENABLE_SSE4_2), "SSE2 not supported");
# elif (SIMD_N_BITS==256)
    std::cout << "set MKL dispacthing to AVX2\n";
    MYASSERT(mkl_enable_instructions(MKL_ENABLE_AVX2), "AVX2 not supported");
# elif (SIMD_N_BITS==512)
    std::cout << "set MKL dispacthing to AVX512\n";
    MYASSERT(mkl_enable_instructions(MKL_ENABLE_AVX512), "AVX512 not supported");
# else
   NOT_IMPLEMENTED;
# endif
#endif
    size_t nRepeat = 1;
    int slow = 1;
    // parse command line arguments
    for (int i = 1; i < argc; i += 2) {
        string key(argv[i]);
        const char* value = argv[i + 1];
        if (key == "-n")
            nRepeat = atoi(value);
        else if (key == "-s")
            slow = atoi(value);
        else {
            usage();
            exit(1);
        }
    }
    std::cout << "nRepeat = " << nRepeat << "\n";
    std::cout << "slow = " << slow << "\n";
    if (slow)
        nRandomPerf *= 1000;

    for (size_t i = 0; i < nRepeat; ++i) {
        std::cout << "Iteration: " << std::setw(2) << i + 1 << ": " 
                  << "generating " << nRandomPerf << " 32-bits random numbers\n";
        {
            size_t m = 0;
            std::cout
                << std::setw(messageSpacing[m++]) << "Generator"
                << std::setw(messageSpacing[m++]) << "WordSize"
                << std::setw(messageSpacing[m++]) << "RegSize"
                << std::setw(messageSpacing[m++]) << "BlkSize"
                << std::setw(messageSpacing[m++]) << "QueryMode"
                << "\n";
        }

        originalPerformance();
        sfmtPerformance<true>(1);
        sfmtPerformance<false>(156 * 128 / 32);
#ifdef TEST_MKL
        MKL_INT mklQrySizeMT[] = { 1, 16, 624, 624 * 4, 624 * 8, 624 * 16 };
        for (auto sz : mklQrySizeMT)
            mklPerformance(VSL_BRNG_MT19937, sz);
        MKL_INT mklQrySizeSFMT[] = { 1, 16, 624, 624 * 2, 624 * 4 };
        for (auto sz : mklQrySizeSFMT)
            mklPerformance(VSL_BRNG_SFMT19937, sz);
#endif
        vRandGenPerformance0<Details::VMT19937Base, 32, 128, 256, 512>();
        vRandGenPerformance0<Details::VSFMT19937Base, 128, 256, 512>();
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
        << std::setw(1+spacing[s++]) << std::right << "tavg/tdev"
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
            << std::setw(spacing[s++]) << std::right << std::fixed << std::setprecision(3) << v.stdev /v.avg*100 << "%"
            << std::setw(spacing[s++]) << std::right << std::fixed << std::setprecision(1) << nRandomPerf / 1.0e6 / v.avg
            << "\n";
    }

    return 0;
}
