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
#include <set>

using namespace std;

const uint32_t seedlength = 4;
const uint32_t seedinit[seedlength] = { 0x123, 0x234, 0x345, 0x456 };

const uint64_t nRandomPerf = uint64_t(624) * 16 * 800000;

extern "C" unsigned long genrand_int32();
extern "C" void init_by_array(unsigned long init_key[], int key_length);

enum Mode {orig, sfmt, mkl_mt, mkl_sfmt, vmt, vsfmt};

const char* modename[] = {"MT19937", "SFMT19937", "VSLMT19937", "VSLSFMT19937", "VMT19937", "VSFMT19937" };

template <typename T>
struct GenTraits;

template <size_t VecLen, VRandGenQueryMode QryMode>
struct GenTraits<VMT19937<VecLen, QryMode>>
{
    static const Mode mode = vmt;
    static const char* name() { return "VMT19937"; }
    // for maximum period, we should select the file based on the number of states
    // but these periods are so large anyway that who do not care!
    static const char* jumpFileName() { return "dat/mt/F19933.bits"; }
};

template <size_t VecLen, VRandGenQueryMode QryMode>
struct GenTraits<VSFMT19937<VecLen, QryMode>>
{
    static const Mode mode = vsfmt;
    static const char* name() { return "VSFMT19937"; }
    // for maximum period, we should select the file based on the number of states
    // but these periods are so large that anyway we do not care!
    static const char* jumpFileName() { return "dat/sfmt/F19935.bits"; }
};

struct Result
{
    Result(Mode _mode, size_t _nb, size_t _blk, size_t _id, double _dt)
        : mode(_mode), nBits(_nb), blkSize(_blk), runId(_id), time(_dt) {}
    Mode mode;
    size_t nBits;
    size_t blkSize;
    mutable size_t runId;
    mutable double time;
    bool operator<(const Result& rhs) const {
        return std::tuple(mode, nBits, blkSize) < std::tuple(rhs.mode, rhs.nBits, rhs.blkSize);
    }
};

template <typename Gen>
Result testPerformance(size_t runId)
{
    const size_t VecLen = Gen::s_regLenBits;
    const VRandGenQueryMode QryMode = Gen::s_queryMode;
    const size_t BlkSize = QryMode == QM_Scalar ? 1 : QryMode == QM_Block16 ? 16 : Gen::s_n32InFullState;
    const Mode mode = GenTraits<Gen>::mode;
    const char* name = GenTraits<Gen>::name();

    std::cout << "Generate " << nRandomPerf << " random numbers with " << name << " length " << VecLen
              << " in blocks of " << BlkSize << " ... ";

    AlignedVector<uint32_t, 64> aligneddst(BlkSize);

    typedef typename Gen::matrix_t matrix_t;
    std::unique_ptr<matrix_t> pjump(new matrix_t(GenTraits<Gen>::jumpFileName()));    // load jump ahead matrix
    Gen mt(seedinit, seedlength, 0, nullptr, pjump.get());
    pjump.reset(nullptr);

    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < nRandomPerf / BlkSize; ++i)
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

    return Result(mode, VecLen, BlkSize, runId, nSeconds);
}


Result originalPerformance(size_t runId)
{
    std::vector<uint32_t> dst(1);

    unsigned long init[seedlength];
    for (size_t i = 0; i < seedlength; ++i)
        init[i] = seedinit[i];

    init_by_array(init, seedlength);
    std::cout << "Generate " << nRandomPerf << " random numbers with original code ... ";
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < nRandomPerf; ++i)
        dst[0] = genrand_int32();
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();
    std::cout << "done in: " << std::fixed << std::setprecision(2) << nSeconds << "s\n";

    return Result(orig, 32, 1, runId, nSeconds);
}

template <size_t BlkSize>
Result sfmtPerformance(size_t runId)
{
    static_assert(BlkSize == 1 || (BlkSize % 4 == 0 && BlkSize >= SFMT_N32), "BlkSize must be a multiple of 4 and >=156*128");
    static_assert(nRandomPerf % BlkSize == 0, "nRandomPerf must be a multiple of BlkSize");
    AlignedVector<uint32_t, 64> aligneddst(BlkSize);

    sfmt_t sfmtgen;
    sfmt_init_gen_rand(&sfmtgen, 12345);

    std::cout << "Generate " << nRandomPerf << " random numbers with SFMT code in blocks of " << BlkSize << " ... ";
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < nRandomPerf / BlkSize; ++i) {
        if constexpr (BlkSize == 1)
            aligneddst[0] = sfmt_genrand_uint32(&sfmtgen);
        else
            sfmt_fill_array32(&sfmtgen, aligneddst.data(), BlkSize);
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();
    std::cout << "done in: " << std::fixed << std::setprecision(2) << nSeconds << "s\n";

    return Result(sfmt, 128, BlkSize, runId, nSeconds);
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

template <MKL_INT GenCode, MKL_INT BlkSize>
Result mklPerformance(size_t runId)
{
    AlignedVector<uint32_t, 64> aligneddst(BlkSize);

    VSLStreamStatePtr stream;
    vslNewStream(&stream, GenCode, 5489);

    std::cout << "Generate " << nRandomPerf << " random numbers with MKL "
              << ((GenCode == VSL_BRNG_MT19937) ? "MT19937" : "SFMT19937")
              << " in blocks of " << BlkSize << " ... ";

    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < nRandomPerf / BlkSize; ++i)
        viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, BlkSize, aligneddst.data());
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();
    std::cout << "done in: " << std::fixed << std::setprecision(2) << nSeconds << "s\n";

    // Deleting the stream
    vslDeleteStream(&stream);

    typedef MKLTraits<GenCode> traits;
    return Result(traits::s_mode, traits::s_wordSize, BlkSize, runId, nSeconds);
}
#endif

void usage()
{
    std::cerr
        << "Invalid command line arguments\n"
        << "Example:\n"
        << "perf [-n nRepeats]\n"
        << "  nRepeats defaults to 1\n"
        ;
    std::exit(-1);
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
        res.push_back(sfmtPerformance<1>(i));
        res.push_back(sfmtPerformance<156*128/32>(i));
#ifdef TEST_MKL
        res.push_back(mklPerformance<VSL_BRNG_MT19937, 1>(i));
        res.push_back(mklPerformance<VSL_BRNG_MT19937, 16>(i));
        res.push_back(mklPerformance<VSL_BRNG_MT19937, 624>(i));
        res.push_back(mklPerformance<VSL_BRNG_MT19937, 624 * 4>(i));
        res.push_back(mklPerformance<VSL_BRNG_MT19937, 624 * 8>(i));
        res.push_back(mklPerformance<VSL_BRNG_MT19937, 624 * 16>(i));
        res.push_back(mklPerformance<VSL_BRNG_SFMT19937, 1>(i));
        res.push_back(mklPerformance<VSL_BRNG_SFMT19937, 16>(i));
        res.push_back(mklPerformance<VSL_BRNG_SFMT19937, 624>(i));
        res.push_back(mklPerformance<VSL_BRNG_SFMT19937, 624 * 2>(i));
        res.push_back(mklPerformance<VSL_BRNG_SFMT19937, 624 * 4>(i));
#endif
        res.push_back(testPerformance<VMT19937<32, QM_Scalar>>(i));
        res.push_back(testPerformance<VMT19937<32, QM_Block16>>(i));
        res.push_back(testPerformance<VMT19937<32, QM_StateSize>>(i));
        res.push_back(testPerformance<VMT19937<128, QM_Scalar>>(i));
        res.push_back(testPerformance<VMT19937<128, QM_Block16>>(i));
        res.push_back(testPerformance<VMT19937<128, QM_StateSize>>(i));
        res.push_back(testPerformance<VSFMT19937<128, QM_Scalar>>(i));
        res.push_back(testPerformance<VSFMT19937<128, QM_Block16>>(i));
        res.push_back(testPerformance<VSFMT19937<128, QM_StateSize>>(i));
#if SIMD_N_BITS >= 256
        res.push_back(testPerformance<VMT19937<256, QM_Scalar>>(i));
        res.push_back(testPerformance<VMT19937<256, QM_Block16>>(i));
        res.push_back(testPerformance<VMT19937<256, QM_StateSize>>(i));
        res.push_back(testPerformance<VSFMT19937<256, QM_Scalar>>(i));
        res.push_back(testPerformance<VSFMT19937<256, QM_Block16>>(i));
        res.push_back(testPerformance<VSFMT19937<256, QM_StateSize>>(i));
#endif
#if SIMD_N_BITS >= 512
        res.push_back(testPerformance<VMT19937<512, QM_Scalar>>(i));
        res.push_back(testPerformance<VMT19937<512, QM_Block16>>(i));
        res.push_back(testPerformance<VMT19937<512, QM_StateSize>>(i));
        res.push_back(testPerformance<VSFMT19937<512, QM_Scalar>>(i));
        res.push_back(testPerformance<VSFMT19937<512, QM_Block16>>(i));
        res.push_back(testPerformance<VSFMT19937<512, QM_StateSize>>(i));
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
    std::cout << std::setw(12) << std::right << "prng"
        << std::setw(8) << std::right << "n-bits"
        << std::setw(8) << std::right << "blksize"
        << std::setw(8) << std::right << "time"
        << "\n";
    for (auto& r : avgresult) {
        std::cout << std::setw(12) << std::right << modename[r.mode]
            << std::setw(8) << std::right << r.nBits
            << std::setw(8) << std::right << r.blkSize
            << std::setw(8) << std::right << std::fixed << std::setprecision(2) << r.time / nRepeat
            << "\n";
    }

    return 0;
}
