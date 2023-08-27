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

template <template <size_t, size_t> class Gen>
struct GenTraits;

template <>
struct GenTraits<Details::VMT19937Base>
{
    static const Mode mode = vmt;
    static const char* name() { return "VMT19937"; }
    // for maximum period, we should select the file based on the number of states
    // but these periods are so large anyway that who do not care!
    static const char* jumpFileName() { return "dat/mt/F19933.bits"; }
    typedef Details::VMT19937Base<32,32>::matrix_t matrix_t;
};

template <>
struct GenTraits<Details::VSFMT19937Base>
{
    static const Mode mode = vsfmt;
    static const char* name() { return "VSFMT19937"; }
    // for maximum period, we should select the file based on the number of states
    // but these periods are so large that anyway we do not care!
    static const char* jumpFileName() { return "dat/sfmt/F19935.bits"; }
    typedef Details::VSFMT19937Base<128, 32>::matrix_t matrix_t;
};

struct Result
{
    Result(Mode _mode, size_t _nb, size_t _ib, size_t _blk, size_t _qryMode, size_t _id, double _dt)
        : mode(_mode), nBits(_nb), nImplBits(_ib), blkSize(_blk), qryMode(_qryMode), runId(_id), time(_dt) {}
    Mode mode;
    size_t nBits, nImplBits;
    size_t blkSize;
    size_t qryMode;
    mutable size_t runId;
    mutable double time;
    bool operator<(const Result& rhs) const {
        return std::tuple(mode, nBits, nImplBits, blkSize, qryMode) < std::tuple(rhs.mode, rhs.nBits, rhs.nImplBits, rhs.blkSize, rhs.qryMode);
    }
};

template <template <size_t, size_t> class Gen, size_t L, size_t I, VRandGenQueryMode QM, typename M>
void vRandGenPerformance3(size_t runId, const M* m, std::vector<Result>& res)
{
    if constexpr (I <= std::min<size_t>(L, SIMD_N_BITS)) {

        typedef Details::VRandGen<Gen<L, I>, QM> gen_t;
        const size_t VecLen = L;
        const VRandGenQueryMode QryMode = QM;
        const size_t BlkSize = QryMode == QM_Scalar ? 1 : QryMode == QM_Block16 ? 16 : gen_t::s_n32InFullState;
        const Mode mode = GenTraits<Gen>::mode;
        const char* name = GenTraits<Gen>::name();

        std::cout << "Generate " << nRandomPerf << " random numbers with " << name << "<" << std::setw(3) << VecLen << "," << std::setw(3) << I << ">"
            << " QrySize=" << BlkSize << " ... ";

        AlignedVector<uint32_t, 64> aligneddst(BlkSize);

        // we provide a jump matrix, although it is redundant for the purpose of just measuring performace
        gen_t mt(seedinit, seedlength, 0, nullptr, m);

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

        res.emplace_back(mode, VecLen, I, BlkSize, BlkSize, runId, nSeconds);
    }
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

    return Result(orig, 32, 32, 1, runId, 1, nSeconds);
}

template <bool ScalarQry>
Result sfmtPerformance(size_t runId, size_t BlkSize)
{
    MYASSERT((BlkSize == 1) || (BlkSize % 4 == 0 && BlkSize >= SFMT_N32), "BlkSize must be a multiple of 4 and >=156*128");
    MYASSERT((nRandomPerf % BlkSize) == 0, "nRandomPerf must be a multiple of BlkSize");
    AlignedVector<uint32_t, 64> aligneddst(BlkSize);

    sfmt_t sfmtgen;
    sfmt_init_gen_rand(&sfmtgen, 12345);

    std::cout << "Generate " << nRandomPerf << " random numbers with SFMT code QrySize=" << BlkSize << " ... ";
    auto start = std::chrono::system_clock::now();
    for (size_t i = 0, n = nRandomPerf / BlkSize; i < n; ++i) {
        if constexpr (ScalarQry)
            aligneddst[0] = sfmt_genrand_uint32(&sfmtgen);
        else
            sfmt_fill_array32(&sfmtgen, aligneddst.data(), BlkSize);
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();
    std::cout << "done in: " << std::fixed << std::setprecision(2) << nSeconds << "s\n";

    return Result(sfmt, 128, 128, BlkSize, BlkSize == 1 ? 1 : 0, runId, nSeconds);
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

template <MKL_INT GenCode>
Result mklPerformance(size_t runId, MKL_INT BlkSize)
{
    MYASSERT(nRandomPerf % BlkSize == 0, "incorrect count");

    AlignedVector<uint32_t, 64> aligneddst(BlkSize);

    VSLStreamStatePtr stream;
    vslNewStream(&stream, GenCode, 5489);

    std::cout << "Generate " << nRandomPerf << " random numbers with MKL "
              << ((GenCode == VSL_BRNG_MT19937) ? "MT19937" : "SFMT19937")
              << " QrySize=" << BlkSize << " ... ";

    auto start = std::chrono::system_clock::now();
    for (size_t i = 0, n = nRandomPerf / BlkSize; i < n; ++i)
        viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, BlkSize, aligneddst.data());
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    double nSeconds = elapsed_seconds.count();
    std::cout << "done in: " << std::fixed << std::setprecision(2) << nSeconds << "s\n";

    // Deleting the stream
    vslDeleteStream(&stream);

    typedef MKLTraits<GenCode> traits;
    return Result(traits::s_mode, traits::s_wordSize, traits::s_wordSize, BlkSize, 0, runId, nSeconds);
}
#endif

template <template <size_t, size_t> class Gen, size_t L, size_t I, VRandGenQueryMode...QMs, typename M>
void vRandGenPerformance2(size_t runId, const M* m, std::vector<Result>& res)
{
    (vRandGenPerformance3<Gen, L, I, QMs>(runId, m, res), ...);
}

template <template <size_t, size_t> class Gen, size_t L, size_t...Is, typename M>
void vRandGenPerformance1(size_t runId, const M* m, std::vector<Result>& res)
{
    (vRandGenPerformance2<Gen, L, Is, QM_Scalar, QM_Block16, QM_StateSize>(runId, m, res), ...);
}

template <template <size_t, size_t> class Gen, size_t...Ls>
void vRandGenPerformance0(size_t runId, std::vector<Result>& res)
{
    typedef typename GenTraits<Gen>::matrix_t matrix_t;
    std::unique_ptr<matrix_t> p(new matrix_t(GenTraits<Gen>::jumpFileName()));
    (vRandGenPerformance1<Gen, Ls, 32, 128, 256, 512>(runId, p.get(), res), ...);
}

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
        res.push_back(sfmtPerformance<true>(i, 1));
        res.push_back(sfmtPerformance<false>(i, 156ul * 128 / 32));
#ifdef TEST_MKL
        MKL_INT mklQrySizeMT[] = { 1, 16, 624, 624 * 4, 624 * 8, 624 * 16 };
        MKL_INT mklQrySizeSFMT[] = { 1, 16, 624, 624 * 2, 624 * 4 };
        for (auto sz : mklQrySizeMT)
            res.push_back(mklPerformance<VSL_BRNG_MT19937>(i, sz));
        for (auto sz : mklQrySizeSFMT)
            res.push_back(mklPerformance<VSL_BRNG_SFMT19937>(i, sz));
#endif
        vRandGenPerformance0<Details::VMT19937Base, 32, 128, 256, 512>(i, res);
        vRandGenPerformance0<Details::VSFMT19937Base, 128, 256, 512>(i, res);
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
        << std::setw(8) << std::right << "i-bits"
        << std::setw(8) << std::right << "blksize"
        << std::setw(8) << std::right << "qrymode"
        << std::setw(8) << std::right << "time"
        << "\n";
    for (auto& r : avgresult) {
        std::cout << std::setw(12) << std::right << modename[r.mode]
            << std::setw(8) << std::right << r.nBits
            << std::setw(8) << std::right << r.nImplBits
            << std::setw(8) << std::right << r.blkSize
            << std::setw(8) << std::right << r.qryMode
            << std::setw(8) << std::right << std::fixed << std::setprecision(2) << r.time / nRepeat
            << "\n";
    }

    return 0;
}
