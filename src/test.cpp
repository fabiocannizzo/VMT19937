#define RANDGEN_TESTING 1

#include "TestUtils.h"
#include "SIMD.h"
#include "cli_args.h"

#define SFMT_MEXP 19937
#include "../SFMT-src-1.5.1/SFMT.h"

#include <utility>
#include <numeric>
#include <random>
#include <fstream>
#include <iomanip>

using namespace xvmt;
using namespace xvmt::details;

const uint32_t g_seedlength = 4;
const uint32_t g_seedinit[g_seedlength] = { 0x123, 0x234, 0x345, 0x456 };

const uint64_t g_nRandomTest = 100;
const uint64_t g_nRandomTest64 = 50;

extern "C" void init_genrand(unsigned long s);
extern "C" unsigned long genrand_int32();
extern "C" void init_by_array(unsigned long init_key[], int key_length);

extern "C" void init_genrand64(unsigned long long seed);
extern "C" unsigned long long genrand64_int64();

enum GenType { VSFMT, VMT, XMT, VMT64, XMT64 };
const char* g_genName[] = { "VSFMT", "VMT", "XMT", "VMT64", "XMT64" };

template <GenType G, size_t L, size_t I, QryMode QM>
struct GenTraits;

template <size_t L, size_t I, QryMode QM>
struct GenTraits<VMT, L, I, QM>
{
    typedef VMT19937<L, QM == QM_Block16, BitLenToIsa<I>::s_isa> gen_t;
};

template <size_t L, size_t I, QryMode QM>
struct GenTraits<XMT, L, I, QM>
{
    static_assert(L == I);
    typedef XMT19937<BitLenToIsa<I>::s_isa, QM == QM_Block16> gen_t;
};

template <size_t L, size_t I, QryMode QM>
struct GenTraits<VSFMT, L, I, QM>
{
    typedef VSFMT19937<L, QM == QM_Block16, BitLenToIsa<I>::s_isa> gen_t;
};

template <size_t L, size_t I, QryMode QM>
struct GenTraits<VMT64, L, I, QM>
{
    typedef VMT19937_64<L, QM == QM_Block16, BitLenToIsa<I>::s_isa> gen_t;
};

template <size_t L, size_t I, QryMode QM>
struct GenTraits<XMT64, L, I, QM>
{
    static_assert(L == I);
    typedef XMT19937_64<BitLenToIsa<I>::s_isa, QM == QM_Block16> gen_t;
};

std::vector<uint32_t> g_benchmark(g_nRandomTest + 10000);
std::vector<uint64_t> g_benchmark64(g_nRandomTest64 + 10000);

void printSome(const std::vector<uint32_t>& v)
{
    std::cout << "\n";
    for (size_t i = 0; i < 16; ++i)
        std::cout << std::setw(10) << v[i] << ((i + 1) % 8 == 0 ? "\n" : " ");
    std::cout << "...\n";
    for (size_t i = 240; i < 240 + 16; ++i)
        std::cout << std::setw(10) << v[i] << ((i + 1) % 8 == 0 ? "\n" : " ");
    std::cout << "...\n";
    for (size_t i = 624 - 16; i < 624; ++i)
        std::cout << std::setw(10) << v[i] << ((i + 1) % 8 == 0 ? "\n" : " ");
    std::cout << "\n";
}

enum EncodeMode {Base64, Hex};

template <size_t nRows, size_t nCols>
void testEncoder(const details::BinaryMatrix<nRows, nCols>& m, EncodeMode enc)
{
    const char* modename = enc == Base64 ? "base64" : "hex";
    details::BinaryMatrix<nRows, nCols> m2;
    std::cout << "saving matrix to " << modename << " stream" << std::endl;
    std::ostringstream os;
    if (enc == Base64) m.toBase64(os); else m.toHex(os);
    std::cout << "first 32 characters of the stream" << std::endl;
    std::string s = os.str();
    for (size_t i = 0; i < 32; ++i) std::cout << s[i];
    std::cout << std::endl;
    std::cout << "reading back the matrix from " << modename << " stream" << std::endl;
    std::istringstream is(os.str());
    if (enc == Base64) m2.fromBase64(is); else m2.fromHex(is);
    std::cout << "compare with original matrix" << std::endl;
    MYASSERT((m == m2), "error in roundtrip");
    std::cout << "completed" << std::endl;
}

template <size_t NBITS>
void testSquare(const details::BinarySquareMatrix<NBITS>& m)
{
    details::BinarySquareMatrix<NBITS> m2, m3;
    for (size_t r = 0; r < NBITS; ++r) {
        for (size_t c = 0; c < NBITS; ++c) {
            size_t s = 0;
            for (size_t k = 0; k < NBITS; ++k) s ^= m.getBit(r, k) && m.getBit(k, c);
            if (s) m2.setBit(r, c);
        }
    }
    const size_t nThreads = 4;
    std::vector<typename details::BinarySquareMatrix<NBITS>::buffer_t> buffers(nThreads);
    m3.square(m, buffers);
    MYASSERT((m2 == m3), "error in square");
}

template <size_t nRows, size_t nCols>
void encodingTests()
{
    details::BinaryMatrix<nRows, nCols> m;
    m.initRand();
    std::cout << "\ngenerated random matrix with size (" << m.s_nBitRows << "x" << m.s_nBitCols << ") with " << m.nnz() << " non zero elements" << std::endl;
    m.printBits(0, 0, 10, 32);
    testEncoder(m, Base64);
    testEncoder(m, Hex);
}

template <size_t NBits>
void squareTest()
{
    std::cout << "testing multiplication with matrices of size: " << NBits << "\n";
    details::BinarySquareMatrix<NBits> m;
    for (size_t i = 0; i < 10; ++i) {
        m.resetZero();
        m.initRand();
        testSquare(m);
    }
}

template <size_t...NBits>
void squareTests(std::index_sequence<NBits...>&&) { (squareTest<NBits>(), ...); }

void generateBenchmark_MT19937()
{
    unsigned long init[g_seedlength];
    for (size_t i = 0; i < g_seedlength; ++i) init[i] = g_seedinit[i];
    std::cout << "Generate MT19937 random numbers with the original C source code ... ";
    init_by_array(init, g_seedlength);
    for (size_t i = 0, n  = g_benchmark.size(); i < n; ++i) g_benchmark[i] = (uint32_t)genrand_int32();
    std::cout << "done!\n";
    printSome(g_benchmark);
}

void generateBenchmark_SFMT19937()
{
    std::cout << "Generate SFMT19937 random numbers with the original C source code ... ";
    sfmt_t sfmtgen;
    sfmt_init_by_array(&sfmtgen, const_cast<uint32_t *>(g_seedinit), g_seedlength);
    for (size_t i = 0, n = g_benchmark.size(); i < n; ++i) g_benchmark[i] = sfmt_genrand_uint32(&sfmtgen);
    std::cout << "done!\n";
    printSome(g_benchmark);
}

void generateBenchmark_MT19937_64()
{
    std::cout << "Generate MT19937-64 random numbers with the original C source code ... ";
    init_genrand64(5489ULL);
    for (size_t i = 0, n = g_benchmark64.size(); i < n; ++i) g_benchmark64[i] = genrand64_int64();
    std::cout << "done!\n";
}

void startTest(const char* name)
{
    std::cout << "\n" << std::setw(40) << std::setfill('*') << "" << "\nTest " << name << "\n" << std::setw(40) << std::setfill('*') << "" << std::endl;
    std::cout << std::setfill(' ');
}

void test_STL_MT19937()
{
    startTest("STL mt19937 equivalence");
    constexpr unsigned long seed = 5489UL;
    constexpr size_t N = 10000;
    init_genrand(seed);
    std::vector<uint32_t> ref(N);
    for (size_t i = 0; i < N; ++i) ref[i] = (uint32_t)genrand_int32();
    std::mt19937 stl(seed);
    for (size_t i = 0; i < N; ++i) { uint32_t v = stl(); MYASSERT(v == ref[i], "FAILED at index " << i); }
    std::cout << "SUCCESS! std::mt19937 matches original MT19937 for " << N << " values\n";
}

void test_XMT19937_64()
{
    startTest("STL mt19937_64 / XMT19937_64 equivalence");
    constexpr unsigned long long seed = 5489ULL;
    constexpr size_t N = 10000;
    init_genrand64(seed);
    std::vector<uint64_t> ref(N);
    for (size_t i = 0; i < N; ++i) ref[i] = genrand64_int64();
    std::mt19937_64 stl(seed);
    for (size_t i = 0; i < N; ++i) { uint64_t v = stl(); MYASSERT(v == ref[i], "FAILED at index " << i); }
    std::cout << "SUCCESS! std::mt19937_64 matches original MT19937-64 for " << N << " values\n";
    XMT19937_64<> xmt;
    xmt.reinit(seed, (const typename XMT19937_64<>::poly_t*)nullptr, nullptr);
    for (size_t i = 0; i < N; ++i) { uint64_t v = xmt.genrand_uint64(); MYASSERT(v == ref[i], "FAILED at index " << i); }
    std::cout << "SUCCESS! XMT19937_64 matches original MT19937-64 for " << N << " values\n";
}

void testEncoding()
{
    startTest("encoding");
    encodingTests<19937, 19937>();
    encodingTests<19937, 1007>();
    encodingTests<1007, 19937>();
    encodingTests<1007, 1007>();
}

void testSquareMatrix()
{
    startTest("square matrix calculation");
    squareTests(std::index_sequence<1, 5, 8, 13, 16, 20, 28, 32, 36, 60, 64, 68, 85, 126, 128, 150>{});
}

template <typename M>
struct JumpMatrix
{
    JumpMatrix() : nSteps(0), exponent(-1) {}
    JumpMatrix(const M* m, size_t nSteps, int exponent = -1) : p(m), nSteps(nSteps), exponent(exponent) {}
    std::unique_ptr<const M> p;
    size_t nSteps;
    int exponent;
};

template <GenType G, size_t L, size_t I, QryMode QM, typename M>
void testEquivalence(size_t nCommonJumpRepeat, const JumpMatrix<M>& commonJump, const JumpMatrix<M>& seqJump)
{
    if constexpr (I <= L && I <= SIMD_N_BITS) {
        using Gen = typename GenTraits<G, L, I, QM>::gen_t;
        using output_word_t = typename Gen::output_word_t;
        using poly_t = typename Gen::poly_t;

        const size_t commonJumpSteps = commonJump.p ? commonJump.nSteps : 0;
        const size_t sequenceJumpSteps = seqJump.p ? seqJump.nSteps : 0;
        const int commonJumpExp = commonJump.p ? commonJump.exponent : -1;
        const int sequenceJumpExp = seqJump.p ? seqJump.exponent : -1;

        constexpr size_t VecLen = Gen::s_regLenBits;
        constexpr QryMode QryMode = QM;
        constexpr size_t s_nStates = Gen::s_nStates;
        constexpr size_t s_n32InOneWord = Gen::s_n32InOneWord;
        constexpr size_t s_nWordInOneWord = s_n32InOneWord * sizeof(uint32_t) / sizeof(output_word_t);

        size_t blkSize;
        switch (QryMode) {
            case QM_Any: blkSize = 0; break;
            case QM_Scalar: blkSize = 1; break;
            case QM_Block16: blkSize = 64 / sizeof(output_word_t); break;
            default: THROW("how did we get here?");
        }

        std::cout << g_genName[G] << "< " << std::setw(3) << VecLen << ", " << std::setw(7) << queryModeName(QryMode) << ", " << std::setw(3) << Gen::s_regLenBitsHw << ">"
            << ", common jump of " << std::setw(4) << commonJumpSteps << " repeated " << nCommonJumpRepeat << " times, sequence jump of " << std::setw(4) << sequenceJumpSteps
            << ", block size " << std::setw(4) << (blkSize > 0 ? std::to_string(blkSize) : "rand") << " ... " << std::flush;

        const size_t nTest = sizeof(output_word_t) == 4 ? g_nRandomTest : g_nRandomTest64;
        std::vector<output_word_t> aligneddst(nTest);
        const M* jumpMat = (s_nStates > 1) ? seqJump.p.get() : nullptr;
        Gen mtPoly;
        Gen mtMatrix;
        bool polyOk = false;
        poly_t cp, sp;

        if (nCommonJumpRepeat <= 1) {
            std::string poly_dir, gen_tag;
            if (G == VSFMT) { poly_dir = "./dat/poly/sfmt/"; gen_tag = "sfmt"; }
            else if (G == VMT || G == XMT) { poly_dir = "./dat/poly/mt32/"; gen_tag = "mt32"; }
            else { poly_dir = "./dat/poly/mt64/"; gen_tag = "mt64"; }

            std::string cPolyFile, sPolyFile;
            if (commonJumpExp >= 0) { std::stringstream ss; ss << poly_dir << "J" << std::setw(5) << std::setfill('0') << commonJumpExp << "." << gen_tag << ".bits"; cPolyFile = ss.str(); }
            if (sequenceJumpExp >= 0) { std::stringstream ss; ss << poly_dir << "J" << std::setw(5) << std::setfill('0') << sequenceJumpExp << "." << gen_tag << ".bits"; sPolyFile = ss.str(); }

            bool cpExists = false, spExists = false;
            if (!cPolyFile.empty()) { std::ifstream is(cPolyFile, std::ios::binary); if (is) { cp.fromBin(is); cpExists = true; } }
            if (!sPolyFile.empty()) { std::ifstream is(sPolyFile, std::ios::binary); if (is) { sp.fromBin(is); spExists = true; } }

            if ((commonJumpExp < 0 || cpExists) && (sequenceJumpExp < 0 || spExists)) {
                if constexpr (sizeof(output_word_t) == 4) mtPoly.reinit(g_seedinit, g_seedlength, cpExists ? &cp : nullptr, spExists ? &sp : nullptr);
                else mtPoly.reinit((output_word_t)5489ULL, cpExists ? &cp : nullptr, spExists ? &sp : nullptr);
                polyOk = true;
            }
        }

        // Always initialize Matrix version
        if constexpr (sizeof(output_word_t) == 4) mtMatrix.reinit(g_seedinit, g_seedlength, nCommonJumpRepeat, commonJump.p.get(), jumpMat);
        else mtMatrix.reinit((output_word_t)5489ULL, nCommonJumpRepeat, commonJump.p.get(), jumpMat);

        auto runCheck = [&](Gen& mt, const char* label) {
            output_word_t* pdst = aligneddst.data();
            if constexpr (QryMode != QM_Any) {
                for (size_t i = 0; i < (size_t)(nTest / blkSize); ++i)
                    if constexpr (QryMode == QM_Scalar) {
                        if constexpr (sizeof(output_word_t) == 4) *pdst++ = mt.genrand_uint32(); else *pdst++ = mt.genrand_uint64();
                    } else {
                        if constexpr (sizeof(output_word_t) == 4) mt.genrand_uint32_blk16(pdst); else mt.genrand_word_blk(pdst);
                        pdst += blkSize;
                    }
            } else {
                size_t n = nTest;
                while (n) {
                    size_t sz = std::min(n, (size_t)(rand() % 100 + 1));
                    if constexpr (sizeof(output_word_t) == 4) mt.genrand_uint32_anySize(pdst, sz); else mt.genrand_word_anySize(pdst, sz);
                    pdst += sz; n -= sz;
                }
            }

            for (size_t i = 0; i < nTest; ++i) {
                output_word_t r2 = aligneddst[i];
                size_t genIndex = (i % (s_nWordInOneWord * s_nStates)) / s_nWordInOneWord;
                size_t seqIndex = (i % s_nWordInOneWord) + (i / (s_nWordInOneWord * s_nStates)) * s_nWordInOneWord;
                size_t benchmarkindex = seqIndex + commonJumpSteps * nCommonJumpRepeat + sequenceJumpSteps * genIndex;
                if constexpr (sizeof(output_word_t) == 4) {
                    if (g_benchmark[benchmarkindex] != r2) {
                        std::cout << "Mismatch at i=" << i << ", genIndex=" << genIndex << ", seqIndex=" << seqIndex << "\n";
                        std::cout << "  benchmarkindex=" << benchmarkindex << ", expected=" << std::hex << g_benchmark[benchmarkindex] << ", got=" << r2 << std::dec << "\n";
                        std::cout << "  Context (g_benchmark): ";
                        for (int k = -2; k <= 2; ++k) if ((int)benchmarkindex + k >= 0) std::cout << std::hex << g_benchmark[benchmarkindex + k] << " ";
                        std::cout << std::dec << std::endl;
                    }
                    MYASSERT(g_benchmark[benchmarkindex] == r2, label << " FAILED!");
                } else {
                    MYASSERT(g_benchmark64[benchmarkindex] == r2, label << " FAILED!");
                }
            }
        };

        if (polyOk) runCheck(mtPoly, "Polynomial");
        runCheck(mtMatrix, "Matrix");
        std::cout << (polyOk ? "BOTH SUCCESS!\n" : "MATRIX SUCCESS!\n");
    }
}

template <GenType G, size_t L, size_t I, QryMode QM, typename M>
void equivalenceTests3(const JumpMatrix<M>& jumpSmall, const JumpMatrix<M>& jumpBig)
{
    JumpMatrix<M> noJump{};
    if constexpr (QM == QM_Scalar) {
        testEquivalence<G, L, I, QM>(0, noJump, noJump);
        testEquivalence<G, L, I, QM>(1, jumpSmall, noJump);
        testEquivalence<G, L, I, QM>(1, jumpBig, noJump);
    }
}

template <GenType G, size_t L, size_t I, QryMode...QMs, typename M> void equivalenceTests2(const JumpMatrix<M>& js, const JumpMatrix<M>& jb) { (equivalenceTests3<G, L, I, QMs>(js, jb), ...); }
template <GenType G, size_t L, size_t...Is, typename M> void equivalenceTests1(const JumpMatrix<M>& js, const JumpMatrix<M>& jb) { (equivalenceTests2<G, L, Is, QM_Scalar, QM_Block16, QM_Any>(js, jb), ...); }
template <GenType G, size_t...Ls, typename M> void equivalenceTests0(const JumpMatrix<M>& js, const JumpMatrix<M>& jb) {
    if constexpr (G == XMT || G == XMT64) (equivalenceTests1<G, Ls, Ls>(js, jb), ...);
    else if constexpr (G == VMT64) (equivalenceTests1<G, Ls, 128, 256, 512>(js, jb), ...);
    else (equivalenceTests1<G, Ls, 32, 128, 256, 512>(js, jb), ...);
}

// Direct polynomial-vs-matrix comparison for a single generator instantiation.
// No g_benchmark is needed: the two jump methods are compared against each other.
// Useful for large exponents (e.g. 2^100, 2^19933) where the absolute stream
// position cannot be stored in a pre-computed g_benchmark array.
template <GenType G, size_t L, size_t I, typename M>
void testPolyVsMatrix(int exp, const M* matPtr)
{
    if constexpr (I <= L && I <= SIMD_N_BITS) {
        using Gen = typename GenTraits<G, L, I, QM_Scalar>::gen_t;
        using output_word_t = typename Gen::output_word_t;
        using poly_t = typename Gen::poly_t;
        constexpr size_t VecLen = Gen::s_regLenBits;

        std::string poly_dir, gen_tag;
        if (G == VSFMT) { poly_dir = "./dat/poly/sfmt/"; gen_tag = "sfmt"; }
        else if (G == VMT || G == XMT) { poly_dir = "./dat/poly/mt32/"; gen_tag = "mt32"; }
        else { poly_dir = "./dat/poly/mt64/"; gen_tag = "mt64"; }

        std::stringstream ss;
        ss << poly_dir << "J" << std::setw(5) << std::setfill('0') << exp << "." << gen_tag << ".bits";
        std::string polyFile = ss.str();

        std::cout << g_genName[G] << "< " << std::setw(3) << VecLen << ", " << std::setw(7) << "Scalar" << ", " << std::setw(3) << Gen::s_regLenBitsHw << ">"
            << ", poly==matrix 2^" << exp << " ... " << std::flush;

        poly_t poly;
        {
            std::ifstream is(polyFile, std::ios::binary);
            if (!is) { std::cout << "SKIPPED (no J" << exp << " file)\n"; return; }
            poly.fromBin(is);
        }

        const size_t nTest = sizeof(output_word_t) == 4 ? g_nRandomTest : g_nRandomTest64;
        Gen genPoly, genMatrix;
        if constexpr (sizeof(output_word_t) == 4) {
            genPoly.reinit(g_seedinit, g_seedlength, &poly, nullptr);
            genMatrix.reinit(g_seedinit, g_seedlength, 1, matPtr, nullptr);
        } else {
            genPoly.reinit((output_word_t)5489ULL, &poly, nullptr);
            genMatrix.reinit((output_word_t)5489ULL, 1, matPtr, nullptr);
        }

        for (size_t i = 0; i < nTest; ++i) {
            output_word_t a, b;
            if constexpr (sizeof(output_word_t) == 4) { a = genPoly.genrand_uint32(); b = genMatrix.genrand_uint32(); }
            else { a = genPoly.genrand_uint64(); b = genMatrix.genrand_uint64(); }
            if (a != b) {
                std::cout << "Mismatch at i=" << i << ", poly=" << std::hex << a << ", matrix=" << b << std::dec << "\n";
                MYASSERT(a == b, "Poly vs Matrix FAILED!");
            }
        }
        std::cout << "SUCCESS!\n";
    }
}

template <GenType G, size_t L, size_t...Is, typename M>
void polyVsMatrixTests2(int exp, const M* matPtr) { (testPolyVsMatrix<G, L, Is>(exp, matPtr), ...); }

template <GenType G, size_t...Ls, typename M>
void polyVsMatrixTests0(int exp, const M* matPtr) {
    if constexpr (G == XMT || G == XMT64) (polyVsMatrixTests2<G, Ls, Ls>(exp, matPtr), ...);
    else if constexpr (G == VMT64) (polyVsMatrixTests2<G, Ls, 128, 256, 512>(exp, matPtr), ...);
    else (polyVsMatrixTests2<G, Ls, 32, 128, 256, 512>(exp, matPtr), ...);
}

void test_XVMT19937()
{
    generateBenchmark_MT19937();
    typedef MT19937Matrix<32> matrix_t;
    typedef JumpMatrix<matrix_t> pmatrix_t;
    pmatrix_t noJump;
    pmatrix_t jumpMatrix1(new matrix_t, 1, 0);
    pmatrix_t jumpMatrix512(new matrix_t(std::string("./dat/matrix/mt32/F00009.mt32.bits")), 512, 9);
    pmatrix_t jumpMatrixPeriod(new matrix_t(std::string("./dat/matrix/mt32/F19937.mt32.bits")), 1, 19937);
    startTest(g_genName[VMT]);
    equivalenceTests0<VMT, 32, 128, 256, 512>(jumpMatrix1, jumpMatrix512);
    std::cout << "VMT19937: a jump of size 2^19937 is equivalent to a jump of size 1\n";
    testEquivalence<VMT, 128, 128, QM_Scalar>(1, jumpMatrixPeriod, noJump);
    startTest(g_genName[XMT]);
    equivalenceTests0<XMT, 32, 128, 256, 512>(jumpMatrix1, jumpMatrix512);
    for (int exp : {9, 100, 19933, 19934, 19935, 19936}) {
        std::stringstream sf; sf << "./dat/matrix/mt32/F" << std::setw(5) << std::setfill('0') << exp << ".mt32.bits";
        matrix_t mat(sf.str());
        polyVsMatrixTests0<VMT, 32, 128, 256, 512>(exp, &mat);
        polyVsMatrixTests0<XMT, 32, 128, 256, 512>(exp, &mat);
    }
}

void test_VSFMT19937()
{
    startTest(g_genName[VSFMT]);
    generateBenchmark_SFMT19937();
    typedef SFMT19937Matrix matrix_t;
    typedef JumpMatrix<matrix_t> pmatrix_t;
    pmatrix_t noJump;
    pmatrix_t jumpMatrix4(new matrix_t, 4, 0);
    pmatrix_t jumpMatrix512(new matrix_t(std::string("./dat/matrix/sfmt/F00009.sfmt.bits")), 512, 7);
    equivalenceTests0<VSFMT, 128, 256, 512>(jumpMatrix4, jumpMatrix512);
    for (int exp : {9, 100, 19933, 19934, 19935, 19936}) {
        std::stringstream sf; sf << "./dat/matrix/sfmt/F" << std::setw(5) << std::setfill('0') << exp << ".sfmt.bits";
        matrix_t mat(sf.str());
        polyVsMatrixTests0<VSFMT, 128, 256, 512>(exp, &mat);
    }
}

void test_VMT19937_64()
{
    generateBenchmark_MT19937_64();
    typedef MT19937Matrix<64> matrix_t;
    typedef JumpMatrix<matrix_t> pmatrix_t;
    pmatrix_t noJump;
    pmatrix_t jumpMatrix1(new matrix_t, 1, 0);
    pmatrix_t jumpMatrix512(new matrix_t(std::string("./dat/matrix/mt64/F00009.mt64.bits")), 512, 9);
    pmatrix_t jumpMatrixPeriod(new matrix_t(std::string("./dat/matrix/mt64/F19937.mt64.bits")), 1, 19937);
    startTest(g_genName[VMT64]);
    equivalenceTests0<VMT64, 128, 256, 512>(jumpMatrix1, jumpMatrix512);
    equivalenceTests1<VMT64, 64, 64>(jumpMatrix1, jumpMatrix512);
    startTest(g_genName[XMT64]);
    equivalenceTests0<XMT64, 64, 128, 256, 512>(jumpMatrix1, jumpMatrix512);
    for (int exp : {9, 100, 19933, 19934, 19935, 19936}) {
        std::stringstream sf; sf << "./dat/matrix/mt64/F" << std::setw(5) << std::setfill('0') << exp << ".mt64.bits";
        matrix_t mat(sf.str());
        polyVsMatrixTests0<VMT64, 128, 256, 512>(exp, &mat);
        polyVsMatrixTests0<XMT64, 64, 128, 256, 512>(exp, &mat);
    }
}

template <typename T> void printReg(std::string&& name, T v) {
    constexpr unsigned n = sizeof(T); alignas(sizeof(T)) unsigned char bytes[n]; std::copy_n((unsigned char*)&v, n, bytes);
    std::cout << std::setw(3) << name << ": "; for (unsigned i = 0; i < n; ++i) std::cout << std::setw(3) << (int)bytes[i]; std::cout << '\n';
}

template <size_t n32, typename T> void testAlignR32(const unsigned char *data, T a, T b) {
    auto c = T::template alignr32<n32>(a, b); const unsigned char* got = (const unsigned char*) &c;
    for (size_t i = 0; i < sizeof(T); ++i) { size_t srcIndex = i + n32 * 4; MYASSERT(got[i] == data[srcIndex], "error in alignR32"); }
}

template <size_t...n32s> void testSimdAlignR32(std::index_sequence<n32s...>&&) {
    constexpr size_t n32 = (sizeof...(n32s) - 1); constexpr size_t nBits = n32 * 32;
    using T = xvmt::details::SimdRegister<nBits, BitLenToIsa<nBits>::s_isa>;
    alignas(64) unsigned char data[128]; std::iota(data, data + 128, 0);
    T v0((const T*)data); T v1(((const T*)data)+1); (testAlignR32<n32s>(data, v0, v1), ...);
}

void test_SIMD_special_methods() {
    using XV = xvmt::details::SimdRegister<128, BitLenToIsa<128>::s_isa>;
    XV a(1, 2, 3, 4); XV b(5, 6, 7, 8);
    MYASSERT(XV::alignr32<1>(a, b).eq(XV(2, 3, 4, 5)), "alignr32 failed");
    std::cout << "SIMD special methods tests passed!\n";
}

int main(int argc, const char** argv)
{
    ArgMap args = parseArgs(argc, argv);
    try {
        test_SIMD_special_methods();
        testEncoding();
        testSquareMatrix();
        test_STL_MT19937();
        test_XMT19937_64();
        test_XVMT19937();
        test_VSFMT19937();
        test_VMT19937_64();
    }
    catch (const std::exception& e) { std::cout << e.what() << "\n"; return -1; }
    return 0;
}
