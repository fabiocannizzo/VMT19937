#define TESTING
#define SIMD_EMULATION

#include "VRandGen.h"

#include "../SFMT-src-1.5.1/SFMT.h"

const uint32_t seedlength = 4;
const uint32_t seedinit[seedlength] = { 0x123, 0x234, 0x345, 0x456 };

const uint64_t nRandomTest = 50ul * 624 * 16;

extern "C" unsigned long genrand_int32();
extern "C" void init_by_array(unsigned long init_key[], int key_length);

template <typename T>
struct GenTraits;

template <size_t VecLen, VRandGenQueryMode QryMode>
struct GenTraits<VMT19937<VecLen, QryMode>>
{
    static const char* name() { return "VMT19937"; }
};

template <size_t VecLen, VRandGenQueryMode QryMode>
struct GenTraits<VSFMT19937<VecLen, QryMode>>
{
    static const char* name() { return "VSFMT19937"; }
};

std::vector<uint32_t> benchmark(nRandomTest + 10000);

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
void testEncoder(const BinaryMatrix<nRows, nCols>& m, EncodeMode enc)
{
    const char* modename = enc == Base64 ? "base64" : "hex";
    
    BinaryMatrix<nRows, nCols> m2;

    std::cout << "saving matrix to " << modename << " stream\n";
    std::ostringstream os;
    if (enc == Base64)
        m.toBase64(os);
    else
        m.toHex(os);

    std::cout << "first 32 characters of the stream\n";
    std::string s = os.str();
    for (size_t i = 0; i < 32; ++i)
        std::cout << s[i];
    std::cout << "\n";

    std::cout << "reading back the matrix from " << modename << " stream\n";
    std::istringstream is(os.str());
    if (enc == Base64)
        m2.fromBase64(is);
    else
        m2.fromHex(is);

    std::cout << "compare with original matrix\n";
    MYASSERT((m == m2), "error in roundtrip");

    std::cout << "completed\n";
}

template <size_t NBITS>
void testSquare(const BinarySquareMatrix<NBITS>& m)
{
    BinarySquareMatrix<NBITS> m2, m3;

    // slow bit by bit multiplication
    //std::cout << "compute matrix multiplication the classical way\n";
    for (size_t r = 0; r < NBITS; ++r) {
        //std::cout << r << "\n";
        for (size_t c = 0; c < NBITS; ++c) {
            size_t s = 0;
            for (size_t k = 0; k < NBITS; ++k) {
                s ^= m.getBit(r, k) && m.getBit(k, c);
            }
            if (s)
                m2.setBit(r, c);
        }
    }

    const size_t nThreads = 4;
    //std::cout << "compute matrix multiplication vectorially\n";
    std::vector<typename BinarySquareMatrix<NBITS>::buffer_t> buffers(nThreads);
    m3.square(m, buffers);

    MYASSERT((m2 == m3), "error in square");

    //std::cout << "SUCCESS\n";
}

template <size_t nRows, size_t nCols>
void encodingTests()
{
    BinaryMatrix<nRows, nCols> m;
    m.initRand();
    std::cout << "\ngenerated random matrix with size (" << m.s_nBitRows << "x" << m.s_nBitCols << ") with " << m.nnz() << " non zero elements\n";
    m.printBits(0, 0, 10, 32);

    testEncoder(m, Base64);
    testEncoder(m, Hex);
}

template <size_t NBits>
void squareTest()
{
    std::cout << "testing multiplication with matrices of size: " << NBits << "\n";
    BinarySquareMatrix<NBits> m;
    for (size_t i = 0; i < 10; ++i) {
        m.resetZero();
        m.initRand();
        testSquare(m);
    }
}

template <size_t...NBits>
void squareTests(std::index_sequence<NBits...>&&)
{
    (squareTest<NBits>(), ...);
}


void generateBenchmark_MT19937()
{
    unsigned long init[seedlength];
    for (size_t i = 0; i < seedlength; ++i)
        init[i] = seedinit[i];

    std::cout << "Generate MT19937 random numbers with the original C source code ... ";
    init_by_array(init, seedlength);
    for (size_t i = 0, n  = benchmark.size(); i < n; ++i)
        benchmark[i] = (uint32_t)genrand_int32();
    std::cout << "done!\n";
    printSome(benchmark);
}

void generateBenchmark_SFMT19937()
{
    std::cout << "Generate SFMT19937 random numbers with the original C source code ... ";
    sfmt_t sfmtgen;
    sfmt_init_by_array(&sfmtgen, const_cast<uint32_t *>(seedinit), seedlength);
    for (size_t i = 0, n = benchmark.size(); i < n; ++i)
        benchmark[i] = sfmt_genrand_uint32(&sfmtgen);
    std::cout << "done!\n";
    printSome(benchmark);
}

void startTest(const char* name)
{
    std::cout << "\n"
              << std::setw(40) << std::setfill('*') << "" << "\n" 
              << "Test " << name << "\n"
              << std::setw(40) << std::setfill('*') << "" << "\n\n"
              << std::setfill(' ');
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
    JumpMatrix() : jumpSize(0) {}
    JumpMatrix(const M* m, size_t jumpSize) : p(m), jumpSize(jumpSize) {}
    std::unique_ptr<const M> p;
    size_t jumpSize;  // jump size (i.e. number of elements skipped)
};

template <typename Gen, typename M>
void testEquivalence(size_t nCommonJumpRepeat, const JumpMatrix<M>& commonJump, const JumpMatrix<M>& seqJump)
{
    const size_t commonJumpSize = commonJump.p ? commonJump.jumpSize : 0;
    const size_t sequenceJumpSize = seqJump.p ? seqJump.jumpSize : 0;

    MYASSERT((!commonJump.p && nCommonJumpRepeat == 0) || (commonJump.p && nCommonJumpRepeat > 0), "(nCommonJumpRepeat>0) <=> commonJumpSize>0 ");

    const size_t VecLen = Gen::s_regLenBits;
    const VRandGenQueryMode QryMode = Gen::s_queryMode;
    const size_t BlkSize = QryMode == QM_Scalar ? 1 : QryMode == QM_Block16 ? 16 : Gen::s_n32InFullState;
    const size_t s_nStates = Gen::s_nStates;
    const size_t s_n32InOneWord = Gen::s_n32InOneWord;


    std::cout << GenTraits<Gen>::name() << ": Testing equivalence of generators with SIMD length " << VecLen
        << ", common jump ahead of " << commonJumpSize << " repeated " << nCommonJumpRepeat << " times, sequence jump size of " << sequenceJumpSize
        << ", block size " << BlkSize << " ... ";

    std::vector<uint32_t> aligneddst(nRandomTest);

    Gen mt(seedinit, seedlength, nCommonJumpRepeat, commonJump.p.get(), seqJump.p.get());
    for (size_t i = 0; i < nRandomTest / BlkSize; ++i)
        if constexpr (QryMode == QM_Scalar)
            aligneddst[i] = mt.genrand_uint32();
        else if constexpr (QryMode == QM_Block16)
            mt.genrand_uint32_blk16(aligneddst.data() + i * BlkSize);
        else if constexpr (QryMode == QM_StateSize)
            mt.genrand_uint32_stateBlk(aligneddst.data() + i * BlkSize);
        else
            NOT_IMPLEMENTED;

    for (size_t i = 0; i < nRandomTest; ++i) {
        uint32_t r2 = aligneddst[i];
        size_t genIndex = (i % (s_n32InOneWord * s_nStates)) / s_n32InOneWord;
        size_t seqIndex = (i % s_n32InOneWord) + (i / (s_n32InOneWord * s_nStates)) * s_n32InOneWord;
        size_t benchmarkindex = seqIndex + commonJumpSize * nCommonJumpRepeat + sequenceJumpSize * genIndex;
        MYASSERT(benchmark[benchmarkindex] == r2, "FAILED!\n"
                << "Difference found: out[" << i << "] = " << r2
                << ", benchmark[" << benchmarkindex  << "] = " << benchmark[benchmarkindex]);
    }

    std::cout << "SUCCESS!\n";
}

template <typename M, size_t L, VRandGenQueryMode QM>
struct GenFromMatrix;

template <size_t L, VRandGenQueryMode QM>
struct GenFromMatrix<MT19937Matrix, L, QM>
{
    typedef VMT19937<L, QM> gen_t;
};

template <size_t L, VRandGenQueryMode QM>
struct GenFromMatrix<SFMT19937Matrix, L, QM>
{
    typedef VSFMT19937<L, QM> gen_t;
};

template <size_t L, typename M>
void equivalenceTestsByMode(size_t nCommonJumpRepeat, const JumpMatrix<M>& commonJump, const JumpMatrix<M>& seqJump)
{
    testEquivalence<typename GenFromMatrix<M, L, QM_Scalar>::gen_t>(nCommonJumpRepeat, commonJump, seqJump);
    testEquivalence<typename GenFromMatrix<M, L, QM_Block16>::gen_t>(nCommonJumpRepeat, commonJump, seqJump);
    testEquivalence<typename GenFromMatrix<M, L, QM_StateSize>::gen_t>(nCommonJumpRepeat, commonJump, seqJump);
}

template <size_t L, bool testSeqJump, typename M>
void equivalenceTests(const JumpMatrix<M>& jumpSmall, const JumpMatrix<M>& jumpBig)
{
    JumpMatrix<M> noJump{};
    equivalenceTestsByMode<L>(0, noJump, noJump);
    equivalenceTestsByMode<L>(1, jumpSmall, noJump);
    equivalenceTestsByMode<L>(2, jumpSmall, noJump);
    equivalenceTestsByMode<L>(1, jumpBig, noJump);
    if constexpr (testSeqJump) {
        equivalenceTestsByMode<L>(1, jumpBig, jumpBig);
        equivalenceTestsByMode<L>(2, jumpBig, jumpBig);
        equivalenceTestsByMode<L>(0, noJump, jumpBig);
    }
}

void test_VMT19937()
{
    startTest("VMT19937");

    generateBenchmark_MT19937();

    typedef MT19937Matrix matrix_t;
    typedef JumpMatrix<matrix_t> pmatrix_t;

    pmatrix_t noJump;
    pmatrix_t jumpMatrix1(new matrix_t, 1);                                          // jump ahead 1 element
    pmatrix_t jumpMatrix512(new matrix_t(std::string("./dat/mt/F00009.bits")), 512);    // jump ahead 2^9 (512) elements
    pmatrix_t jumpMatrixPeriod(new matrix_t(std::string("./dat/mt/F19937.bits")), 1);   // jump ahead 2^19937 elements

    equivalenceTests<32, false>(jumpMatrix1, jumpMatrix512);
    equivalenceTests<128, true>(jumpMatrix1, jumpMatrix512);
    equivalenceTests<256, true>(jumpMatrix1, jumpMatrix512);
    equivalenceTests<512, true>(jumpMatrix1, jumpMatrix512);

    // since the period is 2^19937-1, after applying a jump matrix of 2^19937, we restart the sequence from step 1
    std::cout << "VMT19937: a jump of size 2^19937 is equivalent to a jump of size 1\n";
    testEquivalence<VMT19937<32, QM_Scalar>>(1, jumpMatrixPeriod, noJump);
}

void test_VSFMT19937()
{
    std::cout << "\nTest VSFMT19937\n";

    generateBenchmark_SFMT19937();

    typedef SFMT19937Matrix matrix_t;
    typedef JumpMatrix<matrix_t> pmatrix_t;

    pmatrix_t noJump;
    pmatrix_t jumpMatrix4(new matrix_t, 4);                                              // jump ahead 1 element
    pmatrix_t jumpMatrix512(new matrix_t(std::string("./dat/sfmt/F00009.bits")), 512);   // jump ahead 2^9 (1024) elements

    equivalenceTests<128, false>(jumpMatrix4, jumpMatrix512);
    equivalenceTests<256, true>(jumpMatrix4, jumpMatrix512);
    equivalenceTests<512, true>(jumpMatrix4, jumpMatrix512);
}

int main()
{
    try {
        test_VSFMT19937();
        testEncoding();
        testSquareMatrix();
        test_VMT19937();
    }
    catch (const std::exception& e) {
        std::cout << e.what() << "\n";
        return -1;
    }

    return 0;
}
