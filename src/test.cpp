#define TESTING
#define SIMD_EMULATION

#include "VMT19937.h"

const uint32_t seedlength = 4;
const uint32_t seedinit[seedlength] = { 0x123, 0x234, 0x345, 0x456 };

const uint64_t nRandomTest = 50ul * 624 * 16;
//const uint64_t nRandomTest = 624 * 1;

extern "C" unsigned long genrand_int32();
extern "C" void init_by_array(unsigned long init_key[], int key_length);


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
    m3.square(m, buffers, nullptr);

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

template <template <size_t, VMT19937QueryMode> class Gen, size_t VecLen, VMT19937QueryMode BlkMode>
void testEquivalence(size_t mCommonJumpRepeat, const MT19937Matrix* commonJump, const MT19937Matrix* seqJump, size_t commonJumpSize, size_t sequenceJumpSize)
{
    typedef Gen<VecLen, BlkMode> gen_t;
    static const size_t BlkSize = gen_t::s_qryBlkSize;

    std::cout << "Testing equivalence of generators with SIMD length " << VecLen
        << ", common jump ahead of " << commonJumpSize << " repeated " << mCommonJumpRepeat << " times, sequence jump size of " << sequenceJumpSize
        << ", block size " << BlkSize << " ... ";

    const static size_t s_nStates = gen_t::s_nStates;

    AlignedVector<uint32_t, 64> aligneddst(nRandomTest);

    gen_t mt(seedinit, seedlength, mCommonJumpRepeat, commonJump, seqJump);
    for (size_t i = 0; i < nRandomTest / BlkSize; ++i)
        if constexpr (BlkMode == QM_Scalar)
            aligneddst[i] = mt.genrand_uint32();
        else if constexpr (BlkMode == QM_Block16)
            mt.genrand_uint32_blk16(aligneddst.data() + i * BlkSize);
        else if constexpr (BlkMode == QM_StateSize)
            mt.genrand_uint32_stateBlk(aligneddst.data() + i * BlkSize);
        else
            NOT_IMPLEMENTED;

    for (size_t i = 0; i < nRandomTest; ++i) {
        uint32_t r2 = aligneddst[i];
        size_t seqIndex = i / s_nStates;
        size_t genIndex = i % s_nStates;
        size_t benchmarkindex = seqIndex + commonJumpSize* mCommonJumpRepeat + sequenceJumpSize * genIndex;
        MYASSERT(benchmark[benchmarkindex] == r2, "FAILED!\n"
                << "Difference found: out[" << i << "] = " << r2
                << ", benchmark[" << benchmarkindex  << "] = " << benchmark[benchmarkindex]);
    }

    std::cout << "SUCCESS!\n";
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

void test_VMT19937()
{
    startTest("VMT19937");

    generateBenchmark_MT19937();

    MT19937Matrix jumpMatrix1;                                          // jump ahead 1 element
    MT19937Matrix jumpMatrix512(std::string("./dat/F00009.bits"));      // jump ahead 2^9 (512) elements
    MT19937Matrix jumpMatrix1024(std::string("./dat/F00010.bits"));     // jump ahead 2^10 (1024) elements
    MT19937Matrix jumpMatrixPeriod(std::string("./dat/F19937.bits"));   // jump ahead 2^19937 elements

    testEquivalence<VMT19937, 32, QM_Scalar>(0, nullptr, nullptr, 0, 0);
    testEquivalence<VMT19937, 32, QM_Scalar>(1, &jumpMatrix1024, nullptr, 1024, 0);
    // two jumps of 512 are equivalent to one jump of 1024
    testEquivalence<VMT19937, 32, QM_Scalar>(2, &jumpMatrix512, nullptr, 512, 0);
    // since the period is 2^19937-1, after applying a jump matrix of 2^19937, we restart from the sequence from step 1
    testEquivalence<VMT19937, 32, QM_Scalar>(1, &jumpMatrixPeriod, nullptr, 1, 0);

    testEquivalence<VMT19937, 64, QM_Scalar>(0, nullptr, nullptr, 0, 0);
    testEquivalence<VMT19937, 128, QM_Scalar>(0, nullptr, nullptr, 0, 0);
    testEquivalence<VMT19937, 128, QM_Scalar>(0, nullptr, &jumpMatrix1, 0, 1);
    testEquivalence<VMT19937, 128, QM_Scalar>(0, nullptr, &jumpMatrix1024, 0, 1024);
    testEquivalence<VMT19937, 128, QM_Block16>(0, nullptr, nullptr, 0, 0);
    testEquivalence<VMT19937, 128, QM_Block16>(0, nullptr, &jumpMatrix1, 0, 1);
    testEquivalence<VMT19937, 128, QM_Block16>(0, nullptr, &jumpMatrix1024, 0, 1024);
    testEquivalence<VMT19937, 128, QM_StateSize>(0, nullptr, nullptr, 0, 0);
    testEquivalence<VMT19937, 128, QM_StateSize>(0, nullptr, &jumpMatrix1, 0, 1);
    testEquivalence<VMT19937, 128, QM_StateSize>(0, nullptr, &jumpMatrix1024, 0, 1024);

    testEquivalence<VMT19937, 256, QM_Scalar>(0, nullptr, nullptr, 0, 0);
    testEquivalence<VMT19937, 256, QM_Scalar>(0, nullptr, &jumpMatrix1, 0, 1);
    testEquivalence<VMT19937, 256, QM_Scalar>(0, nullptr, &jumpMatrix1024, 0, 1024);
    testEquivalence<VMT19937, 256, QM_Block16>(0, nullptr, nullptr, 0, 0);
    testEquivalence<VMT19937, 256, QM_Block16>(0, nullptr, &jumpMatrix1, 0, 1);
    testEquivalence<VMT19937, 256, QM_Block16>(0, nullptr, &jumpMatrix1024, 0, 1024);
    testEquivalence<VMT19937, 256, QM_StateSize>(0, nullptr, nullptr, 0, 0);
    testEquivalence<VMT19937, 256, QM_StateSize>(0, nullptr, &jumpMatrix1, 0, 1);
    testEquivalence<VMT19937, 256, QM_StateSize>(0, nullptr, &jumpMatrix1024, 0, 1024);

    testEquivalence<VMT19937, 512, QM_Scalar>(0, nullptr, nullptr, 0, 0);
    testEquivalence<VMT19937, 512, QM_Scalar>(0, nullptr, &jumpMatrix1, 0, 1);
    testEquivalence<VMT19937, 512, QM_Scalar>(0, nullptr, &jumpMatrix1024, 0, 1024);
    testEquivalence<VMT19937, 512, QM_Block16>(0, nullptr, nullptr, 0, 0);
    testEquivalence<VMT19937, 512, QM_Block16>(0, nullptr, &jumpMatrix1, 0, 1);
    testEquivalence<VMT19937, 512, QM_Block16>(0, nullptr, &jumpMatrix1024, 0, 1024);
    testEquivalence<VMT19937, 512, QM_StateSize>(0, nullptr, nullptr, 0, 0);
    testEquivalence<VMT19937, 512, QM_StateSize>(0, nullptr, &jumpMatrix1, 0, 1);
    testEquivalence<VMT19937, 512, QM_StateSize>(0, nullptr, &jumpMatrix1024, 0, 1024);
}

int main()
{
    try {
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
