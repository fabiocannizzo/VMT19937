#define TESTING
#define SIMD_EMULATION

#include "jump_matrix.h"

#include "MSMT19937.h"

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

template <size_t VecLen, MSMT19937QueryMode BlkMode>
void testEquivalence(const MT19937Matrix* commonJump, const MT19937Matrix* seqJump, size_t commonJumpSize, size_t sequenceJumpSize)
{
    typedef MSMT19937<VecLen, BlkMode> gen_t;
    static const size_t BlkSize = gen_t::s_qryBlkSize;

    std::cout << "Testing equivalence of generators with SIMD length " << VecLen
        << " and common jump ahead of " << commonJumpSize << " and sequence jump size of " << sequenceJumpSize
        << " and block size " << BlkSize << " ... ";

    const static size_t s_M = VecLen / 32;

    std::vector<uint32_t> dst(nRandomTest + 64 / sizeof(uint32_t));
    uint32_t* aligneddst = (uint32_t*)((intptr_t)dst.data() + (64 - ((intptr_t)dst.data() % 64)));

    gen_t mt(seedinit, seedlength, commonJump, seqJump);
    for (size_t i = 0; i < nRandomTest / BlkSize; ++i)
        if constexpr (BlkMode == QM_Scalar)
            aligneddst[i] = mt.genrand_uint32();
        else if constexpr (BlkMode == QM_Block16)
            mt.genrand_uint32_blk16(aligneddst + i * BlkSize);
        else if constexpr (BlkMode == QM_StateSize)
            mt.genrand_uint32_stateBlk(aligneddst + i * BlkSize);
        else
            NOT_IMPLEMENTED;

    for (size_t i = 0; i < nRandomTest; ++i) {
        uint32_t r2 = aligneddst[i];
        size_t seqIndex = i / s_M;
        size_t genIndex = i % s_M;
        size_t benchmarkindex = seqIndex + commonJumpSize + sequenceJumpSize * genIndex;
        MYASSERT(benchmark[benchmarkindex] == r2, "FAILED!\n"
                << "Difference found: out[" << i << "] = " << r2
                << ", benchmark[" << benchmarkindex  << "] = " << benchmark[benchmarkindex]);
    }

    std::cout << "SUCCESS!\n";
}

void generateBenchmark()
{
    unsigned long init[seedlength];
    for (size_t i = 0; i < seedlength; ++i)
        init[i] = seedinit[i];

    std::cout << "Generate random numbers with the original C source code ... ";
    init_by_array(init, seedlength);
    for (size_t i = 0, n  = benchmark.size(); i < n; ++i)
        benchmark[i] = (uint32_t)genrand_int32();
    std::cout << "done!\n";
    printSome(benchmark);
}


int main()
{
    try {
        encodingTests<19937, 19937>();
        encodingTests<19937, 1007>();
        encodingTests<1007, 19937>();
        encodingTests<1007, 1007>();

        squareTests(std::index_sequence<1, 5, 8, 13, 16, 20, 28, 32, 36, 60, 64, 68, 85, 126, 128, 150>{});

        generateBenchmark();

        MT19937Matrix jumpMatrix1;
        MT19937Matrix jumpMatrix1024("./dat/F00010.bits");
        MT19937Matrix jumpMatrixPeriod("./dat/F19937.bits");

        testEquivalence<32, QM_Scalar>(nullptr, nullptr, 0, 0);
        testEquivalence<32, QM_Scalar>(&jumpMatrix1024, nullptr, 1024, 0);
        // since the period is 2^19937-1, after applying a jump matrix of 2^19937, we restart from the sequence from step 1
        testEquivalence<32, QM_Scalar>(&jumpMatrixPeriod, nullptr, 1, 0);
        testEquivalence<64, QM_Scalar>(nullptr, nullptr, 0, 0);
        testEquivalence<128, QM_Scalar>(nullptr, nullptr, 0, 0);
        testEquivalence<128, QM_Scalar>(nullptr, &jumpMatrix1, 0, 1);
        testEquivalence<128, QM_Scalar>(nullptr, &jumpMatrix1024, 0, 1024);
        testEquivalence<128, QM_Block16>(nullptr, nullptr, 0, 0);
        testEquivalence<128, QM_Block16>(nullptr, &jumpMatrix1, 0, 1);
        testEquivalence<128, QM_Block16>(nullptr, &jumpMatrix1024, 0, 1024);
        testEquivalence<128, QM_StateSize>(nullptr, nullptr, 0, 0);
        testEquivalence<128, QM_StateSize>(nullptr, &jumpMatrix1, 0, 1);
        testEquivalence<128, QM_StateSize>(nullptr, &jumpMatrix1024, 0, 1024);

        testEquivalence<256, QM_Scalar>(nullptr, nullptr, 0, 0);
        testEquivalence<256, QM_Scalar>(nullptr, &jumpMatrix1, 0, 1);
        testEquivalence<256, QM_Scalar>(nullptr, &jumpMatrix1024, 0, 1024);
        testEquivalence<256, QM_Block16>(nullptr, nullptr, 0, 0);
        testEquivalence<256, QM_Block16>(nullptr, &jumpMatrix1, 0, 1);
        testEquivalence<256, QM_Block16>(nullptr, &jumpMatrix1024, 0, 1024);
        testEquivalence<256, QM_StateSize>(nullptr, nullptr, 0, 0);
        testEquivalence<256, QM_StateSize>(nullptr, &jumpMatrix1, 0, 1);
        testEquivalence<256, QM_StateSize>(nullptr, &jumpMatrix1024, 0, 1024);

        testEquivalence<512, QM_Scalar>(nullptr, nullptr, 0, 0);
        testEquivalence<512, QM_Scalar>(nullptr, &jumpMatrix1, 0, 1);
        testEquivalence<512, QM_Scalar>(nullptr, &jumpMatrix1024, 0, 1024);
        testEquivalence<512, QM_Block16>(nullptr, nullptr, 0, 0);
        testEquivalence<512, QM_Block16>(nullptr, &jumpMatrix1, 0, 1);
        testEquivalence<512, QM_Block16>(nullptr, &jumpMatrix1024, 0, 1024);
        testEquivalence<512, QM_StateSize>(nullptr, nullptr, 0, 0);
        testEquivalence<512, QM_StateSize>(nullptr, &jumpMatrix1, 0, 1);
        testEquivalence<512, QM_StateSize>(nullptr, &jumpMatrix1024, 0, 1024);
    }
    catch (const std::exception& e) {
        std::cout << e.what() << "\n";
        return -1;
    }

    return 0;
}
