#define TESTING

#include "jump_matrix.h"

#include "../include/MT19937-SIMD.h"

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

    std::cout << "first 16 characters of the stream\n";
    std::string s = os.str();
    for (size_t i = 0; i < 16; ++i)
        std::cout << s[i];
    std::cout << "\n";

    std::cout << "reading back the matrix from " << modename << " stream\n";
    std::istringstream is(os.str());
    if (enc == Base64)
        m2.fromBase64(is);
    else
        m2.fromHex(is);

    std::cout << "compare with original matrix\n";
    if (!(m == m2))
        throw std::invalid_argument("error in roundtrip");

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

    if (!(m2 == m3))
        throw std::invalid_argument("error in square");

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

template <size_t VecLen, size_t BlkSize = 1>
void testEquivalence(const BinaryMatrix<19937>* commonJump, const BinaryMatrix<19937>* seqJump, size_t commonJumpSize, size_t sequenceJumpSize)
{
    std::cout << "Testing equivalence of generators with SIMD length " << VecLen
        << " and common jump ahead of " << commonJumpSize << " and sequence jump size of " << sequenceJumpSize
        << " and block size " << BlkSize << " ... ";

    const static size_t s_M = VecLen / 32;

    std::vector<uint32_t> dst(nRandomTest + 64 / sizeof(uint32_t));
    uint32_t* aligneddst = (uint32_t*)((intptr_t)dst.data() + (64 - ((intptr_t)dst.data() % 64)));

    MT19937SIMD<VecLen> mt(seedinit, seedlength, commonJump, seqJump);
    for (size_t i = 0; i < nRandomTest / BlkSize; ++i)
        switch (BlkSize) {
            case 1: aligneddst[i] = mt.genrand_uint32(); break;
            case 16: mt.genrand_uint32_blk16(aligneddst + i * BlkSize); break;
            case (624 * (VecLen / 32)):  mt.genrand_uint32_stateBlk(aligneddst + i * (624 * s_M)); break;
            default: throw std::invalid_argument("not implemented");
        };

    for (size_t i = 0; i < nRandomTest; ++i) {
        uint32_t r2 = aligneddst[i];
        size_t seqIndex = i / s_M;
        size_t genIndex = i % s_M;
        size_t benchmarkindex = seqIndex + commonJumpSize + sequenceJumpSize * genIndex;
        if (benchmark[benchmarkindex] != r2) {
            std::cout << "FAILED!\n"
                << "Difference found: out[" << i << "] = " << r2
                << ", benchmark[" << benchmarkindex  << "] = " << benchmark[benchmarkindex] << "\n";
            throw;
        }
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

void recode()
{
    MT19937Matrix f, g;

    for (size_t i = 19931; i <= 19937; ++i) {

        std::ostringstream fIn; fIn  << "C:/workspace/repos/MT19937-SIMD/dat/F" << i << ".b64";
        std::ostringstream fOut; fOut << "C:/workspace/repos/MT19937-SIMD/dat/F" << i << ".bits";

        {
            std::ifstream is(fIn.str());
            f.fromBase64(is);
            f.printSparsity();
        }
        {
            std::ofstream os(fOut.str(), std::ios::binary);
            f.toBin(os);
        }

        g.fromBinFile(fOut.str().c_str());

        if (!(f == g)) {
            std::cout << "binary round trip error" << "\n";
            std::exit(-1);
        }
    }
    exit(0);
}

int main()
{
    //recode();
    try {
#if 1
        encodingTests<19937, 19937>();
        encodingTests<19937, 1007>();
        encodingTests<1007, 19937>();
        encodingTests<1007, 1007>();

        squareTests(std::index_sequence<1, 5, 8, 13, 16, 20, 28, 32, 36, 60, 64, 68, 85, 126, 128, 150>{});
#endif
        generateBenchmark();

        MT19937Matrix jumpMatrix1, jumpMatrix1024("./dat/F00010.bits"), jumpMatrixPeriod("./dat/F19937.bits");

//        initMT19937(jumpMatrix1);
//        getBinaryMatrix("./dat/F00010.bits", jumpMatrix1024);
//        getBinaryMatrix("./dat/F19937.bits", jumpMatrixPeriod);

        testEquivalence<32>(nullptr, nullptr, 0, 0);
        testEquivalence<32>(&jumpMatrix1024, nullptr, 1024, 0);
        // since the period is 2^19937-1, after applying a jump matrix of 2^19937, we restart from the sequence from step 1
        testEquivalence<32>(&jumpMatrixPeriod, nullptr, 1, 0);
        testEquivalence<64>(nullptr, nullptr, 0, 0);
        testEquivalence<128>(nullptr, nullptr, 0, 0);
        testEquivalence<128>(nullptr, &jumpMatrix1, 0, 1);
        testEquivalence<128>(nullptr, &jumpMatrix1024, 0, 1024);
        testEquivalence<128, 16>(nullptr, nullptr, 0, 0);
        testEquivalence<128, 16>(nullptr, &jumpMatrix1, 0, 1);
        testEquivalence<128, 16>(nullptr, &jumpMatrix1024, 0, 1024);
        testEquivalence<128, 624 * 4>(nullptr, nullptr, 0, 0);
        testEquivalence<128, 624 * 4>(nullptr, &jumpMatrix1, 0, 1);
        testEquivalence<128, 624 * 4>(nullptr, &jumpMatrix1024, 0, 1024);

#if SIMD_N_BITS > 128
        testEquivalence<256>(nullptr, nullptr, 0, 0);
        testEquivalence<256>(nullptr, &jumpMatrix1, 0, 1);
        testEquivalence<256>(nullptr, &jumpMatrix1024, 0, 1024);
        testEquivalence<256, 16>(nullptr, nullptr, 0, 0);
        testEquivalence<256, 16>(nullptr, &jumpMatrix1, 0, 1);
        testEquivalence<256, 16>(nullptr, &jumpMatrix1024, 0, 1024);
        testEquivalence<256, 624 * 8>(nullptr, nullptr, 0, 0);
        testEquivalence<256, 624 * 8>(nullptr, &jumpMatrix1, 0, 1);
        testEquivalence<256, 624 * 8>(nullptr, &jumpMatrix1024, 0, 1024);
#endif
#if SIMD_N_BITS > 256
        testEquivalence<512>(nullptr, nullptr, 0, 0);
        testEquivalence<512>(nullptr, &jumpMatrix1, 0, 1);
        testEquivalence<512>(nullptr, &jumpMatrix1024, 0, 1024);
        testEquivalence<512, 16>(nullptr, nullptr, 0, 0);
        testEquivalence<512, 16>(nullptr, &jumpMatrix1, 0, 1);
        testEquivalence<512, 16>(nullptr, &jumpMatrix1024, 0, 1024);
        testEquivalence<512, 624 * 16>(nullptr, nullptr, 0, 0);
        testEquivalence<512, 624 * 16>(nullptr, &jumpMatrix1, 0, 1);
        testEquivalence<512, 624 * 16>(nullptr, &jumpMatrix1024, 0, 1024);
#endif
    }
    catch (const std::exception& e) {
        std::cout << e.what() << "\n";
        return -1;
    }

    return 0;
}
