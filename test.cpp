#include "jump_matrix.h"

#include "MT19937-SIMD.h"

const uint32_t seedlength = 4;
const uint32_t seedinit[seedlength] = { 0x123, 0x234, 0x345, 0x456 };

const uint64_t nRandomTest = 500000;

extern "C" unsigned long genrand_int32();
extern "C" void init_by_array(unsigned long init_key[], int key_length);


std::vector<uint32_t> benchmark(nRandomTest);

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

template <size_t N>
void testSquare(const BinarySquareMatrix<N>& m)
{
    BinarySquareMatrix<N> m2, m3;

    // slow bit by bit multiplication
    //std::cout << "compute matrix multiplication the classical way\n";
    for (size_t r = 0; r < N; ++r) {
        //std::cout << r << "\n";
        for (size_t c = 0; c < N; ++c) {
            size_t s = 0;
            for (size_t k = 0; k < N; ++k) {
                s ^= m.getBit(r, k) && m.getBit(k, c);
            }
            if (s)
                m2.setBit(r, c);
        }
    }

    const size_t nThreads = 4;
    //std::cout << "compute matrix multiplication vectorially\n";
    std::vector<typename BinarySquareMatrix<N>::buffer_t> buffers(nThreads);
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

template <size_t N>
void squareTest()
{
    std::cout << "testing multiplication with matrices of size: " << N << "\n";
    BinarySquareMatrix<N> m;
    for (size_t i = 0; i < 10; ++i) {
        m.resetZero();
        m.initRand();
        testSquare(m);
    }
}

template <size_t...N>
void squareTests(std::index_sequence<N...>&&)
{
    (squareTest<N>(), ...);
}

template <size_t VecLen>
void testEquivalence(const BinaryMatrix<19937>* jumpMatrix = NULL, size_t nJumps = 0)
{
    std::cout << "Testing equivalence of generators with SIMD length " << VecLen << " and jump ahead of " << nJumps << " ... ";

    const size_t M = VecLen / 32;

    MT19937SIMD<VecLen> mt(seedinit, seedlength, jumpMatrix);

    for (size_t i = 0; i < nRandomTest; ++i) {
        uint32_t r2 = mt.genrand_uint32();
        size_t seqIndex = i / M;
        size_t genIndex = i % M;
        if (benchmark[seqIndex + nJumps * genIndex] != r2) {
            std::cout << "FAILED!\n"
                << "Difference found at index " << i << ": expected " << benchmark[i] << ", but got " << r2 << "\n";
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
    for (size_t i = 0; i < nRandomTest; ++i)
        benchmark[i] = (uint32_t)genrand_int32();
    std::cout << "done!\n";
    printSome(benchmark);
}

int main()
{
    try {
#if 0
        encodingTests<19937, 19937>();
        encodingTests<19937, 1007>();
        encodingTests<1007, 19937>();
        encodingTests<1007, 1007>();

        squareTests(std::index_sequence<1, 5, 8, 13, 16, 20, 28, 32, 36, 60, 64, 68, 85, 126, 128, 150>{});
#endif
        generateBenchmark();

        MT19937Matrix jumpMatrix;
        initMT19937(jumpMatrix);

        testEquivalence<32>();
        testEquivalence<64>();
        testEquivalence<128>();
        testEquivalence<128>(&jumpMatrix, 1);
#if SIMD_N_BITS > 128
        testEquivalence<256>();
        testEquivalence<256>(&jumpMatrix, 1);
#endif
#if SIMD_N_BITS > 256
        testEquivalence<512>();
        testEquivalence<512>(&jumpMatrix, 1);
#endif



    }
    catch (const std::exception& e) {
        std::cout << e.what() << "\n";
        return -1;
    }

    return 0;
}
