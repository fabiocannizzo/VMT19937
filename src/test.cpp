#define VRANDGEN_TESTING 1

#include "TestUtils.h"
#include "SIMD.h"

#define SFMT_MEXP 19937
#include "../SFMT-src-1.5.1/SFMT.h"

#include <utility>
#include <numeric>

const uint32_t seedlength = 4;
const uint32_t seedinit[seedlength] = { 0x123, 0x234, 0x345, 0x456 };

const uint64_t nRandomTest = 50ul * 624 * 16;

extern "C" unsigned long genrand_int32();
extern "C" void init_by_array(unsigned long init_key[], int key_length);

enum GenType { VSFMT, VMT, XMT };
const char* genName[] = { "VSFMT", "VMT", "XMT" };

template <GenType G, size_t L, size_t I, QryMode QM>
struct GenTraits;

template <size_t L, size_t I, QryMode QM>
struct GenTraits<VMT, L, I, QM>
{
    typedef VMT19937<L, QM == QM_Block16, BitLenToIsa<I>::isa> gen_t;
};

template <size_t L, size_t I, QryMode QM>
struct GenTraits<XMT, L, I, QM>
{
    static_assert(L == I);
    typedef XMT19937<L, QM == QM_Block16, BitLenToIsa<I>::isa> gen_t;
};

template <size_t L, size_t I, QryMode QM>
struct GenTraits<VSFMT, L, I, QM>
{
    typedef VSFMT19937<L, QM == QM_Block16, BitLenToIsa<I>::isa> gen_t;
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

template <GenType G, size_t L, size_t I, QryMode QM, typename M>
void testEquivalence(size_t nCommonJumpRepeat, const JumpMatrix<M>& commonJump, const JumpMatrix<M>& seqJump)
{
    using Gen = typename GenTraits<G, L, I, QM>::gen_t;

    const size_t commonJumpSize = commonJump.p ? commonJump.jumpSize : 0;
    const size_t sequenceJumpSize = seqJump.p ? seqJump.jumpSize : 0;

    MYASSERT(((commonJump.p != nullptr) == (nCommonJumpRepeat > 0)), "commnJump matrix should be provided only if nCommonJumpRepeat>0");

    constexpr size_t VecLen = Gen::s_regLenBits;
    constexpr QryMode QryMode = QM;
    size_t blkSize;
    switch (QryMode) {
        case QM_Any: blkSize = 0; break;
        case QM_Scalar: blkSize = 1; break;
        case QM_Block16: blkSize = 16; break;
        default: THROW("how did we get here?");
    }
    constexpr size_t s_nStates = Gen::s_nStates;
    constexpr size_t s_n32InOneWord = Gen::s_n32InOneWord;

    std::cout << genName[G] << "< " << std::setw(3) << VecLen << ", "
        << std::setw(7) << queryModeName(QryMode) << ", " << std::setw(3) << Gen::s_regLenBitsHw << ">"
        << ", common jump of " << std::setw(4) << commonJumpSize << " repeated " << nCommonJumpRepeat << " times, sequence jump of " << std::setw(4) << sequenceJumpSize
        << ", block size " << std::setw(4);
    if (blkSize > 0)
        std::cout << blkSize;
    else
        std::cout << "rand";
    std::cout << " ... ";

    std::vector<uint32_t> aligneddst(nRandomTest);

    const M* jumpMat = nullptr;
    if constexpr (s_nStates > 1)
        jumpMat = seqJump.p.get();
    else
        MYASSERT(!seqJump.p, "sequential jump provided for single-state generator");

    Gen mt(seedinit, seedlength, nCommonJumpRepeat, commonJump.p.get(), jumpMat);

    uint32_t* dst = aligneddst.data();
    if constexpr (QryMode != QM_Any) {
        for (size_t i = 0; i < nRandomTest / blkSize; ++i)
            if constexpr (QryMode == QM_Scalar)
                *dst++ = mt.genrand_uint32();
            else if constexpr (QryMode == QM_Block16) {
                mt.genrand_uint32_blk16(dst);
                dst += blkSize;
            }
            else
                NOT_IMPLEMENTED;
    }
    else { // QryMode == QM_Any
        size_t n = nRandomTest;
        uint32_t* dst = aligneddst.data();
        const size_t sz[] = { 1, 2, 3, 4, 7, 8, 9, 15, 16, 17, 31, 32, 33, 127, 128, 129, 623, 624, 625, 800
                            , 624 * 2 - 1, 624 * 2, 624 * 2 + 1
                            , 624 * 4 - 1, 624 * 4, 624 * 4 + 1
                            , 624 * 8 - 1, 624 * 8, 624 * 8 + 1
                            };
        while (n) {
            size_t szi = std::rand() % (sizeof(sz) / sizeof(*sz));
            size_t m = std::min<size_t>(n, sz[szi]);
            mt.genrand_uint32_anySize(dst, m);
            dst += m;
            n -= m;
        }
    }

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


template <GenType G, size_t L, size_t I, QryMode QM, typename M>
void equivalenceTests3(const JumpMatrix<M>& jumpSmall, const JumpMatrix<M>& jumpBig)
{
    JumpMatrix<M> noJump{};

    using Gen = typename GenTraits<G, L, I, QM>::gen_t;

    if constexpr (I <= L && I <= SIMD_N_BITS) {
        if constexpr (QM != QM_Any) {
            testEquivalence<G, L, I, QM>(0, noJump, noJump);
            testEquivalence<G, L, I, QM>(1, jumpSmall, noJump);
            testEquivalence<G, L, I, QM>(2, jumpSmall, noJump);
            testEquivalence<G, L, I, QM>(1, jumpBig, noJump);
            if constexpr (L > 32 && Gen::s_nStates > 1) {
                testEquivalence<G, L, I, QM>(1, jumpBig, jumpBig);
                testEquivalence<G, L, I, QM>(2, jumpBig, jumpBig);
                testEquivalence<G, L, I, QM>(0, noJump, jumpBig);
            }
        }
        else {
            // we repeat this test multiple times, as there are random number involved
            for (size_t i = 0; i < 10; ++i)
                testEquivalence<G, L, I, QM>(0, noJump, noJump);
        }
    }
}

template <GenType G, size_t L, size_t I, QryMode...QMs, typename M>
void equivalenceTests2(const JumpMatrix<M>& jumpSmall, const JumpMatrix<M>& jumpBig)
{
    (equivalenceTests3<G, L, I, QMs>(jumpSmall, jumpBig), ...);
}

template <GenType G, size_t L, size_t...Is, typename M>
void equivalenceTests1(const JumpMatrix<M>& jumpSmall, const JumpMatrix<M>& jumpBig)
{
    (equivalenceTests2<G, L, Is, QM_Scalar, QM_Block16, QM_Any>(jumpSmall, jumpBig), ...);
}

template <GenType G, size_t...Ls, typename M>
void equivalenceTests0(const JumpMatrix<M>& jumpSmall, const JumpMatrix<M>& jumpBig)
{
    if constexpr (G == XMT)
        (equivalenceTests1<G, Ls, Ls>(jumpSmall, jumpBig), ...);
    else
        (equivalenceTests1<G, Ls, 32, 128, 256, 512>(jumpSmall, jumpBig), ...);
}

void test_XVMT19937()
{
    generateBenchmark_MT19937();

    typedef MT19937Matrix matrix_t;
    typedef JumpMatrix<matrix_t> pmatrix_t;

    pmatrix_t noJump;
    pmatrix_t jumpMatrix1(new matrix_t, 1);                                          // jump ahead 1 element
    pmatrix_t jumpMatrix512(new matrix_t(std::string("./dat/mt/F00009.bits")), 512);    // jump ahead 2^9 (512) elements
    pmatrix_t jumpMatrixPeriod(new matrix_t(std::string("./dat/mt/F19937.bits")), 1);   // jump ahead 2^19937 elements

    // test VMT generator
    startTest(genName[VMT]);
    equivalenceTests0<VMT, 32, 128, 256, 512>(jumpMatrix1, jumpMatrix512);
    // since the period is 2^19937-1, after applying a jump matrix of 2^19937, we restart the sequence from step 1
    std::cout << "VMT19937: a jump of size 2^19937 is equivalent to a jump of size 1\n";
    testEquivalence<VMT, 32, 32, QM_Scalar>(1, jumpMatrixPeriod, noJump);

    // test XMT generator
    startTest(genName[XMT]);
    equivalenceTests0<XMT, 32, 128, 256, 512>(jumpMatrix1, jumpMatrix512);
    // since the period is 2^19937-1, after applying a jump matrix of 2^19937, we restart the sequence from step 1
    std::cout << "XMT19937: a jump of size 2^19937 is equivalent to a jump of size 1\n";
    testEquivalence<XMT, 32, 32, QM_Scalar>(1, jumpMatrixPeriod, noJump);
    testEquivalence<XMT, 128, 128, QM_Scalar>(1, jumpMatrixPeriod, noJump);
#if SIMD_N_BITS>=256
    testEquivalence<XMT, 256, 256, QM_Scalar>(1, jumpMatrixPeriod, noJump);
#endif
#if SIMD_N_BITS>=512
    testEquivalence<XMT, 512, 512, QM_Scalar>(1, jumpMatrixPeriod, noJump);
#endif
}

void test_VSFMT19937()
{
    startTest(genName[VSFMT]);

    generateBenchmark_SFMT19937();

    typedef SFMT19937Matrix matrix_t;
    typedef JumpMatrix<matrix_t> pmatrix_t;

    pmatrix_t noJump;
    pmatrix_t jumpMatrix4(new matrix_t, 4);                                              // jump ahead 1 element
    pmatrix_t jumpMatrix512(new matrix_t(std::string("./dat/sfmt/F00009.bits")), 512);   // jump ahead 2^9 (1024) elements

    equivalenceTests0<VSFMT, 128, 256, 512>(jumpMatrix4, jumpMatrix512);
}

template <typename T>
void printReg(std::string&& name, T v)
{
    constexpr unsigned n = sizeof(T);
    alignas(sizeof(T)) unsigned char bytes[n];
    std::copy_n((unsigned char*)&v, n, bytes);
    std::cout << std::setw(3) << name << ": ";
    for (unsigned i = 0; i < n; ++i)
        std::cout << std::setw(3) << (int)bytes[i];
    std::cout << '\n';
}

#define MYASSERT_XV_EQ(a, b, msg) MYASSERT(a.eq(b), msg)

template <size_t n32, typename T>
void testAlignR32(const unsigned char *data, T a, T b)
{
    auto c = T::template alignr32<n32>(a, b);
    const unsigned char* got = (const unsigned char*)  &c;
    printReg(std::to_string(n32), c);
    for (size_t i = 0; i < sizeof(T); ++i) {
        size_t srcIndex = i + n32 * 4;
        unsigned char expected = data[srcIndex];
        MYASSERT(got[i] == expected, "error in alignR32<" << n32 << ">: got[" << i << "]=" << (int)got[i] << ", expected=" << (int)expected);
    }
}

template <size_t...n32s>
void testSimdAlignR32(std::index_sequence<n32s...>&&)
{
    constexpr size_t n32 = (sizeof...(n32s) - 1);
    constexpr size_t nBits = n32 * 32;
    using T = Details::SimdRegister<nBits, BitLenToIsa<nBits>::isa>;
    std::cout << "\nTest SimdRegister<" << nBits << ", " << (int)BitLenToIsa<nBits>::isa << ">::alignr32\n";
    alignas(64) unsigned char data[128];
    std::iota(data, data + 128, 0);

    T v0((const T*)data);
    T v1(((const T*)data)+1);

    printReg("v0", v0);
    printReg("v1", v1);

    (testAlignR32<n32s>(data, v0, v1), ...);
}

void test_SIMD_special_methods()
{
    std::cout << "\n--- SIMD special methods tests ---\n";

    using XV = Details::SimdRegister<128, BitLenToIsa<128>::isa>;

    // Test alignr32
    {
        XV a(1, 2, 3, 4);
        XV b(5, 6, 7, 8);

        XV r0 = XV::alignr32<0>(a, b);
        MYASSERT(r0.eq(a), "alignr32<0> failed");

        XV r1 = XV::alignr32<1>(a, b);
        MYASSERT(r1.eq(XV(2, 3, 4, 5)), "alignr32<1> failed");

        XV r2 = XV::alignr32<2>(a, b);
        MYASSERT(r2.eq(XV(3, 4, 5, 6)), "alignr32<2> failed");

        XV r3 = XV::alignr32<3>(a, b);
        MYASSERT(r3.eq(XV(4, 5, 6, 7)), "alignr32<3> failed");

        XV r4 = XV::alignr32<4>(a, b);
        MYASSERT(r4.eq(b), "alignr32<4> failed");
    }

    // Test shl128 / shr128
    {
        XV a(0x01020304u, 0x05060708u, 0x090A0B0Cu, 0x0D0E0F10u);

        XV l4 = XV::shl128<4>(a);
        MYASSERT(l4.eq(XV(0, 0x01020304u, 0x05060708u, 0x090A0B0Cu)), "shl128<4> failed");

        XV r4 = XV::shr128<4>(a);
        MYASSERT(r4.eq(XV(0x05060708u, 0x090A0B0Cu, 0x0D0E0F10u, 0)), "shr128<4> failed");

        XV l8 = XV::shl128<8>(a);
        MYASSERT(l8.eq(XV(0, 0, 0x01020304u, 0x05060708u)), "shl128<8> failed");

        XV r8 = XV::shr128<8>(a);
        MYASSERT(r8.eq(XV(0x090A0B0Cu, 0x0D0E0F10u, 0, 0)), "shr128<8> failed");
    }

    // Test ifOddCst32ElseZero
    {
        XV a(1, 2, 3, 4); // odd, even, odd, even
        XV cst(0xAAAAAAAAu);
        XV res = a.ifOddCst32ElseZero(cst);
        MYASSERT(res.eq(XV(0xAAAAAAAAu, 0, 0xAAAAAAAAu, 0)), "ifOddCst32ElseZero failed");
    }


    // Test parity
    {
        XV a(1, 0, 0, 0);
        MYASSERT(a.parity() == 1, "parity failed (1)");

        XV b(1, 1, 0, 0);
        MYASSERT(b.parity() == 0, "parity failed (0)");

        XV c(0x12345678u, 0x87654321u, 0x11223344u, 0x44332211u);
        uint32_t p = popcnt(0x12345678u) ^ popcnt(0x87654321u) ^ popcnt(0x11223344u) ^ popcnt(0x44332211u);
        MYASSERT(c.parity() == (uint8_t)(p & 1), "parity failed (complex)");
    }

    std::cout << "SIMD special methods tests passed!\n";
}

int main()
{
    try {
        test_SIMD_special_methods();
#if 0
        testSimdAlignR32(std::make_index_sequence<128 / 32 + 1>{});
#if SIMD_N_BITS>=256
        testSimdAlignR32(std::make_index_sequence<256 / 32 + 1>{});
#endif
#if SIMD_N_BITS>=512
        testSimdAlignR32(std::make_index_sequence<512 / 32 + 1>{});
#endif
        testEncoding();
        testSquareMatrix();
#endif
        test_XVMT19937();
        test_VSFMT19937();
    }
    catch (const std::exception& e) {
        std::cout << e.what() << "\n";
        return -1;
    }

    return 0;
}
