#include "bit_matrix.h"

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
        throw std::exception("error in roundtrip");

    std::cout << "completed\n";
}

template <size_t N>
void testSquare(const BinaryMatrix<N, N>& m)
{
    BinaryMatrix<N, N> m2, m3;

    // slow bit by bit multiplication
    for (size_t r = 0; r < N; ++r) {
        for (size_t c = 0; c < N; ++c) {
            size_t s = 0;
            for (size_t k = 0; k < N; ++c) {
                s += m.getBit(r, k) & m.getBit(k, c);
            }
            if (s & 1)
                m2.setBit(r, c);
        }
    }


}

template <size_t nRows, size_t nCols>
void runTests()
{
    BinaryMatrix<nRows, nCols> m;
    m.initRand();
    std::cout << "\ngenerated random matrix with size (" << m.s_nBitCols << "x" << m.s_nBitCols << ") with " << m.nnz() << " non zero elements\n";
    m.printBits(0, 0, 10, 32);

    testEncoder(m, Base64);
    testEncoder(m, Hex);

    if constexpr (nRows == nCols)
        testSquare(m);
}

int main()
{
    try {
        runTests<19937, 19937>();
        runTests<19937, 1007>();
        runTests<1007, 19937>();

        //testMatrixSquare;
        //testMatrixByVectorMult;

        //testMatlabEncoder;
        //testF0Matrix;
    }
    catch (const std::exception& e) {
        std::cout << e.what() << "\n";
        return -1;
    }

    return 0;
}