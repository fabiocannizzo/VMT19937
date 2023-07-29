#include "bit_matrix.h"

#include <thread>
#include <atomic>
#include <memory>
#include <sstream>

template <size_t N>
struct BinarySquareMatrix : BinaryMatrix<N, N>
{
    typedef BinaryMatrix<N,N> base_t;

    typedef uint32_t word_t;

    static const size_t s_nBits = N;
    static const size_t s_nWordBits = sizeof(word_t) * 8;

    static const size_t s_nBitCols = base_t::s_nBitCols;
    static const size_t s_nBitColsPerBlk = 8;
    const size_t s_nColBlks = s_nBits / s_nBitColsPerBlk;
    typedef BinaryMatrix<s_nBitColsPerBlk, base_t::s_nBitRows> buffer_t;

    template <size_t NRows, size_t nColumns, size_t...Is>
    static void transpose16x8(buffer_t& buffer, const base_t& src, size_t rowBitIndex, size_t colBitIndex,  std::index_sequence<Is...>&&)
    {
        const size_t simdBits = 128;
        const size_t nRowsPerBlk = simdBits / s_nBitColsPerBlk;

        static_assert(s_nBitColsPerBlk == 8);

        const size_t n = sizeof...(Is) - 1;
        auto cs = _mm_set_epi8((n - Is < NRows ? src.template getWordU<uint8_t>(rowBitIndex + (n - Is), colBitIndex) : 0)...);
        
        for (size_t cb = 0; cb < nColumns; ++cb) {
            auto mask = _mm_movemask_epi8(_mm_slli_epi64(cs, (int)(s_nBitColsPerBlk - 1 - cb)));
            uint16_t* pc = (uint16_t*)buffer.rowBegin(cb);
            pc[rowBitIndex / (8 * sizeof(uint16_t))] = (uint16_t)mask;
        }
    }

    // transpose columns from colBitIndex to colBitIndex+nBitColsPerBlk-1 and copies then to colBuffer
    // so that then is it more efficient to compute the square of the matrix, where we have to multiply rows by columns
    template <size_t nColumns>
    static void transposeColumns(buffer_t& colBuffer, const base_t& src, size_t colBitIndex)
    {
        const size_t simdBits = 128;
        const size_t nRowsPerBlk = simdBits / s_nBitColsPerBlk;
        const size_t nRowBlks = s_nBits / nRowsPerBlk;
        const size_t nResidualRows = base_t::s_nBitRows - nRowBlks * nRowsPerBlk;

        colBuffer.resetZero();

        for (size_t i = 0; i < nRowBlks; ++i)
            transpose16x8<nRowsPerBlk, nColumns>(colBuffer, src, i * nRowsPerBlk, colBitIndex, std::make_index_sequence<nRowsPerBlk>{});

        if (nResidualRows)
            transpose16x8<nResidualRows, nColumns>(colBuffer, src, nRowBlks * nRowsPerBlk, colBitIndex, std::make_index_sequence<nRowsPerBlk>{});
    }

    template <size_t nColumns, typename DST>
    void multiplyBlock(DST& colBuffer, const base_t& src, size_t colBitIndex, const base_t *target)
    {
        // transpose columns from colBitIndex to colBitIndex+nBitColsPerBlk-1 and copies then to colBuffer
        // so that then is it more efficient to compute the square of the matrix, where we have to multiply rows by columns
        transposeColumns<nColumns>(colBuffer, src, colBitIndex);

#if 0 || defined(_DEBUG)
        static bool firstError = true;
        for (size_t r = 0; r < base_t::s_nBitRows; ++r) {
            if (!firstError)
                break;
            for (size_t c = 0; c < nColumns; ++c)
                if (src.getBit(r, c + colBitIndex) != colBuffer.getBit(c, r)) {
                    std::cout << "error at coord " << r << "," << c << "\n";
                    src.printBits(r / 16 * 16, (c + colBitIndex) / 8 * 8, 16, 8);
                    colBuffer.printBits(0, r / 16 * 16, 8, 16);
                    firstError = true;
                    break;
                }
        }
#endif

        const uint8_t* pcs[nColumns];
        for (size_t c = 0; c < nColumns; ++c)
            pcs[c] = colBuffer.rowBegin(c);

        for (size_t r = 0; r < s_nBits; ++r) {
            auto pr = src.rowBegin(r);
            uint8_t active = BinaryVectorMultiplier<SIMD_N_BITS>::multiply8<nColumns, base_t::s_nBitCols, base_t::s_nBitColsPadded>(pr, pcs);
            this->rowBegin(r)[colBitIndex / 8] = active;
        }
    }

    void squareBlock(const base_t& src, buffer_t& colBuffer, std::atomic<size_t>& nextJobIndex, const base_t* target)
    {
        while (true) {
            const size_t curJobIndex = nextJobIndex++;
            const size_t colBitIndex = curJobIndex * s_nBitColsPerBlk;
            if (curJobIndex < s_nColBlks) {
                multiplyBlock<s_nBitColsPerBlk>(colBuffer, src, colBitIndex, target);
            }
            else if (curJobIndex == s_nColBlks) { // in the last block we have less than s_nBitColsPerBlk columns to process
                constexpr size_t nColumns = s_nBitCols % s_nBitColsPerBlk;
                if constexpr (nColumns)
                    multiplyBlock<nColumns>(colBuffer, src, colBitIndex, target);
            }
            else
                break;
        }
    }

    void square(const base_t& src, std::vector<buffer_t>& buffers, const base_t* target)
    {
        std::vector< std::shared_ptr<std::thread>> threads(buffers.size() - 1);

        std::atomic<size_t> nextJobIndex(0);

        // launch participating threads
        for (size_t i = 0; i < threads.size(); ++i)
            threads[i].reset(new std::thread(&BinarySquareMatrix::squareBlock, this, std::cref(src), std::ref(buffers[i]), std::ref(nextJobIndex), target));

        // main thread participate to the loop
        squareBlock(src, buffers.back(), nextJobIndex, target);

        for (auto& th : threads)
            th->join();
    }

    void printSparsity()
    {
        auto [n, frac] = this->sparseIndex();
        std::cout << "nnz: " << n << ", sparse index: " << std::fixed << std::setprecision(4) << 100.0 * frac << "%\n";
    }
};

struct MT19937Matrix : BinarySquareMatrix<19937>
{
    typedef BinarySquareMatrix<19937> base_t;

    // Initialize the matrix as per MT19937 32 bit generator transition matrix
    // This is equivalent to a jump ahead of 1 random number
    MT19937Matrix()
    {
        init1();
    }

    MT19937Matrix(const char* filename)
    {
        fromBinFile(filename);
    }

    void init1()
    {
        static const size_t s_nBits = base_t::s_nBitRows;
        static const size_t s_nWordBits = base_t::s_nWordBits;
        static const uint32_t s_matA = 0x9908B0DF;
        static const uint32_t s_M = 397;

        // from row 0 to to row nBits - 32, state bits are just shifted left by 32 bits
        for (uint32_t r = 0; r < s_nBits - s_nWordBits; ++r)
            setBit(r, r + s_nWordBits);

        // the new state element is composed of element which was in position M
        for (uint32_t i = 0; i < s_nWordBits; i++)
            setBit(s_nBits - s_nWordBits + i, 1 + (s_M - 1) * s_nWordBits + i);
        for (uint32_t i = 0; i < s_nWordBits; ++i)
            if (s_matA & (uint32_t(1) << i))
                setBit(s_nBits - s_nWordBits + i, 1);
        setBit(s_nBits - 2, 0);
        for (uint32_t i = 0; i < s_nWordBits - 2; ++i)
            setBit(s_nBits - s_nWordBits + i, 2 + i);
    }

    // initialize from a binary file saved with the toBin method
    void fromBinFile(const char* filename)
    {
        std::ifstream is(filename, std::ios::binary);
        MYASSERT(is.is_open(), "error opening binary file: " << filename);
        fromBin(is);
#ifdef TESTING
        std::cout << "loaded matrix from file: " << filename << "\n";
        printSparsity();
#endif
    }
};
