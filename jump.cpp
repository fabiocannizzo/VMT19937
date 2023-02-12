#include "bit_matrix.h"

#include <thread>
#include <memory>
#include <sstream>


void wait()
{
    std::cout << "press a key...";
    char c;
    std::cin >> c;
}



struct F : BinaryMatrix<19937>
{
    typedef BinaryMatrix<19937> base_t;

    typedef uint32_t word_t;

    static const uint32_t s_matA = 0x9908B0DF;
    static const uint32_t s_M = 397;
    static const size_t s_nBits = 19937;
    static const size_t s_nWordBits = sizeof(word_t) * 8;

    void init()
    {
        // from row 0 to to row nBits - 32, state bits are just shifted left by 32 bits
        for (uint32_t r = 0; r < s_nBits - s_nWordBits; ++r)
            setBit(r, r + s_nWordBits);

        // the new state element is composed of element which was in position M
        // ...
        for (uint32_t i = 0; i < s_nWordBits; i++)
            setBit(s_nBits - s_nWordBits + i, 1 + (s_M - 1) * s_nWordBits + i);
        for (uint32_t i = 0; i < s_nWordBits; ++i)
            if (s_matA & (uint32_t(1) << i))
                setBit(s_nBits - s_nWordBits + i, 1);
        setBit(s_nBits - 2, 0);
        for (uint32_t i = 0; i < s_nWordBits - 2; ++i)
            setBit(s_nBits - s_nWordBits + i, 2 + i);
    }

    static const size_t s_nBitColsPerBlk = 8;
    const size_t s_nColBlks = s_nBits / s_nBitColsPerBlk;
    typedef BinaryMatrix<s_nBitColsPerBlk, base_t::s_nBitRows> buffer_t;

    template <size_t NRows, size_t nColumns, size_t...Is>
    static void transpose16x8(buffer_t& buffer, const base_t& src, size_t rowBitIndex, size_t colBitIndex,  std::index_sequence<Is...>&&)
    {
        const size_t simdBits = 128;
        const size_t nRowsPerBlk = simdBits / s_nBitColsPerBlk;
        //const size_t nRowBlks = s_nBits / nRowsPerBlk;

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

    //static bool multiplySlow(const uint8_t* pr, const uint8_t* pc)
    //{
    //    uint64_t cnt = 0;
    //    for (size_t j = 0; j < base_t::s_nBitCols; ++j) {
    //        auto [ind, off] = bitPos(j);
    //        cnt += (bool)((pr[ind] & bitmask(off)) & (pc[ind] & bitmask(off)));
    //    }
    //    bool result = cnt & 1;
    //    return result;
    //}

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
        //std::cout << colBuffer.countAllBits() << ",";

        const uint8_t* pcs[nColumns];
        for (size_t c = 0; c < nColumns; ++c)
            pcs[c] = colBuffer.rowBegin(c);

        for (size_t r = 0; r < s_nBits; ++r) {
            auto pr = src.rowBegin(r);
            uint8_t active = BinaryVectorMultiplier<SIMD_N_BITS>::multiply8<nColumns, base_t::s_nBitCols, base_t::s_nBitColsPadded>(pr, pcs);
            rowBegin(r)[colBitIndex / 8] = active;
        }
/*

#if 1
                if (target) {
                    bool check = target->getBit(r, c + colBitIndex);
                    if (check != active) {
                        //bool check2 = multiplySlow(pr, pc);
                        bool check3 = src.multiplyRowByCol(r, c + colBitIndex);
                        std::cout << "error at (" << r << "," << c << "): matlab " << check
                            << ", got " << active
                            //<< ", slow " << check2
                            << ", raw " << check3
                            << "\n";
                        wait();
                    }
                }
#endif
                if (active)
                    setBit(r, colBitIndex + c);
            }
        }
        */
    }

    void squareBlock(const base_t& f, buffer_t& colBuffer, std::atomic<size_t>& nextJobIndex, const base_t* target)
    {
        while (true) {
            const size_t curJobIndex = nextJobIndex++;
            const size_t colBitIndex = curJobIndex * s_nBitColsPerBlk;
            if (curJobIndex < s_nColBlks) {
                multiplyBlock<s_nBitColsPerBlk>(colBuffer, f, colBitIndex, target);
            }
            else if (curJobIndex == s_nColBlks) { // in the last block we have less than s_nBitColsPerBlk columns to process
                if (s_nBitCols % s_nBitColsPerBlk)
                    multiplyBlock<s_nBitCols% s_nBitColsPerBlk>(colBuffer, f, colBitIndex, target);
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
            threads[i].reset(new std::thread(&F::squareBlock, this, std::cref(src), std::ref(buffers[i]), std::ref(nextJobIndex), target));

        // main thread participate to the loop
        squareBlock(src, buffers.back(), nextJobIndex, target);

        for (auto& th : threads)
            th->join();
    }

};

void printSparsity(const F& src)
{
    auto [n, frac] = src.sparseIndex();
    std::cout << "nnz: " << n << ", sparse index: " << std::fixed << std::setprecision(4) << 100.0 * frac << "%\n";
}

void square(const F& src, F& dst, std::vector<F::buffer_t>& buffers, const F* target)
{
    auto start = std::chrono::system_clock::now();
    dst.square(src, buffers, target);
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "done in: " << std::fixed << std::setprecision(2) << elapsed_seconds.count() << "s" << std::endl;
}

int main(int narg, const char *args[])
{
    F f[2];
    f[0].init();
    f[0].printBits(F::s_nBits - 32, 0, 32, 32);
    std::cout << "f0: ";  printSparsity(f[0]);

    F m[3];
    for (size_t i = 0; i < 3; ++i) {
        char n[] = { 'm', (char) (i + '0'), 0 };
        std::ostringstream os;
        os << "r" << i << ".txt";
        m[i].fromMatlabSparseFile(os.str().c_str());
        std::cout << n << ": ";
        printSparsity(m[i]);
    }

    if (!(f[0] == m[0]))
        std::cout << "base matrix is different\n";

    //f[0].base64Out("./", "F00000");

    const size_t nThreads = 4;

    std::vector<F::buffer_t> buffers(nThreads);

    //wait();

    //std::cout << "computing f1\n";
    //square(f[0], f[1], buffers, &m[1]);
    //std::cout << "f1: ";  printSparsity(f[1]);
    //if (!(f[1] == m[1])) {
    //    std::cout << "matlab matrix is different\n";
    //    wait();
    //}

    //wait();

    //while (true) {
    //    f[0].resetZero();
    //    std::cout << "computing f2\n";
    //    square(f[1], f[0], buffers, &m[2]);
    //    if (!(f[0] == m[2])) {
    //        std::cout << "matlab matrix is different\n";
    //        wait();
    //    }
    //    std::cout << "f2: ";  printSparsity(f[0]);
    //}

    for (size_t i = 0; i < 19936; ++i) {
        std::cout << "computing F^(2^" << i + 1 << ")\n";
        size_t in = i % 2;
        size_t out = (i + 1) % 2;
        f[out].resetZero();
//        if (i == 0)
            square(f[in], f[out], buffers, nullptr);
//        else
//            square(f[in], f[out], buffers, &target);
        //if (i==0)
        //    f[1].matlabOut("xR2");
        if (i < 2)
            if (!(f[out] == m[i+1]))
                std::cout << "base matrix is different\n";
        printSparsity(f[out]);

        std::ofstream os("temp");
        f[out].toBase64(os);
        os.close();
        std::ifstream is("temp");
        f[in].resetZero();
        f[in].fromBase64(is);
        is.close();
        if (!(f[in] == f[out]))
            std::cout << "read save problem\n";
/*
        if ((i % 200) == 0 || i > 19930) {
            std::ostringstream os;
            os << "F" << std::setw(5) << std::setfill('0') << i+1;
            f[out].base64Out("", os.str().c_str());
        }
*/

    }

    return 0;
}
