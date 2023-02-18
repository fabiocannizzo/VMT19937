#pragma once

#include "SIMD.h"
#include "base64.h"

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstring>

inline constexpr uint8_t bitmask(size_t b)
{
    return (uint8_t)(uint8_t(1) << b);
}

// simd implementation
template <size_t SimdBits>
struct BinaryVectorMultiplier
{
    // multyply the row pointed by pr by all the nRows pointed by pc
    // all rows must have size nBitCols and be padded to the next multiple of the simd vector length
    template <size_t nRows, size_t nBitCols, size_t nBitColsPadded>
    static uint8_t multiply8(const uint8_t* _pr, const uint8_t** _pc)
    {
        typedef SimdRegister<SimdBits> simd_t;

        const size_t nSimdBytes = sizeof(simd_t);
        const size_t nSimdBits = 8 * nSimdBytes;
        const size_t nSimdStep = nBitCols / nSimdBits + ((nBitCols % nSimdBits) != 0);
        static_assert(nSimdStep * nSimdBits <= nBitColsPadded);

        typename simd_t::Konst accumulator[nRows];
        for (size_t i = 0; i < nSimdStep; ++i) {
            simd_t a(_pr + i * nSimdBytes);
            for (size_t c = 0; c < nRows; ++c) {
                simd_t b(_pc[c] + i * nSimdBytes);
                simd_t p = a & b;
                accumulator[c] ^= p;
            }
        }

        uint8_t result = 0;
        static_assert(nRows <= 8 * sizeof(result));

        for (size_t c = 0; c < nRows; ++c)
            result |= accumulator[c].parity() << c;

        return result;
    }
};


template <size_t _nBitRows, size_t _nBitCols = _nBitRows>
class BinaryMatrix
{
    static const size_t s_nAlignBits = 512;  // suitable for AVX512
    static const size_t s_nAlignBytes = s_nAlignBits / 8;  // suitable for AVX512
    static_assert(s_nAlignBits == 8 * s_nAlignBytes);

    std::vector<uint8_t> m_storage;
    uint8_t* m_data;

public:
    static const size_t s_nBitRows = _nBitRows;
    static const size_t s_nBitCols = _nBitCols;
    static const size_t s_nPaddingBits = (s_nAlignBits - (s_nBitCols % s_nAlignBits)) % s_nAlignBits;
    static const size_t s_nBitColsPadded = s_nBitCols + s_nPaddingBits;
    static const size_t s_nBytesPerRow = s_nBitCols / 8+ (s_nBitCols % 8 != 0);
    static const size_t s_nBytesPerPaddedRow = s_nBitColsPadded / 8;  // inclusive of padding
    static const size_t s_nUsedBytes = s_nBitRows * s_nBytesPerPaddedRow;

private:
    uint8_t& getByte(size_t rowIndex, size_t colByteIndex) { return m_data[rowIndex * s_nBytesPerPaddedRow + colByteIndex]; }
    const uint8_t& getByte(size_t rowIndex, size_t colByteIndex) const { return m_data[rowIndex * s_nBytesPerPaddedRow + colByteIndex]; }

    // compute the index in the uint8_t matrix and the offset within the uint8_t element
    static auto bitPos(size_t bitColIndex)
    {
        size_t byteIndex = bitColIndex / 8;
        uint8_t bitOffset = bitColIndex % 8;
        return std::make_pair(byteIndex, bitOffset);
    }

public:


public:

    BinaryMatrix()
        : m_storage(s_nUsedBytes + (s_nAlignBytes - 1), 0)
    {
        // align pointer
        auto addr = reinterpret_cast<std::uintptr_t>(m_storage.data());
        auto moveFwdBy = (s_nAlignBytes - (addr % s_nAlignBytes)) % s_nAlignBytes;
        //if (moveFwdBy >= s_nAlignBytes) {
        //    std::cout << "something wrong with alignment\n";
        //    throw -1;
        //}
        m_data = (uint8_t*)(addr + moveFwdBy);
    }

    bool operator==(const BinaryMatrix& rhs) const
    {
        if (0 == memcmp(m_data, rhs.m_data, s_nUsedBytes))
            return true;
        bool eq = true;
        for (size_t i = 0; i < s_nBitRows && eq; ++i)
            for (size_t j = 0; j < s_nBitCols && eq; ++j) {
                eq &= getBit(i, j) == rhs.getBit(i, j);
                if (!eq) {
                    std::cout << "diff at: (" << i << "," << j << "): " << getBit(i, j) << " vs " << rhs.getBit(i, j) << "\n";
                }
            }
        return eq;
    }

    // set bit at row r and column c
    // r and c are bit indices
    void setBit(size_t bitRowIndex, size_t bitColIndex)
    {
        auto [byteIndex, bitOffset] = bitPos(bitColIndex);
        getByte(bitRowIndex, byteIndex) |= bitmask(bitOffset);
    }

    // get bit at row r and column c
    // r and c are bit indices
    bool getBit(size_t bitRowIndex, size_t bitColIndex) const
    {
        auto [byteIndex, bitOffset] = bitPos(bitColIndex);
        return getByte(bitRowIndex, byteIndex) & bitmask(bitOffset);
    }

    void printBits(size_t bitRowIndex, size_t bitColIndex, size_t nRows, size_t nCols, bool matlabIndexStyle = true) const
    {
        std::cout << "bit columns " << bitColIndex + matlabIndexStyle << " to " << bitColIndex + nCols - (1-matlabIndexStyle) << "\n";
        for (size_t i = 0; i < nRows; ++i) {
            std::cout << std::setw(5) << bitRowIndex + i + matlabIndexStyle << ": ";
            for (size_t j = 0; j < nCols; ++j)
                std::cout << getBit(bitRowIndex + i, bitColIndex + j);
            std::cout << std::endl;
        }
    }

    // count non-zero bit in a range
    size_t nnz(size_t bitRowIndex, size_t bitColIndex, size_t nRows, size_t nCols) const
    {
        size_t n = 0;
        for (size_t i = 0; i < nRows; ++i) {
            for (size_t j = 0; j < nCols; ++j)
                n += getBit(bitRowIndex + i, bitColIndex + j);
        }
        return n;
    }

    void toMatlab(const char *path, const char *name) const
    {
        std::ostringstream fn;
        fn << path << name << ".m";
        std::ofstream os(fn.str());
        std::vector<size_t> I, J, Z;
        os << "function y=" << name << "()\n";
        for (size_t i = 0; i < s_nBitRows; ++i) {
            for (size_t j = 0; j < s_nBitCols; ++j)
                if (getBit(i, j)) {
                    I.push_back(i + 1);
                    J.push_back(j + 1);
                    Z.push_back(1);
                }
        }
        os << "I=[";
        for (auto i : I) os << i << " ";
        os << "];\n";
        os << "J=[";
        for (auto j : J) os << j << " ";
        os << "];\n";
        os << "Z=[";
        for (auto z : Z) os << z << "==1 ";
        os << "];\n";
        os << "y = sparse(I, J, Z, " << s_nBitRows << ", " << s_nBitCols << ");\n";
        os << "end\n";
    }

    template <typename OS>
    void txtRowEncoder(OS& os, void (*encoder)(std::string&, const std::string&)) const
    {
        // FIXME: too many copies
        std::string binStr;
        binStr.resize(s_nBytesPerRow * s_nBitRows);
        char* pout = &binStr[0];
        for (size_t r = 0; r < s_nBitRows; ++r) {
            auto prow = (const char*)rowBegin(r);
            std::copy(prow, prow + s_nBytesPerRow, pout);
            pout += s_nBytesPerRow;
        }
        
        std::string encodedStr;
        (*encoder)(encodedStr, binStr);
        os << encodedStr;
    }

    template <typename IS>
    void txtRowDecoder(IS& is, void (*decoder)(std::string&, const std::string&))
    {
        std::string encodedStr, binStr;
        std::getline(is, encodedStr);
        (*decoder)(binStr, encodedStr);

        if (binStr.length() != s_nBytesPerRow * s_nBitRows) {
            //std::cout << "row " << r << " has bad length stream, got " << binStr.length() << ", expect " << s_nBytesPerRow << "\n";
            throw std::invalid_argument("bad stream length");
        }

        const uint8_t* p = (const uint8_t*)&binStr[0];
        for (size_t r = 0; r < s_nBitRows; ++r, p += s_nBytesPerRow)
            std::copy(p, p + s_nBytesPerRow, (char *) rowBegin(r));
    }

    template <typename OS>
    void toBase64(OS& os) const
    {
        txtRowEncoder(os, &Encoder::textToBase64);
    }

    template <typename OS>
    void fromBase64(OS& os)
    {
        txtRowDecoder(os, &Encoder::base64ToText);
    }


    template <typename OS>
    void toHex(OS& os) const
    {
        txtRowEncoder(os, &Encoder::textToHex);
    }

    template <typename OS>
    void fromHex(OS& os)
    {
        txtRowDecoder(os, &Encoder::hexToText);
    }

    size_t nnz() const
    {
        size_t n = 0;
        static_assert(s_nUsedBytes % sizeof(uint64_t) == 0);
        for (const int64_t *p = (const int64_t*)m_data, * const pend = p + (s_nUsedBytes / sizeof(*p)); p != pend; ++p)
            n += popcnt(*p);
        return n;
    }

    auto sparseIndex() const
    {
        size_t n = nnz();
        return std::make_pair(n, double(n) / (double(s_nBitCols) * double(s_nBitRows)));
    }

    // returns the word of type U containing the bit in position colBitIndex
    template <typename U>
    U getWordU(size_t rowBitIndex, size_t colBitIndex) const
    {
        // requires that colBitIndex is a multiple of sizeof(T)
        auto p = (const U*)rowBegin(rowBitIndex);
        return p[colBitIndex / (8 * sizeof(U))];
    }

    const uint8_t* rowBegin(size_t rowIndex) const { return m_data + rowIndex * s_nBytesPerPaddedRow; }
    uint8_t* rowBegin(size_t rowIndex) { return m_data + rowIndex * s_nBytesPerPaddedRow; }

    void resetZero() { std::fill(m_storage.begin(), m_storage.end(), 0); }

    //bool multiplyRowByCol(size_t r, size_t c) const
    //{
    //    bool result = false;
    //    for (size_t j = 0; j < s_nBitCols; ++j)
    //        result ^= getBit(r, j) && getBit(j, c);
    //    return result;
    //}

    // initialize from a file saved in matlab with the commands:
    //    [I,J,V] = find(sparse(M));
    //    U = [I,J];
    //    save('filename','U', '-ascii');
    void fromMatlabSparseFile(const char* filename)
    {
        std::ifstream is(filename);
        if (!is.good()) {
            std::cout << "failure opening " << filename << "\n";
            std::exit(-1);
        }
        double r, c;
        while (true) {
            is >> r >> c;
            if (!is.eof())
                setBit((size_t)r - 1, (size_t)c - 1);
            else
                break;
        }
    }

    void initRand()
    {
        for (size_t r = 0; r < s_nBitRows; ++r) {
            uint8_t* pr = rowBegin(r);
            size_t c = 0;
            {
                typedef uint16_t word_t;
                for (; c + sizeof(word_t) * 8 < s_nBitCols; c += sizeof(word_t) * 8) {
                    word_t rnd = word_t(rand()) & word_t(-1);
                    ((word_t*)rowBegin(r))[c / (sizeof(word_t) * 8)] = rnd;
                }
            }
            {
                typedef uint8_t word_t;
                for (; c + sizeof(word_t) * 8 < s_nBitCols; c += sizeof(word_t) * 8) {
                    word_t rnd = word_t(rand()) & word_t(-1);
                    ((word_t*)rowBegin(r))[c / (sizeof(word_t) * 8)] = rnd;
                }
            }
            for (;  c < s_nBitCols; ++c) {
                auto rnd = rand() % 2;
                if (rnd)
                    setBit(r, c);
            }
        }
    }

};
