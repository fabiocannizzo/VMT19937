#pragma once

#include "SIMD.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

namespace xvmt {
namespace details {

/**
 * @brief Represents a binary polynomial over GF(2) with a fixed maximum degree.
 * Specifically tuned for MT19937 characteristic polynomial of degree 19937.
 * MaxBits must be a multiple of HwBitLen and MaxBits/HwBitLen must be a power of 2.
 */
template <size_t MaxBits = 32768, ISA Isa = SIMD_ISA>
struct Polynomial {
    static constexpr size_t s_maxBits = MaxBits;
    static constexpr size_t s_n32 = MaxBits / 32;
    using Reg = SimdRegister<MaxBits, Isa>;

    Reg m_data;

    Polynomial() : m_data(Reg::zero()) {}
    explicit Polynomial(const Reg& r) : m_data(r) {}

    static Polynomial zero() { return Polynomial(); }

    void resetZero() { m_data = Reg::zero(); }

    bool isZero() const {
        return m_data.eq(Reg::zero());
    }

    bool operator==(const Polynomial& other) const {
        return m_data.eq(other.m_data);
    }

    Polynomial operator^(const Polynomial& other) const {
        return Polynomial(m_data ^ other.m_data);
    }

    Polynomial& operator^=(const Polynomial& other) {
        m_data = m_data ^ other.m_data;
        return *this;
    }

    Polynomial operator<<(size_t n) const {
        return Polynomial(m_data << n);
    }

    Polynomial operator>>(size_t n) const {
        return Polynomial(m_data >> n);
    }

    /**
     * @brief Get 32-bit word at index.
     */
    uint32_t getWord(size_t i) const {
        if (i >= s_n32) return 0;
        alignas(64) uint32_t words[s_n32];
        m_data.store(words);
        return words[i];
    }

    /**
     * @brief Set 32-bit word at index.
     */
    void setWord(size_t i, uint32_t val) {
        if (i >= s_n32) return;
        alignas(64) uint32_t words[s_n32];
        m_data.store(words);
        words[i] = val;
        m_data = Reg(words);
    }

    /**
     * @brief Load polynomial from hex string (SFMT format: nibbles, LSB first).
     */
    void fromString(const std::string& str) {
        resetZero();
        alignas(64) uint32_t words[s_n32] = {0};
        
        size_t p = 0;
        for (char c : str) {
            uint32_t nibble;
            if (c >= 'a' && c <= 'f') nibble = c - 'a' + 10;
            else if (c >= 'A' && c <= 'F') nibble = c - 'A' + 10;
            else if (c >= '0' && c <= '9') nibble = c - '0';
            else continue;

            for (int j = 0; j < 4; ++j) {
                if (nibble & (1 << j)) {
                    size_t wordIdx = (p + j) / 32;
                    size_t bitIdx = (p + j) % 32;
                    if (wordIdx < s_n32) {
                        words[wordIdx] |= (1U << bitIdx);
                    }
                }
            }
            p += 4;
        }
        m_data = Reg(words);
    }

    /**
     * @brief Convert polynomial to hex string (SFMT format).
     */
    std::string toString() const {
        alignas(64) uint32_t words[s_n32];
        m_data.store(words);

        std::stringstream ss;
        ss << std::hex;
        
        // Find actual degree to truncate output if needed, but SFMT usually prints full length
        size_t lastWord = s_n32;
        while (lastWord > 0 && words[lastWord - 1] == 0) lastWord--;
        if (lastWord == 0) return "0";

        for (size_t i = 0; i < s_n32; i += 1) {
            uint32_t w = words[i];
            for (int j = 0; j < 8; ++j) {
                ss << ((w >> (j * 4)) & 0xF);
            }
        }
        std::string res = ss.str();
        // Trim trailing zeros from the hex string (which are leading in the polynomial)
        size_t lastNonZero = res.find_last_not_of('0');
        if (lastNonZero == std::string::npos) return "0";
        return res.substr(0, lastNonZero + 1);
    }

    /**
     * @brief Save to binary file (~4 KB for 32768 bits).
     */
    void toBin(std::ostream& os) const {
        alignas(64) uint32_t words[s_n32];
        m_data.store(words);
        os.write(reinterpret_cast<const char*>(words), sizeof(words));
    }

    /**
     * @brief Load from binary file.
     */
    void fromBin(std::istream& is) {
        alignas(64) uint32_t words[s_n32];
        is.read(reinterpret_cast<char*>(words), sizeof(words));
        m_data = Reg(words);
    }

    /**
     * @brief Return the degree (index of highest set bit), or 0 for the zero polynomial.
     */
    size_t degree() const {
        alignas(64) uint32_t words[s_n32];
        m_data.store(words);
        for (int i = (int)s_n32 - 1; i >= 0; --i) {
            if (words[i] != 0) {
                uint32_t w = words[i];
                int bit = 31;
                while (!((w >> bit) & 1)) --bit;
                return (size_t)i * 32 + (size_t)bit;
            }
        }
        return 0;
    }

    /**
     * @brief Get bit at index.
     */
    bool getBit(size_t i) const {
        if (i >= MaxBits) return false;
        alignas(64) uint32_t words[s_n32];
        m_data.store(words);
        return (words[i / 32] >> (i % 32)) & 1;
    }

    /**
     * @brief Set bit at index.
     */
    void setBit(size_t i, bool val) {
        if (i >= MaxBits) return;
        alignas(64) uint32_t words[s_n32];
        m_data.store(words);
        if (val) words[i / 32] |= (1U << (i % 32));
        else words[i / 32] &= ~(1U << (i % 32));
        m_data = Reg(words);
    }
};

/**
 * @brief Helper for vectorized polynomial operations.
 */
struct PolyOps {
    /**
     * @brief Vectorized squaring of a polynomial A(x) -> A(x^2).
     * The result has twice the degree.
     */
    template <size_t N, size_t N2, ISA Isa>
    static void square(const Polynomial<N, Isa>& in, Polynomial<N2, Isa>& out) {
        static_assert(N2 >= 2 * N, "Output polynomial must be large enough");
        using Reg128 = SimdRegister<128, Isa>;
        
        alignas(64) uint32_t inWords[N / 32];
        in.m_data.store(inWords);

        alignas(64) uint32_t outWords[N2 / 32] = {0};

        for (size_t i = 0; i < N / 64; ++i) {
            uint64_t w = reinterpret_cast<const uint64_t*>(inWords)[i];
            Reg128 r_w(w);
            // Squaring 64-bit polynomial results in 128-bit polynomial
            Reg128 r_sq = Reg128::template clmul<0x00>(r_w, r_w);
            r_sq.template store<false>(outWords + i * 4);
        }
        out.m_data = typename Polynomial<N2, Isa>::Reg(outWords);
    }

    /**
     * @brief Modular reduction of A(x) by P(x).
     * P(x) is assumed to have degree DEG.
     */
    template <size_t DEG, size_t N, size_t PN, ISA Isa>
    static void reduce(Polynomial<N, Isa>& A, const Polynomial<PN, Isa>& P) {
        alignas(64) uint32_t aWords[N / 32];
        A.m_data.store(aWords);

        alignas(64) uint32_t pWords[PN / 32];
        P.m_data.store(pWords);

        constexpr size_t nA = N / 32;
        constexpr size_t nP = PN / 32;

        // Process bits from highest down to DEG, XOR-ing P shifted by (bitIdx-DEG) into A.
        // The Polynomial::operator<< is a per-lane 32-bit shift (no carry across word boundaries),
        // so we implement the correct bit shift using a scalar carry-chain loop.
        for (int i = (int)nA - 1; i >= (int)(DEG / 32); --i) {
            uint32_t w = aWords[i];
            if (w == 0) continue;
            const int jEnd = (i == (int)(DEG / 32)) ? (int)(DEG % 32) : 0;
            for (int j = 31; j >= jEnd; --j) {
                if (!((w >> j) & 1)) continue;
                const size_t shift = (size_t)i * 32 + (size_t)j - DEG;
                const size_t wOff = shift >> 5;
                const int bOff = (int)(shift & 31);
                if (bOff == 0) {
                    for (size_t pi = 0; pi < nP; ++pi)
                        if (pi + wOff < nA) aWords[pi + wOff] ^= pWords[pi];
                } else {
                    uint32_t carry = 0;
                    for (size_t pi = 0; pi < nP; ++pi) {
                        uint32_t val = (pWords[pi] << bOff) | carry;
                        carry = pWords[pi] >> (32 - bOff);
                        if (pi + wOff < nA) aWords[pi + wOff] ^= val;
                    }
                    if (wOff + nP < nA) aWords[wOff + nP] ^= carry;
                }
                w = aWords[i];
            }
        }
        A.m_data = typename Polynomial<N, Isa>::Reg(aWords);
    }
};

/**
 * @brief Applier that uses a polynomial bitmask to advance generator state.
 */
template <typename Generator>
struct PolynomialJumpApplier {
    using output_word_t = typename Generator::output_word_t;
    static constexpr size_t s_n32InFullState = Generator::s_n32InFullState;

    /**
     * @brief Advance state by jump polynomial g(x).
     * Uses stateToVector/vectorToState to match the binary vector representation.
     */
    static void apply(Generator& gen, const Polynomial<32768, Generator::s_isa>& poly, size_t stateIndex = 0) {
        Generator tempGen = gen;
        // nWords = state size in uint32 per logical state (one per stateIndex slot).
        // For MT: 624 * 1 = 624; for SFMT: 156 * 4 = 624. Both give 624 = 19968 bits.
        // The loop runs 19968 times to cover polynomials of degree up to 19968.
        const size_t nWords = Generator::s_n32InOneState;
        alignas(64) uint32_t resultVec[624 * 2] = {0};
        alignas(64) uint32_t currentVec[624 * 2];

        using XV = SimdRegister<128, Generator::s_isa>;

        for (size_t i = 0; i < nWords * 32; ++i) {
            if (poly.getBit(i)) {
                tempGen.stateToVector(stateIndex, currentVec);
                for (size_t s = 0; s < nWords; s += 4) {
                    XV r_res = XV::template load<false>(resultVec + s);
                    XV r_cur = XV::template load<false>(currentVec + s);
                    (r_res ^ r_cur).template store<false>(resultVec + s);
                }
            }
            tempGen.step();
        }
        gen.vectorToState(stateIndex, resultVec);
    }
};

} // namespace details
} // namespace xvmt
