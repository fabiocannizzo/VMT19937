#include "RandGen.h"
#include "cli_args.h"
#include "polynomial_jump.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <iomanip>

/**
 * characteristic_poly_finder - Derives the characteristic polynomial P(x) for GF(2) linear generators.
 * 
 * This utility uses the Berlekamp-Massey algorithm to find the shortest linear feedback shift register
 * (LFSR) that generates a given sequence. For a generator with state size N, observing 2*N bits
 * is sufficient to uniquely identify the characteristic polynomial.
 * 
 * References:
 * [1] Elwyn Berlekamp (1967). "Nonbinary BCH decoding". International Symposium on Information Theory.
 * [2] James L. Massey (1969). "Shift-register synthesis and BCH decoding". IEEE Transactions on Information Theory.
 * [3] Matsumoto, M.; Nishimura, T. (1998). "Mersenne twister: a 623-dimensionally equidistributed 
 *     uniform pseudo-random number generator". ACM Transactions on Modeling and Computer Simulation.
 */

using namespace std;
using namespace xvmt;
using namespace xvmt::details;

// Thin subclass that exposes one-bit-per-block sampling for BM on SFMT.
// SFMT advances one 128-bit block per step(), so we must read exactly one bit
// per step() call to give BM a true LFSR sequence.
template <ISA Isa = SIMD_ISA>
struct SFMTBitGen : public SFMT19937Base<128, Isa> {
    using Base = SFMT19937Base<128, Isa>;
    static constexpr int s_N_local = Base::s_N;

    void init(const uint32_t* seeds, uint32_t n) {
        Base::reinitMainState(seeds, n);
        Base::reinitPointers();
    }

    // Advance by one 128-bit block and return bit 0 of word 0 of the updated block.
    int stepBit() {
        Base::step();
        int prevIdx = ((int)Base::m_step_idx - 1 + s_N_local) % s_N_local;
        return Base::m_state[prevIdx * 4] & 1;
    }
};

void print_usage() {
    cout << "Usage: characteristic_poly_finder.exe -g=<mt32|mt64|sfmt> [-o=<output_file>]" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "  -g  Generator type: mt32 (MT19937-32), mt64 (MT19937-64), or sfmt (SFMT19937)" << endl;
    cout << "  -o  Output filename (default: characteristic.<type>.hex)" << endl;
    cout << "  -h  Show this help message" << endl;
}

// Bitstream generation logic - observing the sequence of a linear projection of the state.
template <typename Gen>
std::vector<uint64_t> generate_bitstream(Gen& gen, size_t n_bits) {
    size_t n_words = (n_bits + 63) / 64;
    std::vector<uint64_t> stream(n_words, 0);
    
    for (size_t i = 0; i < n_bits; ++i) {
        // We observe the least significant bit of the generator's word output.
        // Since the generator is F2-linear, any linear combination of state bits
        // follows the same characteristic polynomial.
        uint32_t out = gen.genrand_uint32(); 
        if (out & 1) {
            stream[i / 64] |= (1ULL << (i % 64));
        }
    }
    return stream;
}

// Berlekamp-Massey algorithm implementation over GF(2).
// Returns the actual degree L of the found connection polynomial.
size_t berlekamp_massey(const std::vector<uint64_t>& stream, size_t n_bits, std::vector<uint64_t>& poly) {
    const size_t N = n_bits;
    const size_t n_words = (N + 64) / 64;
    
    std::vector<uint64_t> C(n_words, 0);
    std::vector<uint64_t> B(n_words, 0);
    C[0] = 1; // C(x) = 1
    B[0] = 1; // B(x) = 1
    
    int L = 0;
    int m = 1;

    for (int n = 0; n < (int)N; ++n) {
        // Calculate discrepancy d
        int d = (stream[n / 64] >> (n % 64)) & 1;
        for (int i = 1; i <= L; ++i) {
            if ((C[i / 64] >> (i % 64)) & 1) {
                int stream_idx = n - i;
                if ((stream[stream_idx / 64] >> (stream_idx % 64)) & 1) {
                    d ^= 1;
                }
            }
        }

        if (d == 0) {
            m++;
        } else {
            std::vector<uint64_t> T = C;
            
            // C(x) = C(x) + B(x) * x^m
            int word_shift = m / 64;
            int bit_shift = m % 64;
            
            if (bit_shift == 0) {
                for (int i = word_shift; i < (int)n_words; ++i) {
                    C[i] ^= B[i - word_shift];
                }
            } else {
                uint64_t carry = 0;
                for (int i = word_shift; i < (int)n_words; ++i) {
                    uint64_t val = B[i - word_shift];
                    C[i] ^= (val << bit_shift) | carry;
                    carry = val >> (64 - bit_shift);
                }
            }

            if (2 * L <= n) {
                L = n + 1 - L;
                B = T;
                m = 1;
            } else {
                m++;
            }
        }
    }
    poly = C;
    return (size_t)L;
}

int main(int argc, const char** argv) {
    ArgMap args = parseArgs(argc, argv);
    if (consumeArg(args, "h")) {
        print_usage();
        return 0;
    }

    string gentype = "mt32";
    if (!consumeArg(args, "g", false, gentype)) {
        cerr << "Error: Generator type (-g) is required." << endl;
        print_usage();
        return -1;
    }
    
    string outfile = "characteristic." + gentype + ".hex";
    consumeArg(args, "o", false, outfile);

    cout << "Finding characteristic polynomial for " << gentype << " ..." << endl;
    
    size_t N = 19937;
    size_t n_stream_bits = 2 * N + 128; // Extra bits for safety
    
    std::vector<uint64_t> stream;
    
    if (gentype == "mt32") {
        VMT19937<32> gen(1234, 0, nullptr, nullptr);
        stream = generate_bitstream(gen, n_stream_bits);
    } else if (gentype == "mt64") {
        // Must use genrand_uint64() — one call per recurrence step — so BM sees a single
        // GF(2)-linear projection. genrand_uint32() on a 64-bit generator returns alternate
        // 32-bit halves of each 64-bit word, interleaving two projections and causing BM
        // to find the wrong degree (same problem as SFMT, fixed analogously).
        XMT19937_64<> gen;
        gen.reinit(1234ULL, (const XMT19937_64<>::poly_t*)nullptr, nullptr);
        size_t n_words_stream = (n_stream_bits + 63) / 64;
        stream.assign(n_words_stream, 0);
        for (size_t i = 0; i < n_stream_bits; ++i) {
            if (gen.genrand_uint64() & 1)
                stream[i / 64] |= (1ULL << (i % 64));
        }
    } else if (gentype == "sfmt") {
        // SFMT advances one 128-bit block per recurrence step. To give BM a true
        // LFSR sequence we must sample ONE bit per block step, not per uint32 word.
        // Using genrand_uint32() rotates through 4 different linear projections per
        // block, producing a non-LFSR sequence that causes BM to find the wrong degree.
        SFMTBitGen<> gen;
        uint32_t seeds[] = {1, 2, 3, 4};
        gen.init(seeds, 4);
        size_t n_bm_words = (n_stream_bits + 63) / 64;
        stream.assign(n_bm_words, 0);
        for (size_t i = 0; i < n_stream_bits; ++i) {
            if (gen.stepBit())
                stream[i / 64] |= (1ULL << (i % 64));
        }
    } else {
        cerr << "Unknown generator type: " << gentype << endl;
        return -1;
    }
    
    cout << "Bitstream generated. Running Berlekamp-Massey..." << endl;
    auto start = chrono::high_resolution_clock::now();
    
    std::vector<uint64_t> poly_words;
    size_t L = berlekamp_massey(stream, n_stream_bits, poly_words);

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    cout << "Berlekamp-Massey finished in " << elapsed.count() << "s" << endl;
    cout << "BM degree found: " << L << " (expected " << N << ")" << endl;

    // BM returns connection polynomial C(x). We need characteristic polynomial P(x) = reverse(C).
    // P[k] = C[L-k] for k=0..L where L is the actual degree found.
    size_t nPWords = (L + 64) / 64;
    std::vector<uint64_t> p_words(nPWords, 0);
    for (size_t k = 0; k <= L; ++k) {
        size_t srcBit = L - k;
        uint64_t cBit = (poly_words[srcBit / 64] >> (srcBit % 64)) & 1;
        if (cBit) {
            p_words[k / 64] |= (1ULL << (k % 64));
        }
    }

    // Compute Hamming weight for validation (should be odd for irreducible polynomial)
    size_t hw = 0;
    for (uint64_t w : p_words) {
        uint64_t v = w;
        while (v) { hw += v & 1; v >>= 1; }
    }
    cout << "Characteristic polynomial Hamming weight: " << hw
         << (hw % 2 == 1 ? " (odd - ok)" : " (even - suspicious)") << endl;

    // Convert to hex string (LSB first)
    stringstream ss;
    ss << hex << setfill('0');
    for (uint64_t w : p_words) {
        for (int i = 0; i < 16; ++i) {
            ss << ((w >> (i * 4)) & 0xF);
        }
    }
    string full_hex = ss.str();
    // Trim trailing zeros
    size_t last = full_hex.find_last_not_of('0');
    if (last != string::npos) full_hex = full_hex.substr(0, last + 1);

    cout << "Characteristic polynomial length: " << full_hex.length() * 4 << " bits" << endl;
    cout << "First 20 nibbles: " << full_hex.substr(0, 20) << endl;
    cout << "Last  20 nibbles: " << full_hex.substr(full_hex.size() > 20 ? full_hex.size() - 20 : 0) << endl;

    ofstream ofs(outfile);
    if (!ofs) {
        cerr << "Error: Could not open output file " << outfile << endl;
        return -1;
    }
    ofs << full_hex << endl;
    cout << "Saved to " << outfile << endl;
    
    return 0;
}
