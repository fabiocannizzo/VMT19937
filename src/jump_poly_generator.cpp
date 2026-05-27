#include "polynomial_jump.h"
#include "cli_args.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <thread>
#include <mutex>
#include <chrono>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <filesystem>

/**
 * jump_poly_generator - Calculates jump polynomials g(x) = x^J mod P(x) for GF(2) generators.
 *
 * Multi-target mode (-t) chains squarings from the nearest existing file, so computing
 * J_{N+1} from J_N costs exactly one squaring (exploiting x^{2^N} = (x^{2^{N-1}})^2).
 *
 * References:
 * [1] Haramoto, H., Matsumoto, M., L'Ecuyer, P. (2008). "A Fast Jump Ahead Algorithm for
 *     Linear Recurrences in a Polynomial Space". SETA 2008.
 * [2] Knuth, D. E. (1997). "The Art of Computer Programming, Volume 2: Seminumerical Algorithms".
 */

using namespace std;
using namespace xvmt;
using namespace xvmt::details;

typedef Polynomial<32768> Poly19937;
typedef Polynomial<65536> PolySquared;

void print_usage() {
    cout << "Usage:\n"
         << "  Multi-target: jump_poly_generator.exe -g=<gentype> -p=<char_poly_hex> [-t=<n1,n2,...>] [-d=<dir>] [-j=<nthreads>]\n"
         << "  Single:       jump_poly_generator.exe [-n=<exp> | -s=<step>] -g=<gentype> -p=<char_poly_hex> [-o=<file>]\n"
         << "  Chain:        jump_poly_generator.exe -chain=<k> -i=<input> -p=<char_poly_hex> [-o=<file>]\n"
         << "\nOptions:\n"
         << "  -g     Generator type: mt32, mt64, sfmt (default: mt32)\n"
         << "         SFMT has 4 output words per recurrence step (power2=2);\n"
         << "         -n=N means 2^N output words, stored as x^{2^(N-2)} mod P internally.\n"
         << "  -p     Characteristic polynomial file (.hex format)\n"
         << "  -t     Comma-separated target exponents N (output-word units); default 9,100,19933,19934,19935,19936\n"
         << "         Chains squarings optimally: J_{N+1} = square(J_N), one squaring per step.\n"
         << "         Scans -d directory for existing J-files to resume from.\n"
         << "  -d     Output directory for multi-target mode (default: ./dat/<gentype>/)\n"
         << "  -n     Single-target: J = 2^n output words\n"
         << "  -s     Single-target: J = decimal step (output words; must be divisible by 4 for sfmt)\n"
         << "  -chain Apply k additional squarings to the polynomial loaded from -i\n"
         << "  -i     Input polynomial for -chain mode (.bits)\n"
         << "  -o     Output file for single-target / chain mode (.bits)\n"
         << "  -j     Threads (default: hardware concurrency)\n"
         << "  -h     Show this help\n";
}

Poly19937 poly_sq_mod(const Poly19937& in, const Poly19937& P, size_t deg) {
    PolySquared squared;
    PolyOps::square(in, squared);
    if (deg == 19937) PolyOps::reduce<19937>(squared, P);
    else if (deg == 19968) PolyOps::reduce<19968>(squared, P);
    else { cerr << "Unsupported polynomial degree " << deg << "\n"; exit(1); }

    alignas(64) uint32_t words[65536 / 32];
    squared.m_data.store(words);
    Poly19937 res;
    res.m_data = Poly19937::Reg(words);
    return res;
}

Poly19937 poly_mul_mod(const Poly19937& a, const Poly19937& b, const Poly19937& P, size_t deg) {
    Poly19937 res;
    res.resetZero();
    Poly19937 current_b = b;
    for (size_t i = 0; i < deg; ++i) {
        if (a.getBit(i)) res ^= current_b;
        bool shift_out = current_b.getBit(deg - 1);
        current_b = current_b << 1;
        if (shift_out) current_b ^= P;
    }
    return res;
}

string mkJFileName(const string& dir, size_t n, const string& gentype) {
    ostringstream os;
    os << dir << "J" << setw(5) << setfill('0') << n << "." << gentype << ".bits";
    return os.str();
}

Poly19937 load_char_poly(const string& char_file, size_t& deg) {
    ifstream ifs(char_file);
    if (!ifs) { cerr << "Error: Cannot open characteristic polynomial file: " << char_file << "\n"; exit(1); }
    string line;
    while (getline(ifs, line) && (line.empty() || line[0] == '#'));
    Poly19937 P;
    P.fromString(line);
    deg = P.degree();
    cout << "Polynomial degree: " << deg << "\n";
    return P;
}

// Multi-target mode: compute all targets by chaining squarings from the nearest existing file.
// Targets and all on-disk N values are in output-word units; internal poly exponent = N - power2.
// One squaring maps J_N -> J_{N+1} regardless of power2 (the offset cancels out in differences).
void run_multi_target(const string& gentype, size_t power2, const Poly19937& P, size_t deg,
                      const string& dir, const set<size_t>& targets) {
    // Seed onDisk with the implicit base (x^1 = J_{power2}, never stored)
    set<size_t> onDisk;
    onDisk.insert(power2);

    const string ext = "." + gentype + ".bits";
    for (const auto& entry : filesystem::directory_iterator(dir)) {
        if (!filesystem::is_regular_file(entry)) continue;
        string name = entry.path().filename().string();
        if (name.size() < 1 + 5 + ext.size() || name[0] != 'J') continue;
        if (name.substr(name.size() - ext.size()) != ext) continue;
        size_t n = (size_t)atoi(name.substr(1, 5).c_str());
        if (n >= power2) onDisk.insert(n);
    }

    size_t currentN = SIZE_MAX;
    Poly19937 current;
    auto wallStart = chrono::high_resolution_clock::now();
    size_t nSquarings = 0;

    for (size_t t : targets) {
        // Find closest existing N <= t
        auto it = onDisk.upper_bound(t);
        if (it == onDisk.begin()) { cerr << "No reachable starting point for target " << t << "\n"; continue; }
        --it;
        size_t s = *it;

        if (s == t) { cout << "J" << t << " already exists, skipping.\n"; continue; }

        // Load starting point if different from what's already in memory
        if (currentN != s) {
            if (s == power2) {
                current.resetZero();
                current.setBit(1, true); // x^1 = J_{power2}
                cout << "Starting from base polynomial x^1\n";
            } else {
                string fn = mkJFileName(dir, s, gentype);
                cout << "Loading J" << s << " from " << fn << "\n";
                ifstream is(fn, ios::binary);
                if (!is) { cerr << "Error: Cannot open " << fn << "\n"; continue; }
                current.fromBin(is);
            }
            currentN = s;
        }

        // Chain squarings from s to t
        for (size_t i = s + 1; i <= t; ++i) {
            auto t0 = chrono::high_resolution_clock::now();
            current = poly_sq_mod(current, P, deg);
            auto t1 = chrono::high_resolution_clock::now();
            double sq_s = chrono::duration<double>(t1 - t0).count();
            ++nSquarings;
            double elapsed = chrono::duration<double>(t1 - wallStart).count();
            cout << "  J" << i;
            if (i == t) cout << " [TARGET]";
            cout << "  squaring: " << fixed << setprecision(2) << sq_s << "s"
                 << "  elapsed: " << fixed << setprecision(1) << elapsed << "s\n";
            currentN = i;
        }

        string fn = mkJFileName(dir, t, gentype);
        ofstream os(fn, ios::binary);
        if (!os) { cerr << "Error: Cannot write " << fn << "\n"; continue; }
        current.toBin(os);
        onDisk.insert(t);
        cout << "Saved: " << fn << "\n";
    }

    double total = chrono::duration<double>(chrono::high_resolution_clock::now() - wallStart).count();
    cout << "Done. Total squarings: " << nSquarings << "  Total time: " << fixed << setprecision(1) << total << "s\n";
}

int main(int argc, const char** argv) {
    ArgMap args = parseArgs(argc, argv);
    if (consumeArg(args, "h") || argc == 1) { print_usage(); return 0; }

    // Generator type and power2
    size_t power2 = 0;
    string gentype = "mt32";
    string gentypeArg;
    if (consumeArg(args, "g", false, gentypeArg)) {
        gentype = gentypeArg;
        if (gentype == "sfmt") power2 = 2;
        else if (gentype == "mt32" || gentype == "mt64") power2 = 0;
        else { cerr << "Error: Unknown generator type '" << gentype << "'. Use mt32, mt64, or sfmt.\n"; return -1; }
    }

    string char_file;
    if (!consumeArg(args, "p", true, char_file)) {
        cerr << "Error: Characteristic polynomial file (-p) is required.\n"; return -1;
    }

    // --- Multi-target mode ---
    string targets_str = "9,100,19933,19934,19935,19936";
    string dir_arg;
    bool has_t = consumeArg(args, "t", false, targets_str);
    bool has_d = consumeArg(args, "d", false, dir_arg);

    // Detect multi-target mode: activated when -t or -d is given, or neither -n/-s/-chain present
    size_t exponent_peek = 0;
    string step_peek;
    size_t chain_peek = 0;
    bool has_n = consumeArg(args, "n", false, exponent_peek);
    bool has_s = consumeArg(args, "s", false, step_peek);
    bool has_chain = consumeArg(args, "chain", false, chain_peek);

    bool multi_mode = has_t || has_d || (!has_n && !has_s && !has_chain);

    if (multi_mode) {
        string dir = has_d ? dir_arg : ("./dat/" + gentype + "/");
        if (!dir.empty() && dir.back() != '/') dir += '/';
        if (!filesystem::exists(dir)) {
            cerr << "Error: Output directory does not exist: " << dir << "\n"; return -1;
        }

        set<size_t> targets;
        {
            stringstream ss(targets_str);
            string tok;
            while (getline(ss, tok, ',')) {
                size_t n = (size_t)atoi(tok.c_str());
                if (n < power2) { cerr << "Warning: target " << n << " < power2=" << power2 << ", skipping.\n"; continue; }
                targets.insert(n);
            }
        }
        if (targets.empty()) { cerr << "Error: No valid targets.\n"; return -1; }

        cout << "Multi-target mode: generator=" << gentype << "  dir=" << dir << "  targets=";
        for (auto it = targets.begin(); it != targets.end(); ++it)
            cout << *it << (next(it) == targets.end() ? "" : ",");
        cout << "\n";

        size_t deg;
        Poly19937 P = load_char_poly(char_file, deg);
        run_multi_target(gentype, power2, P, deg, dir, targets);
        return 0;
    }

    // --- Single / chain mode ---
    size_t deg;
    Poly19937 P = load_char_poly(char_file, deg);

    string outfile;
    consumeArg(args, "o", false, outfile);

    Poly19937 result;

    if (has_chain) {
        string chain_input;
        if (!consumeArg(args, "i", false, chain_input)) {
            cerr << "Error: -chain requires -i=<input_file>\n"; return -1;
        }
        cout << "Chain mode: applying " << chain_peek << " squaring(s) to " << chain_input << " ...\n";
        auto t0 = chrono::high_resolution_clock::now();
        ifstream ichain(chain_input, ios::binary);
        if (!ichain) { cerr << "Error: Cannot open " << chain_input << "\n"; return -1; }
        result.fromBin(ichain);
        for (size_t i = 0; i < chain_peek; ++i)
            result = poly_sq_mod(result, P, deg);
        double elapsed = chrono::duration<double>(chrono::high_resolution_clock::now() - t0).count();
        cout << "Done in " << fixed << setprecision(3) << elapsed << "s\n";
    } else if (has_n) {
        size_t exponent = exponent_peek;
        if (exponent < power2) {
            cerr << "Error: -n=" << exponent << " too small for -g=" << gentype << " (min " << power2 << ")\n";
            return -1;
        }
        size_t poly_exp = exponent - power2;
        cout << "Computing J = 2^" << exponent << " output words"
             << (power2 > 0 ? " (poly exp = 2^" + to_string(poly_exp) + " block steps)" : "")
             << " ...\n";
        auto t0 = chrono::high_resolution_clock::now();
        result.setBit(1, true);
        for (size_t i = 0; i < poly_exp; ++i)
            result = poly_sq_mod(result, P, deg);
        double elapsed = chrono::duration<double>(chrono::high_resolution_clock::now() - t0).count();
        cout << "Done in " << fixed << setprecision(3) << elapsed << "s\n";
    } else {
        // -s mode
        string step_str = step_peek;
        if (power2 > 0) {
            int div = 1 << power2;
            int rem = 0;
            for (char c : step_str) rem = (rem * 10 + (c - '0')) % div;
            if (rem != 0) {
                cerr << "Error: step '" << step_str << "' not divisible by " << div
                     << " (required for -g=" << gentype << ")\n";
                return -1;
            }
            for (size_t p = 0; p < power2; ++p) {
                int carry = 0;
                string next;
                for (char c : step_str) {
                    int val = carry * 10 + (c - '0');
                    next += char('0' + val / 2);
                    carry = val % 2;
                }
                size_t first = next.find_first_not_of('0');
                step_str = (first == string::npos) ? "0" : next.substr(first);
            }
        }
        cout << "Computing J = " << step_peek
             << (power2 > 0 ? " output words (" + step_str + " block steps)" : " output words")
             << " ...\n";
        auto t0 = chrono::high_resolution_clock::now();

        vector<uint8_t> step_bits;
        string current = step_str;
        while (current != "0" && !current.empty()) {
            int rem = 0;
            string next;
            for (char c : current) {
                int val = rem * 10 + (c - '0');
                next += to_string(val / 2);
                rem = val % 2;
            }
            step_bits.push_back((uint8_t)rem);
            size_t first = next.find_first_not_of('0');
            current = (first == string::npos) ? "0" : next.substr(first);
        }

        Poly19937 base;
        base.setBit(1, true);
        result.setBit(0, true);
        for (size_t i = 0; i < step_bits.size(); ++i) {
            if (step_bits[i]) result = poly_mul_mod(result, base, P, deg);
            if (i < step_bits.size() - 1) base = poly_sq_mod(base, P, deg);
        }
        double elapsed = chrono::duration<double>(chrono::high_resolution_clock::now() - t0).count();
        cout << "Done in " << fixed << setprecision(3) << elapsed << "s\n";
    }

    if (!outfile.empty()) {
        ofstream ofs(outfile, ios::binary);
        if (!ofs) { cerr << "Error: Cannot write " << outfile << "\n"; return -1; }
        result.toBin(ofs);
        cout << "Saved to " << outfile << "\n";
    } else {
        cout << "Jump polynomial (hex):\n" << result.toString() << "\n";
    }

    return 0;
}
