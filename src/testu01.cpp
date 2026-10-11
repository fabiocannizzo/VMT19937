extern "C" {
#   include "TestU01.h"
}

#define SIMD_EMULATION

#include "RandGen.h"
#include "TestUtils.h"
#include "cli_args.h"

#include <cstdint>
#include <cstdio>
#include <functional>
#include <memory>
#include <sstream>
#include <string>
#include <iostream>
#include <cmath>

using namespace std;
using namespace xvmt;
using namespace xvmt::details;

const uint32_t g_seedlength = 4;
const uint32_t g_seedinit[g_seedlength] = { 0x123, 0x234, 0x345, 0x456 };
const uint64_t g_seedinit64[g_seedlength] = { 0x123, 0x234, 0x345, 0x456 };

enum Modes { SmallCrush = 0, Crush = 1, BigCrush = 2 };

// All generators use 128-bit hardware registers, whatever NBITS: wider virtual registers
// are emulated, which gives the same sequence whatever instruction set the compiler enables.
constexpr ISA g_isa = BitLenToIsa<128>::s_isa;

// log2 of the number of states M (a power of 2)
size_t log2States(size_t m)
{
    size_t log2m = 0;
    while ((size_t(1) << log2m) < m)
        ++log2m;
    return log2m;
}

// run the selected battery on a generator of 32-bit words
void runBattery(const std::string& name, unsigned int (*getRnd)(void), size_t mode)
{
    char genName[64];
    snprintf(genName, sizeof(genName), "%s", name.c_str());

    // create TestU01 generator wrapper
    unif01_Gen* gen = unif01_CreateExternGenBits(genName, getRnd);

    // Run the tests.
    switch (mode) {
        case SmallCrush:
            bbattery_SmallCrush(gen);
            break;
        case Crush:
            bbattery_Crush(gen);
            break;
        case BigCrush:
            bbattery_BigCrush(gen);
            break;
        default:
            MYASSERT(false, "Invalid mode " << mode);
    }

    // Clean up.
    unif01_DeleteExternGenBits(gen);
}

// VMT19937 (32 bit), states de-phased via jump matrices (jump ahead by 2^19937/M)
template <size_t NBITS>
struct TestRunner
{
    using gen_t = VMT19937<NBITS, QM_Scalar, g_isa>;

    static constexpr size_t s_m = NBITS / 32;

    static gen_t* s_genptr;

    static unsigned int getRnd()
    {
        return s_genptr->genrand_uint32();
    }

    static void run(size_t mode, const std::string& datDir)
    {
        std::unique_ptr<MT19937Matrix<32>> jump;
        if (s_m > 1) {
            std::ostringstream os;
            os << datDir << "/matrix/mt32/F" << (19937 - log2States(s_m)) << ".mt32.bits";
            jump = std::make_unique<MT19937Matrix<32>>(os.str());
        }

        gen_t g(g_seedinit, g_seedlength, 0, nullptr, jump.get());
        s_genptr = &g;

        std::ostringstream name;
        name << "VMT19937 (M=" << s_m << ")";
        runBattery(name.str(), getRnd, mode);
    }
};

template <size_t NBITS>
typename TestRunner<NBITS>::gen_t* TestRunner<NBITS>::s_genptr = nullptr;

// feeds TestU01 with the words of a generator with 32-bit output
template <typename Gen>
struct Feed32
{
    static Gen* s_genptr;

    static unsigned int getRnd()
    {
        return s_genptr->genrand_uint32();
    }
};

template <typename Gen>
Gen* Feed32<Gen>::s_genptr = nullptr;

// feeds TestU01 with the words of a generator with 64-bit output,
// each split in two 32-bit words: most significant half first
template <typename Gen>
struct Feed64
{
    static Gen* s_genptr;
    static uint64_t s_word;
    static bool s_hasLow;

    static unsigned int getRnd()
    {
        if (s_hasLow) {
            s_hasLow = false;
            return uint32_t(s_word);
        }
        s_word = s_genptr->genrand_uint64();
        s_hasLow = true;
        return uint32_t(s_word >> 32);
    }
};

template <typename Gen>
Gen* Feed64<Gen>::s_genptr = nullptr;
template <typename Gen>
uint64_t Feed64<Gen>::s_word = 0;
template <typename Gen>
bool Feed64<Gen>::s_hasLow = false;

// generators with states de-phased via jump polynomials (jump ahead by 2^19937/M)
template <typename Gen, template <typename> class Feed, typename SeedT>
void runPoly(const char* label, const char* tag, const SeedT* seeds, size_t mode, const std::string& datDir)
{
    using poly_t = typename Gen::poly_t;
    const size_t m = Gen::s_nStates;

    poly_t jump;
    if (m > 1) {
        std::ostringstream os;
        os << datDir << "/poly/" << tag << "/J" << (19937 - log2States(m)) << "." << tag << ".bits";
        jump.fromBinFile(os.str());
    }

    Gen g(seeds, g_seedlength, (const poly_t*)nullptr, (const poly_t*)(m > 1 ? &jump : nullptr));
    Feed<Gen>::s_genptr = &g;

    std::ostringstream name;
    name << label << " (M=" << m << ")";
    runBattery(name.str(), Feed<Gen>::getRnd, mode);
}

void usage()
{
    std::cerr
        << "Invalid command line arguments\n"
        << "Syntax:\n"
        << "   testu01 -b=<nbits> -m=<mode{0,1,2}> [-g=<gen{mt32,mt64,sfmt}>] [-d=<datdir>] [-wait]\n"
        << "   nbits determines the number of states M; with one state we get the base generator:\n"
        << "      mt32: 32, 128, 256, 512\n"
        << "      mt64: 64, 128, 256, 512\n"
        << "      sfmt: 128, 256, 512\n"
        << "Example:\n"
        << "   testu01 -b=128 -m=0\n"
        << "   testu01 -g=mt64 -b=512 -m=2\n";
    THROW("");
}


int main(int argc, const char** argv)
{
    ArgMap args = parseArgs(argc, argv);
    waitForDebugger(args);

    // parse command line arguments
    size_t nBits = 0;
    size_t mode = 0;
    std::string genType = "mt32";
    std::string datDir = "./dat";
    try {
        consumeArg(args, "b", true, nBits);
        consumeArg(args, "m", false, mode);
        consumeArg(args, "g", false, genType);
        consumeArg(args, "d", false, datDir);

        if (!args.empty())
            usage();
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        usage();
    }

    try {
        if (genType == "mt32") {
            switch (nBits) {
                case 32:
                    TestRunner<32>::run(mode, datDir);
                    break;
                case 128:
                    TestRunner<128>::run(mode, datDir);
                    break;
                case 256:
                    TestRunner<256>::run(mode, datDir);
                    break;
                case 512:
                    TestRunner<512>::run(mode, datDir);
                    break;
                default:
                    MYASSERT(false, "Invalid number of bits " << nBits);
            }
        }
        else if (genType == "mt64") {
            switch (nBits) {
                case 64: // single state: numerically identical to the original MT19937-64
                    runPoly<XMT19937_64<g_isa>, Feed64>("MT19937-64", "mt64", g_seedinit64, mode, datDir);
                    break;
                case 128:
                    runPoly<VMT19937_64<128, false, g_isa>, Feed64>("VMT19937-64", "mt64", g_seedinit64, mode, datDir);
                    break;
                case 256:
                    runPoly<VMT19937_64<256, false, g_isa>, Feed64>("VMT19937-64", "mt64", g_seedinit64, mode, datDir);
                    break;
                case 512:
                    runPoly<VMT19937_64<512, false, g_isa>, Feed64>("VMT19937-64", "mt64", g_seedinit64, mode, datDir);
                    break;
                default:
                    MYASSERT(false, "Invalid number of bits " << nBits);
            }
        }
        else if (genType == "sfmt") {
            switch (nBits) {
                case 128: // single state: numerically identical to the original SFMT19937
                    runPoly<VSFMT19937<128, false, g_isa>, Feed32>("SFMT19937", "sfmt", g_seedinit, mode, datDir);
                    break;
                case 256:
                    runPoly<VSFMT19937<256, false, g_isa>, Feed32>("VSFMT19937", "sfmt", g_seedinit, mode, datDir);
                    break;
                case 512:
                    runPoly<VSFMT19937<512, false, g_isa>, Feed32>("VSFMT19937", "sfmt", g_seedinit, mode, datDir);
                    break;
                default:
                    MYASSERT(false, "Invalid number of bits " << nBits);
            }
        }
        else {
            MYASSERT(false, "Invalid generator " << genType);
        }
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return -1;
    }
    catch (...) {
        std::cerr << "Error: unknown exception\n";
        return -1;
    }

    return 0;
}
