extern "C" {
#   include "TestU01.h"
}

#define SIMD_EMULATION

#include "RandGen.h"
#include "TestUtils.h"
#include "cli_args.h"

#include <cstdint>
#include <functional>
#include <string>
#include <iostream>
#include <cmath>

using namespace std;
using namespace xvmt;
using namespace xvmt::details;

MT19937Matrix<32> g_j19933(std::string("./dat/matrix/mt32/F19933.mt32.bits"));
MT19937Matrix<32> g_j19934(std::string("./dat/matrix/mt32/F19934.mt32.bits"));
MT19937Matrix<32> g_j19935(std::string("./dat/matrix/mt32/F19935.mt32.bits"));

MT19937Matrix<32>* g_pjump[4] = { nullptr, &g_j19935, &g_j19934, &g_j19933 };
char g_genNames[4][64] = { "VMT19937 (M=1)", "VMT19937 (M=4)", "VMT19937 (M=8)", "VMT19937 (M=16)"};

const uint32_t g_seedlength = 4;
const uint32_t g_seedinit[g_seedlength] = { 0x123, 0x234, 0x345, 0x456 };

enum Modes { SmallCrush = 0, Crush = 1, BigCrush = 2 };

template <size_t NBITS>
struct TestRunner
{
    using gen_t = VMT19937<NBITS, QM_Scalar>;

    static constexpr size_t s_m = NBITS / 32;
    static constexpr size_t s_log2 = s_m == 16 ? 4 : // log2(M) is not constexpr until C++26
                                   s_m == 8  ? 3 :
                                   s_m == 4  ? 2 :
                                   s_m == 2  ? 1 : 0;
    static constexpr size_t s_arrayIndex = s_m == 1 ? 0 : s_log2 - 1;

    static gen_t* s_genptr;

    static uint32_t getRnd()
    {
        return s_genptr->genrand_uint32();
    }

    static void run(size_t mode)
    {
        gen_t g(g_seedinit, g_seedlength, 0, nullptr, g_pjump[s_arrayIndex]);
        s_genptr = &g;

        // create TestU01 generator wrapper
        unif01_Gen* gen = unif01_CreateExternGenBits(g_genNames[s_arrayIndex], getRnd);

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
};

template <size_t NBITS>
typename TestRunner<NBITS>::gen_t* TestRunner<NBITS>::s_genptr = nullptr;

void usage()
{
    std::cerr
        << "Invalid command line arguments\n"
        << "Syntax:\n"
        << "   testu01 -b=<nbits> -m=<mode{0,1,2}> [-wait]\n"
        << "Example:\n"
        << "   testu01 -b=128 -m=0\n";
    THROW("");
}


int main(int argc, const char** argv)
{
    ArgMap args = parseArgs(argc, argv);
    waitForDebugger(args);

    // parse command line arguments
    size_t nBits = 0;
    size_t mode = 0;
    try {
        consumeArg(args, "b", true, nBits);
        consumeArg(args, "m", false, mode);

        if (!args.empty())
            usage();
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        usage();
    }

    try {
        switch (nBits) {
            case 32:
                TestRunner<32>::run(mode);
                break;
            case 128:
                TestRunner<128>::run(mode);
                break;
            case 256:
                TestRunner<256>::run(mode);
                break;
            case 512:
                TestRunner<512>::run(mode);
                break;
            default:
                MYASSERT(false, "Invalid number of bits " << nBits);
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
