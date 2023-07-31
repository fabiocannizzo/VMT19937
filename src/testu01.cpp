extern "C" {
#   include "TestU01.h"
}

#define SIMD_EMULATION

#include "MSMT19937.h"

#include <cstdint>
#include <functional>
#include <string>
#include <iostream>
#include <cmath>

using namespace std;
//using namespace std::placeholders;

MT19937Matrix j19933(std::string("./dat/F19933.bits"));
MT19937Matrix j19934(std::string("./dat/F19934.bits"));
MT19937Matrix j19935(std::string("./dat/F19935.bits"));

MT19937Matrix* pjump[4] = { nullptr, &j19935, &j19934, &j19933 };
char genNames[4][64] = { "MS-MT1937 (M=1)", "MS-MT1937 (M=4)", "MS-MT1937 (M=8)", "MS-MT1937 (M=16)"};

const uint32_t seedlength = 4;
const uint32_t seedinit[seedlength] = { 0x123, 0x234, 0x345, 0x456 };

enum Modes { SmallCrush = 0, Crush = 1, BigCrush = 2 };

template <size_t NBITS>
struct TestRunner
{
    typedef MSMT19937<NBITS, QM_Scalar> gen_t;

    static const size_t M = NBITS / 32;
    static const size_t ArrayIndex = M == 1 ? 0 : std::log2(M) - 1;

    static gen_t* s_genptr;

    static uint32_t getRnd()
    {
        return s_genptr->genrand_uint32();
    }

    static void run(size_t mode)
    {
        gen_t g(seedinit, seedlength, 0, nullptr, pjump[ArrayIndex]);
        s_genptr = &g;

        // create TestU01 generator wrapper
        unif01_Gen* gen = unif01_CreateExternGenBits(genNames[ArrayIndex], getRnd);

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
        << "   testu01 -b nbits -m mode{0,1,2}\n"
        << "Example:\n"
        << "   testu01 -b 128\n";
    THROW("");
}


int main(int argc, const char** argv)
{
    // parse command line arguments
    size_t nBits = 0;
    size_t mode = 0;
    if (argc % 2 == 0)
        usage();
    for (int i = 1; i < argc; i += 2) {
        string key(argv[i]);
        string value(argv[i + 1]);
        if (key == "-b")
            nBits = atoi(value.c_str());
        else if (key == "-m")
            mode = atoi(value.c_str());
        else {
            usage();
        }
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
