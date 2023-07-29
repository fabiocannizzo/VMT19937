extern "C" {
#   include "TestU01.h"
}

#include "jump_matrix.h"
#include "MSMT19937.h"

#include <cstdint>
#include <functional>
#include <string>
#include <iostream>

using namespace std;
//using namespace std::placeholders;

MT19937Matrix j19933("./dat/F19933.bits");
MT19937Matrix j19934("./dat/F19934.bits");
MT19937Matrix j19935("./dat/F19935.bits");

const uint32_t seedlength = 4;
const uint32_t seedinit[seedlength] = { 0x123, 0x234, 0x345, 0x456 };

std::function<uint32_t(void)> s_getrnd;

uint32_t getrnd()
{
    return s_getrnd();
}

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


    char name[64];
    std::unique_ptr<MSMT19937<32>> p32;
    std::unique_ptr<MSMT19937<128>> p128;

    try {

        switch (nBits) {
            case 32:
            {
                p32.reset(new MSMT19937<32>(seedinit, seedlength, nullptr, nullptr));
                strcpy(name, "MS-1-MT1937");
                s_getrnd = bind(&MSMT19937<32>::genrand_uint32, p32.get());
                break;
            }
            case 128:
            {
                p128.reset(new MSMT19937<128>(seedinit, seedlength, nullptr, &j19935));
                strcpy(name, "MS-4-MT1937");
                s_getrnd = bind(&MSMT19937<128>::genrand_uint32, p128.get());
                break;
            }
            default:
                MYASSERT(false, "Invalid number of bits " << nBits);
        };

        unif01_Gen* gen = unif01_CreateExternGenBits(name, getrnd);

        // Run the tests.
        switch (mode) {
            case 0:
                bbattery_SmallCrush(gen);
                break;
            case 1:
                bbattery_Crush(gen);
                break;
            case 2:
                bbattery_BigCrush(gen);
                break;
            default:
                MYASSERT(false, "Invalid mode " << mode);
        }

        // Clean up.
        unif01_DeleteExternGenBits(gen);
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
