extern "C" {
#   include "TestU01.h"
}

#include "jump_matrix.h"
#include "MSMT19937.h"

MT19937Matrix j19933("./dat/F19933.bits");
MT19937Matrix j19934("./dat/F19934.bits");
MT19937Matrix j19935("./dat/F19935.bits");

// Example PRNG: Xorshift 32

static unsigned int y = 2463534242U;

unsigned int xorshift (void)
{
    y ^= (y << 13);
    y ^= (y >> 17);
    return y ^= (y << 5);
}

int main()
{
    char name[] = "Xorshift 32";

    // Create TestU01 PRNG object for our generator
    unif01_Gen* gen = unif01_CreateExternGenBits(name, xorshift);

    // Run the tests.
    bbattery_SmallCrush(gen);

    // Clean up.
    unif01_DeleteExternGenBits(gen);

    return 0;
}
