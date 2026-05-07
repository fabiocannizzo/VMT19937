/* CPUID probe: prints detected ISA to stdout. Used by Makefile for MSVC native detection. */
#include <stdio.h>
#ifdef _MSC_VER
#  include <intrin.h>
static void do_cpuid(unsigned info[4], unsigned leaf, unsigned subleaf)
{ __cpuidex((int*)info, (int)leaf, (int)subleaf); }
#else
#  include <cpuid.h>
static void do_cpuid(unsigned info[4], unsigned leaf, unsigned subleaf)
{ __cpuid_count(leaf, subleaf, info[0], info[1], info[2], info[3]); }
#endif

int main(void)
{
    unsigned info[4];
    do_cpuid(info, 7, 0);
    if (info[1] & (1u << 16)) {      /* AVX-512F */
        puts("avx512vl");            /* all current AVX-512F CPUs also have VL */
    } else if (info[1] & (1u << 5)) { /* AVX2 */
        puts("avx2");
    } else {
        do_cpuid(info, 1, 0);
        if (info[2] & (1u << 20))   /* SSE4.2 */
            puts("sse42");
        else
            puts("scalar");
    }
    return 0;
}
