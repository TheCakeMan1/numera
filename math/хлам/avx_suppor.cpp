#include "avx_support.h"

bool is_avx_supported() {
    // CPUID.01H:ECX.AVX[bit 28] = 1
    int cpuInfo[4];
    __cpuid(cpuInfo, 1);

    bool avxSupported = (cpuInfo[2] & (1 << 28)) != 0;
    return avxSupported;
}