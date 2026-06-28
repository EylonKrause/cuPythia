#pragma once
#include <cstdint>

// Portability: MSVC's <cmath> does not define M_PI unless _USE_MATH_DEFINES was
// set before <cmath> was included. Provide a fallback so every kernel builds
// identically on Windows (MSVC) and Linux (glibc). All kernels include this header.
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Portability: hardware double-precision atomicAdd is native on Pascal (sm_60) and up — which is
// cuPythia's minimum target. This CAS fallback is compiled ONLY for older arches (Maxwell sm_50,
// Kepler sm_3x), so the generator also builds for pre-Pascal GPUs without touching the sm_60+ path
// (where the native instruction is used). No effect on host or on sm_60+ device code.
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
__device__ inline double atomicAdd(double* address, double val) {
  unsigned long long int* a = (unsigned long long int*)address;
  unsigned long long int old = *a, assumed;
  do { assumed = old;
    old = atomicCAS(a, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
  } while (assumed != old);              // NaN-aware: loops until the CAS sticks
  return __longlong_as_double(old);
}
#endif

// Host/device-identical SplitMix64 counter-based RNG (no cuRAND dependency).
// Counter-based so each GPU thread derives an independent stream from (seed,tid)
// with no per-thread state to store, and CPU/GPU produce bit-identical draws —
// which is what lets cuPythia validate GPU kernels against a CPU reference exactly.
__host__ __device__ inline uint64_t splitmix64(uint64_t x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}

// uint64 -> double in [0,1), using the top 53 bits (one ULP = 2^-53).
__host__ __device__ inline double u01(uint64_t u) {
  return (u >> 11) * (1.0 / 9007199254740992.0); // 1 / 2^53
}
