#pragma once
#include <cstdint>
#include <cstdlib>

// Multi-GPU sharding: each GPU runs the SAME binary on a disjoint CONTIGUOUS slice of the global event
// index space. The launcher (run.sh) sets CUPYTHIA_SHARD (this GPU's index i in 0..G-1) and
// CUPYTHIA_SHARD_N (events per shard, M). cupythia_shard_offset() returns i*M; the kernel ADDS it to
// its local event index e before deriving ANY seed (eg = e + offset), so shard i computes exactly
// global events [i*M, i*M+cnt). Since every seed of a given event is a pure function of that global
// index, the merged output is BIT-IDENTICAL to one GPU running all events -- same seeds, same events,
// same refragment-drops, same order (structurally disjoint, no overlap, no gap, no probabilistic
// collision). With the env vars absent (default / single-GPU) the offset is 0, so eg==e and committed
// results stay byte-identical. Host-only (called from main(); pass the result into the kernel).
inline uint64_t cupythia_shard_offset(){
  const char* si = getenv("CUPYTHIA_SHARD");
  const char* sn = getenv("CUPYTHIA_SHARD_N");
  if(!si || !sn) return 0;                                     // not sharded -> no offset (e unchanged)
  return (uint64_t)strtoull(si,nullptr,10) * (uint64_t)strtoull(sn,nullptr,10);
}

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
