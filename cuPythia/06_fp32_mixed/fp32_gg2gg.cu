// cuPythia kernel 06 — FP32 vs FP64: breaking the consumer-GPU FP64 ceiling.
//
// Kernel 03 (fused resident gg->gg) topped out ~7-9x, FP64-division-bound:
// consumer Blackwell (RTX 5050) runs FP64 at ~1/64 of FP32. This runs the SAME
// fused MC in BOTH double and float (identical RNG samples), so the comparison
// is apples-to-apples. Because the Monte-Carlo statistical error (~1/sqrt(N))
// dwarfs FP32 rounding, float gives the SAME physics accuracy at much higher
// throughput -- the honest way to break the ceiling.
//
// Build: nvcc -O3 -arch=sm_120 -o fp32_gg2gg fp32_gg2gg.cu
// Run:   ./fp32_gg2gg [trialsPerThread=4000]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <cuda_runtime.h>
#include "../common/rng.cuh"

template <typename T>
__host__ __device__ inline T gg2gg_sigmaT(T sH, T tH, T uH, T alpS) {
  T sH2 = sH * sH, tH2 = tH * tH, uH2 = uH * uH;
  T a = (T)(9.0 / 4.0) * (tH2 / sH2 + (T)2 * tH / sH + (T)3 + (T)2 * sH / tH + sH2 / tH2);
  T b = (T)(9.0 / 4.0) * (uH2 / sH2 + (T)2 * uH / sH + (T)3 + (T)2 * sH / uH + sH2 / uH2);
  T c = (T)(9.0 / 4.0) * (tH2 / uH2 + (T)2 * tH / uH + (T)3 + (T)2 * uH / tH + uH2 / tH2);
  return ((T)M_PI / sH2) * alpS * alpS * (T)0.5 * (a + b + c);
}

template <typename T>
__global__ void fusedT(uint64_t seed, uint64_t kPer, T s, T alpS, T cMax, T* gSum) {
  uint64_t tid = blockIdx.x * (uint64_t)blockDim.x + threadIdx.x;
  uint64_t ctr = seed + tid * 0x100000001B3ULL;
  T local = (T)0;
  for (uint64_t i = 0; i < kPer; ++i) {
    T c = (T)(2.0 * u01(splitmix64(ctr++)) - 1.0) * cMax; // same RNG draw, cast to T
    T t = (T)(-0.5) * s * ((T)1 - c);
    local += gg2gg_sigmaT<T>(s, t, -s - t, alpS);
  }
  atomicAdd(gSum, local);
}

static double simpson(int N, double s, double alpS, double cMax) {
  double a = -cMax, h = (2.0 * cMax) / N, sum = 0.0;
  for (int i = 0; i <= N; ++i) {
    double cc = a + i * h, t = -0.5 * s * (1.0 - cc);
    double w = (i == 0 || i == N) ? 1.0 : (i % 2 ? 4.0 : 2.0);
    sum += w * gg2gg_sigmaT<double>(s, t, -s - t, alpS);
  }
  return (sum * h / 3.0) * (s / 2.0);
}

using Clock = std::chrono::steady_clock;
static double msec(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}
#define CK(call) do { cudaError_t e_=(call); if(e_!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); return 1; } } while(0)

template <typename T>
static int runPrec(uint64_t seed, uint64_t kPer, int blocks, int threads,
                   double s, double alpS, double cMax, double* sigma_out, double* ms_out) {
  T* d; CK(cudaMalloc(&d, sizeof(T)));
  T zero = (T)0; CK(cudaMemcpy(d, &zero, sizeof(T), cudaMemcpyHostToDevice));
  fusedT<T><<<blocks, threads>>>(seed, 1, (T)s, (T)alpS, (T)cMax, d); // warmup
  CK(cudaDeviceSynchronize());
  CK(cudaMemcpy(d, &zero, sizeof(T), cudaMemcpyHostToDevice));
  auto t0 = Clock::now();
  fusedT<T><<<blocks, threads>>>(seed, kPer, (T)s, (T)alpS, (T)cMax, d);
  CK(cudaDeviceSynchronize());
  auto t1 = Clock::now();
  T h = (T)0; CK(cudaMemcpy(&h, d, sizeof(T), cudaMemcpyDeviceToHost)); cudaFree(d);
  uint64_t total = (uint64_t)blocks * threads * kPer;
  *sigma_out = (2.0 * cMax) * (s / 2.0) * ((double)h / (double)total) * 0.3893793721;
  *ms_out = msec(t0, t1);
  return 0;
}

int main(int argc, char** argv) {
  uint64_t kPer = (argc > 1) ? strtoull(argv[1], nullptr, 10) : 4000ULL;
  const int blocks = 1024, threads = 256;
  double s = 100.0 * 100.0, alphaS = 0.118, cMax = 0.9, seedS = 0;
  uint64_t seed = 0xABCDEF01ULL;

  double sigD, msD, sigF, msF;
  if (runPrec<double>(seed, kPer, blocks, threads, s, alphaS, cMax, &sigD, &msD)) return 1;
  if (runPrec<float >(seed, kPer, blocks, threads, s, alphaS, cMax, &sigF, &msF)) return 1;
  double ref = simpson(2000000, s, alphaS, cMax) * 0.3893793721;
  uint64_t total = (uint64_t)blocks * threads * kPer;
  double relD = fabs(sigD - ref) / ref, relF = fabs(sigF - ref) / ref;
  (void)seedS;

  printf("FP32 vs FP64 fused g g -> g g  (%.3e trials, sHat=(100 GeV)^2)\n", (double)total);
  printf("  Simpson reference = %.6e mb\n", ref);
  printf("  FP64: sigma=%.6e mb  relerr=%.2e  %.1f ms  rate=%.3e/s\n",
         sigD, relD, msD, total / (msD / 1000.0));
  printf("  FP32: sigma=%.6e mb  relerr=%.2e  %.1f ms  rate=%.3e/s\n",
         sigF, relF, msF, total / (msF / 1000.0));
  printf("  FP32/FP64 speedup = %.1fx   (same MC accuracy: stat error >> FP32 rounding)\n",
         msD / msF);
  bool ok = (relD < 1e-3) && (relF < 5e-3);
  printf("VALIDATION: %s (both match quadrature; FP32 within MC tolerance)\n", ok ? "PASS" : "FAIL");
  return ok ? 0 : 2;
}
