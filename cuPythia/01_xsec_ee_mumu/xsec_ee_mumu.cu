// cuPythia kernel 01 — GPU Monte-Carlo integration of sigma(e+ e- -> mu+ mu-).
//
// Tree-level QED (s-channel photon, massless muons, no Z):
//     dsigma/dOmega = alpha^2 / (4 s) * (1 + cos^2 theta)
//     sigma_total   = 4 pi alpha^2 / (3 s)            <-- analytic check
//
// MC: sample cos(theta) uniformly in [-1,1] (uniform in solid angle together
// with phi, total measure 4 pi), average the integrand, multiply by 4 pi.
// Validated against the closed form above; speedup measured vs 1 CPU thread.
// Reuses the host/device SplitMix64 RNG from the toolchain check.
//
// Build: nvcc -O3 -arch=sm_120 -o xsec xsec_ee_mumu.cu
// Run:   ./xsec [sqrt_s_GeV=10] [samplesPerThread=20000]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <cuda_runtime.h>
#include "../common/rng.cuh"

// alpha^2/(4s) * (1+cos^2) integrated over 4pi solid angle, MC-averaged.
__global__ void xsecKernel(uint64_t seed, uint64_t nPerThread,
                           double pref, double* gSum) {
  uint64_t tid = blockIdx.x * (uint64_t)blockDim.x + threadIdx.x;
  uint64_t ctr = seed + tid * 0x100000001B3ULL;
  double local = 0.0;
  for (uint64_t i = 0; i < nPerThread; ++i) {
    double c = 2.0 * u01(splitmix64(ctr++)) - 1.0; // cos theta in [-1,1]
    local += pref * (1.0 + c * c);                 // dsigma/dOmega
  }
  atomicAdd(gSum, local);
}

static double cpuSum(uint64_t seed, uint64_t total, double pref) {
  uint64_t ctr = seed;
  double sum = 0.0;
  for (uint64_t i = 0; i < total; ++i) {
    double c = 2.0 * u01(splitmix64(ctr++)) - 1.0;
    sum += pref * (1.0 + c * c);
  }
  return sum;
}

#define CK(call) do { cudaError_t e_=(call); if(e_!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); return 1; } } while(0)

int main(int argc, char** argv) {
  const double alpha   = 1.0 / 137.035999084;
  const double GeVm2_nb = 0.3893793721e6; // 1 GeV^-2 = 0.3894 mb = 3.8938e5 nb
  double sqrtS = (argc > 1) ? atof(argv[1]) : 10.0;     // GeV
  uint64_t perThread = (argc > 2) ? strtoull(argv[2], nullptr, 10) : 20000ULL;
  double s = sqrtS * sqrtS;
  double pref = alpha * alpha / (4.0 * s);              // GeV^-2 per steradian
  double sigmaAnalytic_nb = (4.0 * M_PI * alpha * alpha / (3.0 * s)) * GeVm2_nb;

  const int blocks = 1024, threads = 256;
  uint64_t totalGpu = (uint64_t)blocks * threads * perThread;
  uint64_t seed = 0xBEEF1234ULL;

  double* dSum; CK(cudaMalloc(&dSum, sizeof(double)));
  CK(cudaMemset(dSum, 0, sizeof(double)));
  xsecKernel<<<blocks, threads>>>(seed, 1, pref, dSum); // warmup
  CK(cudaDeviceSynchronize());

  CK(cudaMemset(dSum, 0, sizeof(double)));
  auto t0 = std::chrono::high_resolution_clock::now();
  xsecKernel<<<blocks, threads>>>(seed, perThread, pref, dSum);
  CK(cudaDeviceSynchronize());
  auto t1 = std::chrono::high_resolution_clock::now();

  double hSum = 0.0; CK(cudaMemcpy(&hSum, dSum, sizeof(double), cudaMemcpyDeviceToHost));
  double sigmaGpu_nb = 4.0 * M_PI * (hSum / (double)totalGpu) * GeVm2_nb;
  double gpuMs = std::chrono::duration<double, std::milli>(t1 - t0).count();

  uint64_t totalCpu = totalGpu < 200000000ULL ? totalGpu : 200000000ULL;
  auto c0 = std::chrono::high_resolution_clock::now();
  double cSum = cpuSum(seed, totalCpu, pref);
  auto c1 = std::chrono::high_resolution_clock::now();
  double sigmaCpu_nb = 4.0 * M_PI * (cSum / (double)totalCpu) * GeVm2_nb;
  double cpuMs = std::chrono::duration<double, std::milli>(c1 - c0).count();

  double relErr = fabs(sigmaGpu_nb - sigmaAnalytic_nb) / sigmaAnalytic_nb;
  double gpuRate = totalGpu / (gpuMs / 1000.0);
  double cpuRate = totalCpu / (cpuMs / 1000.0);

  printf("e+e- -> mu+mu-  at sqrt(s) = %.1f GeV\n", sqrtS);
  printf("  analytic sigma = %.6f nb   (4 pi alpha^2 / 3s)\n", sigmaAnalytic_nb);
  printf("  GPU MC   sigma = %.6f nb   relerr=%.2e  (%.3e samples, %.2f ms, %.3e/s)\n",
         sigmaGpu_nb, relErr, (double)totalGpu, gpuMs, gpuRate);
  printf("  CPU MC   sigma = %.6f nb   (%.3e samples, %.2f ms, %.3e/s)\n",
         sigmaCpu_nb, (double)totalCpu, cpuMs, cpuRate);
  printf("  throughput speedup (GPU/CPU, 1 thread): %.1fx\n", gpuRate / cpuRate);
  bool ok = relErr < 1e-3;
  printf("VALIDATION: %s (relerr < 1e-3 vs analytic)\n", ok ? "PASS" : "FAIL");
  cudaFree(dSum);
  return ok ? 0 : 2;
}
