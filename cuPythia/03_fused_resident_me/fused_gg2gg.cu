// cuPythia kernel 03 — FUSED, GPU-RESIDENT phase-space + ME for g g -> g g.
//
// The point of the project, demonstrated. Kernel 02 transferred pre-generated
// (sHat,tHat,uHat) arrays and got only 1.3x end-to-end (PCIe + FP64 bound).
// Here each thread instead GENERATES its own trials on-device (counter-based
// RNG, many trials per thread, all in registers), evaluates the SAME verbatim
// Pythia gg->gg matrix element, and accumulates a Monte-Carlo integral. Nothing
// crosses PCIe except one scalar. That is how you keep the speedup.
//
// Physics: MC-integrate dσ/dt̂ over a cosθ cut [-cMax, cMax] (avoids the t,u->0
// Rutherford poles). Validate against a deterministic CPU Simpson quadrature of
// the identical integrand (independent of the RNG) + a CPU MC at the same RNG.
//
// Build: nvcc -O3 -arch=sm_120 -o fused_gg2gg fused_gg2gg.cu
// Run:   ./fused_gg2gg [trialsPerThread=4000]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <cuda_runtime.h>
#include "../common/rng.cuh"

__host__ __device__ inline double pow2(double x) { return x * x; }

// verbatim Pythia 8.317 Sigma2gg2gg::sigmaKin (src/SigmaQCD.cc:115-129) = dσ/dt̂
__host__ __device__ inline double gg2gg_sigma(double sH, double tH, double uH, double alpS) {
  double sH2 = sH * sH, tH2 = tH * tH, uH2 = uH * uH;
  double sigTS = (9. / 4.) * (tH2 / sH2 + 2. * tH / sH + 3. + 2. * sH / tH + sH2 / tH2);
  double sigUS = (9. / 4.) * (uH2 / sH2 + 2. * uH / sH + 3. + 2. * sH / uH + sH2 / uH2);
  double sigTU = (9. / 4.) * (tH2 / uH2 + 2. * tH / uH + 3. + 2. * uH / tH + uH2 / tH2);
  return (M_PI / sH2) * pow2(alpS) * 0.5 * (sigTS + sigUS + sigTU);
}

// each thread: kPer trials, generate cosθ -> (t,u) -> ME, accumulate. Resident.
__global__ void fusedKernel(uint64_t seed, uint64_t kPer, double s, double alpS,
                            double cMax, double* gSum) {
  uint64_t tid = blockIdx.x * (uint64_t)blockDim.x + threadIdx.x;
  uint64_t ctr = seed + tid * 0x100000001B3ULL;
  double local = 0.0;
  for (uint64_t i = 0; i < kPer; ++i) {
    double c = (2.0 * u01(splitmix64(ctr++)) - 1.0) * cMax; // cosθ in [-cMax,cMax]
    double t = -0.5 * s * (1.0 - c);
    local += gg2gg_sigma(s, t, -s - t, alpS);               // dσ/dt̂
  }
  atomicAdd(gSum, local);
}

static double cpuMC(uint64_t seed, uint64_t total, double s, double alpS, double cMax) {
  uint64_t ctr = seed; double sum = 0.0;
  for (uint64_t i = 0; i < total; ++i) {
    double c = (2.0 * u01(splitmix64(ctr++)) - 1.0) * cMax;
    double t = -0.5 * s * (1.0 - c);
    sum += gg2gg_sigma(s, t, -s - t, alpS);
  }
  // <f> over cosθ width (2 cMax), times Jacobian dt̂/dcosθ = s/2
  return (2.0 * cMax) * (s / 2.0) * (sum / (double)total);
}

// deterministic reference: Simpson quadrature of dσ/dt̂ over the same cut.
static double simpson(int N, double s, double alpS, double cMax) {
  double a = -cMax, b = cMax, h = (b - a) / N, sum = 0.0;
  for (int i = 0; i <= N; ++i) {
    double c = a + i * h, t = -0.5 * s * (1.0 - c);
    double f = gg2gg_sigma(s, t, -s - t, alpS);
    double w = (i == 0 || i == N) ? 1.0 : (i % 2 ? 4.0 : 2.0);
    sum += w * f;
  }
  return (sum * h / 3.0) * (s / 2.0); // include Jacobian dt̂/dcosθ = s/2
}

#define CK(call) do { cudaError_t e_=(call); if(e_!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); return 1; } } while(0)
using Clock = std::chrono::steady_clock;
static double msec(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

int main(int argc, char** argv) {
  uint64_t kPer = (argc > 1) ? strtoull(argv[1], nullptr, 10) : 4000ULL;
  const int blocks = 1024, threads = 256;
  double s = 100.0 * 100.0, alphaS = 0.118, cMax = 0.9;
  double GeVm2_mb = 0.3893793721; // 1 GeV^-2 = 0.3894 mb
  uint64_t seed = 0xABCDEF01ULL;
  uint64_t totalGpu = (uint64_t)blocks * threads * kPer;

  double* dSum; CK(cudaMalloc(&dSum, sizeof(double)));
  CK(cudaMemset(dSum, 0, sizeof(double)));
  fusedKernel<<<blocks, threads>>>(seed, 1, s, alphaS, cMax, dSum); // warmup
  CK(cudaDeviceSynchronize());

  CK(cudaMemset(dSum, 0, sizeof(double)));
  auto g0 = Clock::now();
  fusedKernel<<<blocks, threads>>>(seed, kPer, s, alphaS, cMax, dSum);
  CK(cudaDeviceSynchronize());
  auto g1 = Clock::now();
  double hSum = 0.0; CK(cudaMemcpy(&hSum, dSum, sizeof(double), cudaMemcpyDeviceToHost));
  double sigGpu = (2.0 * cMax) * (s / 2.0) * (hSum / (double)totalGpu) * GeVm2_mb;
  double gpuMs = msec(g0, g1);

  uint64_t totalCpu = totalGpu < 200000000ULL ? totalGpu : 200000000ULL;
  auto c0 = Clock::now();
  double sigCpu = cpuMC(seed, totalCpu, s, alphaS, cMax) * GeVm2_mb;
  auto c1 = Clock::now();
  double cpuMs = msec(c0, c1);

  double sigRef = simpson(2000000, s, alphaS, cMax) * GeVm2_mb; // deterministic
  double relGpu = fabs(sigGpu - sigRef) / sigRef;
  double gpuRate = totalGpu / (gpuMs / 1000.0), cpuRate = totalCpu / (cpuMs / 1000.0);

  printf("FUSED resident g g -> g g  (sHat=(100 GeV)^2, alphaS=0.118, |cosθ|<%.1f)\n", cMax);
  printf("  Simpson quadrature  sigma = %.6e mb   (deterministic reference)\n", sigRef);
  printf("  GPU fused MC        sigma = %.6e mb   relerr=%.2e  (%.3e trials, %.2f ms)\n",
         sigGpu, relGpu, (double)totalGpu, gpuMs);
  printf("  CPU MC (same RNG)   sigma = %.6e mb   (%.3e trials, %.2f ms)\n",
         sigCpu, (double)totalCpu, cpuMs);
  printf("  GPU rate = %.3e/s   CPU rate = %.3e/s   -> speedup %.1fx\n",
         gpuRate, cpuRate, gpuRate / cpuRate);
  bool ok = relGpu < 2e-3;
  printf("VALIDATION: %s (GPU MC vs deterministic quadrature, relerr < 2e-3)\n",
         ok ? "PASS" : "FAIL");
  cudaFree(dSum);
  return ok ? 0 : 2;
}
