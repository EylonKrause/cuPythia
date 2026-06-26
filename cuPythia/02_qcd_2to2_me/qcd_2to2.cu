// cuPythia kernel 02 — batched QCD 2->2 matrix element  g g -> g g.
//
// Verbatim port of Pythia 8.317 Sigma2gg2gg::sigmaKin (src/SigmaQCD.cc:115-129),
// the study's #1 GPU target: pure branchless double arithmetic on
// (sHat, tHat, uHat, alphaS), one CUDA thread per trial.
//
//   sigma = (pi/sH^2) * alpS^2 * 0.5 * (sigTS + sigUS + sigTU),  each = (9/4)(...)
//
// Validation (two independent checks):
//   (a) GPU vs the SAME formula on the CPU  -> proves the port is bit-correct.
//   (b) the formula vs the textbook analytic gg->gg differential cross section
//          dσ/dt̂ = (pi/sH^2) alpS^2 (9/2)(3 - tu/s^2 - su/t^2 - st/u^2)
//       -> proves the physics (independent expression, non-circular).
//
// Reports BOTH kernel-only speedup and transfer-inclusive speedup, because the
// honest end-to-end story is PCIe-bound, not kernel-bound (Amdahl).
//
// Build: nvcc -O3 -arch=sm_120 -o qcd_2to2 qcd_2to2.cu
// Run:   ./qcd_2to2 [nTrials=10000000]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <vector>
#include <cuda_runtime.h>
#include "../common/rng.cuh"

__host__ __device__ inline double pow2(double x) { return x * x; }

// --- verbatim Pythia 8.317 Sigma2gg2gg::sigmaKin arithmetic (SigmaQCD.cc:115-129)
__host__ __device__ inline double gg2gg_sigma(double sH, double tH, double uH,
                                              double alpS) {
  double sH2 = sH * sH, tH2 = tH * tH, uH2 = uH * uH;
  double sigTS = (9. / 4.) * (tH2 / sH2 + 2. * tH / sH + 3. + 2. * sH / tH + sH2 / tH2);
  double sigUS = (9. / 4.) * (uH2 / sH2 + 2. * uH / sH + 3. + 2. * sH / uH + sH2 / uH2);
  double sigTU = (9. / 4.) * (tH2 / uH2 + 2. * tH / uH + 3. + 2. * uH / tH + uH2 / tH2);
  double sigSum = sigTS + sigUS + sigTU;
  return (M_PI / sH2) * pow2(alpS) * 0.5 * sigSum;
}

// --- independent textbook cross-check (Ellis-Stirling-Webber / PDG).
// The bare amplitude gives (9/2)(3 - tu/s^2 - su/t^2 - st/u^2); Pythia's sigmaKin
// folds in the 1/2 identical-gluon factor (SigmaQCD.cc:126), so the matching
// differential cross section carries 9/4, not 9/2. (Verified by hand: Pythia's
// rearranged bracket sum B equals exactly 2*(3 - tu/s^2 - su/t^2 - st/u^2).)
__host__ inline double gg2gg_textbook(double sH, double tH, double uH, double alpS) {
  double s2 = sH * sH;
  double br = 3.0 - tH * uH / (sH * sH) - sH * uH / (tH * tH) - sH * tH / (uH * uH);
  return (M_PI / s2) * pow2(alpS) * (9.0 / 4.0) * br; // 9/4 includes identical-gluon 1/2
}

__global__ void meKernel(const double* sH, const double* tH, const double* uH,
                         const double* alpS, double* out, int n) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) out[i] = gg2gg_sigma(sH[i], tH[i], uH[i], alpS[i]);
}

#define CK(call) do { cudaError_t e_=(call); if(e_!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); return 1; } } while(0)
using Clock = std::chrono::high_resolution_clock;
static double ms(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

int main(int argc, char** argv) {
  int n = (argc > 1) ? atoi(argv[1]) : 10000000;
  double s = 100.0 * 100.0;   // sHat = (100 GeV)^2
  double alphaS = 0.118;

  std::vector<double> hS(n), hT(n), hU(n), hA(n), hRef(n), hTb(n), hGpu(n);
  uint64_t ctr = 0x1234567ULL;
  for (int i = 0; i < n; ++i) {
    // cos(theta) kept away from the collinear poles t,u -> 0
    double c = (2.0 * u01(splitmix64(ctr++)) - 1.0) * 0.98;
    double t = -0.5 * s * (1.0 - c);
    hS[i] = s; hT[i] = t; hU[i] = -s - t; hA[i] = alphaS; // massless: s+t+u=0
  }

  // CPU reference (verbatim Pythia formula), timed; plus textbook cross-check.
  auto c0 = Clock::now();
  for (int i = 0; i < n; ++i) hRef[i] = gg2gg_sigma(hS[i], hT[i], hU[i], hA[i]);
  auto c1 = Clock::now();
  double cpuMs = ms(c0, c1);
  for (int i = 0; i < n; ++i) hTb[i] = gg2gg_textbook(hS[i], hT[i], hU[i], hA[i]);

  double *dS, *dT, *dU, *dA, *dO;
  size_t bytes = (size_t)n * sizeof(double);
  CK(cudaMalloc(&dS, bytes)); CK(cudaMalloc(&dT, bytes)); CK(cudaMalloc(&dU, bytes));
  CK(cudaMalloc(&dA, bytes)); CK(cudaMalloc(&dO, bytes));
  int threads = 256, blocks = (n + threads - 1) / threads;
  meKernel<<<blocks, threads>>>(dS, dT, dU, dA, dO, n); // warmup
  CK(cudaDeviceSynchronize());

  auto h0 = Clock::now();
  CK(cudaMemcpy(dS, hS.data(), bytes, cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dT, hT.data(), bytes, cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dU, hU.data(), bytes, cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dA, hA.data(), bytes, cudaMemcpyHostToDevice));
  auto h1 = Clock::now();
  meKernel<<<blocks, threads>>>(dS, dT, dU, dA, dO, n);
  CK(cudaDeviceSynchronize());
  auto h2 = Clock::now();
  CK(cudaMemcpy(hGpu.data(), dO, bytes, cudaMemcpyDeviceToHost));
  auto h3 = Clock::now();

  double h2dMs = ms(h0, h1), kernMs = ms(h1, h2), d2hMs = ms(h2, h3);

  double maxRelGpu = 0.0, maxRelTb = 0.0;
  for (int i = 0; i < n; ++i) {
    double r = hRef[i];
    maxRelGpu = fmax(maxRelGpu, fabs(hGpu[i] - r) / fabs(r));
    maxRelTb  = fmax(maxRelTb,  fabs(hTb[i]  - r) / fabs(r));
  }

  printf("QCD g g -> g g  (Sigma2gg2gg::sigmaKin, verbatim Pythia 8.317)\n");
  printf("  trials            = %d   (sHat=(100 GeV)^2, alphaS=0.118)\n", n);
  printf("  max relerr GPU vs CPU-Pythia formula   = %.3e   (port correctness)\n", maxRelGpu);
  printf("  max relerr formula vs textbook analytic= %.3e   (physics check)\n", maxRelTb);
  printf("  CPU loop            = %8.2f ms\n", cpuMs);
  printf("  GPU kernel only     = %8.2f ms   -> kernel speedup     %.1fx\n",
         kernMs, cpuMs / kernMs);
  printf("  GPU incl. transfer  = %8.2f ms   (H2D %.1f + kern %.2f + D2H %.1f)\n",
         h2dMs + kernMs + d2hMs, h2dMs, kernMs, d2hMs);
  printf("                                       -> end-to-end speedup %.1fx (PCIe-bound)\n",
         cpuMs / (h2dMs + kernMs + d2hMs));
  bool ok = (maxRelGpu < 1e-12) && (maxRelTb < 1e-10);
  printf("VALIDATION: %s (GPU==CPU to <1e-12 AND formula==textbook to <1e-10)\n",
         ok ? "PASS" : "FAIL");

  cudaFree(dS); cudaFree(dT); cudaFree(dU); cudaFree(dA); cudaFree(dO);
  return ok ? 0 : 2;
}
