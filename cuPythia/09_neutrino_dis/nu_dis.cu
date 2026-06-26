// cuPythia kernel 09 — neutrino-quark deep-inelastic scattering (parton model).
//
// HONEST FRAMING: neutrino-NUCLEUS interactions (nuclear effects, Fermi motion,
// final-state interactions) are the domain of SPECIALIZED generators — GENIE,
// NuWro, GiBUU — NOT Pythia. Pythia's role in neutrino physics is the DIS
// final-state HADRONIZATION. cuPythia cannot and does not add nuclear physics.
//
// What it CAN add honestly is the PARTONIC neutrino DIS cross section — the
// textbook electroweak result that neutrino experiments rest on, and the famous
// signature that revealed valence quarks:
//     dσ/dy(ν q)    ∝ 1          (flat in inelasticity y)
//     dσ/dy(ν qbar) ∝ (1 - y)^2
//     dσ/dy = (G_F^2 s / π) * shape(y)
// MC-integrated on GPU and validated against the analytic integrals
//     ∫₀¹ 1 dy = 1,   ∫₀¹ (1-y)^2 dy = 1/3,   ratio = 3.
// (The x-dependence needs real PDFs / LHAPDF; only the analytic y-structure is here.)
//
// Build: nvcc -O3 -arch=sm_120 -o nu_dis nu_dis.cu
// Run:   ./nu_dis [trialsPerThread=4000]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <cuda_runtime.h>
#include "../common/rng.cuh"

// MC-integrate shape(y) over y in [0,1]; antiquark=1 -> (1-y)^2, else flat.
__global__ void nuKernel(uint64_t seed, uint64_t nPer, int antiquark, double* gSum) {
  uint64_t tid = blockIdx.x * (uint64_t)blockDim.x + threadIdx.x;
  uint64_t ctr = seed + tid * 0x100000001B3ULL;
  double local = 0.0;
  for (uint64_t i = 0; i < nPer; ++i) {
    double y = u01(splitmix64(ctr++));
    local += antiquark ? (1.0 - y) * (1.0 - y) : 1.0;
  }
  atomicAdd(gSum, local);
}

#define CK(call) do { cudaError_t e_=(call); if(e_!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); return 1; } } while(0)

static int runFlavour(uint64_t seed, uint64_t nPer, int blocks, int threads,
                      int antiquark, double* meanShape) {
  double* d; CK(cudaMalloc(&d, sizeof(double))); CK(cudaMemset(d, 0, sizeof(double)));
  nuKernel<<<blocks, threads>>>(seed, nPer, antiquark, d);
  CK(cudaDeviceSynchronize());
  double h = 0.0; CK(cudaMemcpy(&h, d, sizeof(double), cudaMemcpyDeviceToHost)); cudaFree(d);
  *meanShape = h / ((double)blocks * threads * nPer);   // ∫₀¹ shape dy
  return 0;
}

int main(int argc, char** argv) {
  uint64_t nPer = (argc > 1) ? strtoull(argv[1], nullptr, 10) : 4000ULL;
  const int blocks = 1024, threads = 256;
  uint64_t seed = 0x4E55ULL; // 'NU'
  double GF = 1.1663787e-5;                 // GeV^-2
  double Enu = 100.0, mN = 0.938;           // GeV: 100 GeV neutrino on a nucleon
  double s = 2.0 * mN * Enu;                // partonic s (x=1 reference)
  double GeVm2_cm2 = 0.3893793721e-27;      // 1 GeV^-2 = 0.3894 mb = 0.3894e-27 cm^2

  double mq, mqbar;
  if (runFlavour(seed, nPer, blocks, threads, 0, &mq))    return 1; // ν q
  if (runFlavour(seed, nPer, blocks, threads, 1, &mqbar)) return 1; // ν qbar
  double pre = GF * GF * s / M_PI;          // dσ/dy normalisation
  double sig_q    = pre * mq    * GeVm2_cm2;
  double sig_qbar = pre * mqbar * GeVm2_cm2;
  double ratio = mq / mqbar;

  printf("Neutrino-quark DIS (parton model), E_nu=%.0f GeV, partonic s=%.1f GeV^2\n", Enu, s);
  printf("  <shape> nu q     = %.5f   (analytic 1)      sigma_q    = %.3e cm^2\n", mq, sig_q);
  printf("  <shape> nu qbar  = %.5f   (analytic 1/3)    sigma_qbar = %.3e cm^2\n", mqbar, sig_qbar);
  printf("  sigma(nu q)/sigma(nu qbar) = %.4f   (analytic 3 -> valence-quark signature)\n", ratio);
  bool ok = fabs(mq - 1.0) < 2e-3 && fabs(mqbar - 1.0/3.0) < 2e-3 && fabs(ratio - 3.0) < 2e-2;
  printf("VALIDATION: %s (y-distributions match flat & (1-y)^2; ratio = 3)\n", ok ? "PASS" : "FAIL");
  return ok ? 0 : 2;
}
