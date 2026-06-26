// cuPythia toolchain check — GPU Monte Carlo estimate of pi.
//
// Purpose: prove the CUDA toolchain + RTX 5050 work end-to-end, and establish
// the host/device-identical counter-based RNG that later cuPythia kernels reuse.
// Known answer (pi) makes correctness unambiguous; throughput is compared
// fairly as samples/second so GPU and CPU need not run the same sample count.
//
// Build: nvcc -O3 -arch=native -o mc_pi mc_pi.cu
// Run:   ./mc_pi [samplesPerThread]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <cuda_runtime.h>
#include "../common/rng.cuh"

__global__ void mcPiKernel(uint64_t seed, uint64_t samplesPerThread,
                           unsigned long long* globalInside) {
  uint64_t tid = blockIdx.x * (uint64_t)blockDim.x + threadIdx.x;
  uint64_t ctr = seed + tid * 0x100000001B3ULL; // distinct stream per thread
  unsigned long long inside = 0;
  for (uint64_t i = 0; i < samplesPerThread; ++i) {
    double x = u01(splitmix64(ctr++));
    double y = u01(splitmix64(ctr++));
    if (x * x + y * y <= 1.0) ++inside;
  }
  atomicAdd(globalInside, inside);
}

static double cpuPi(uint64_t seed, uint64_t total) {
  uint64_t ctr = seed;
  uint64_t inside = 0;
  for (uint64_t i = 0; i < total; ++i) {
    double x = u01(splitmix64(ctr++));
    double y = u01(splitmix64(ctr++));
    if (x * x + y * y <= 1.0) ++inside;
  }
  return 4.0 * (double)inside / (double)total;
}

#define CK(call) do { cudaError_t e_=(call); if(e_!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); return 1; } } while(0)

int main(int argc, char** argv) {
  int dev = 0; cudaDeviceProp prop;
  CK(cudaGetDevice(&dev));
  CK(cudaGetDeviceProperties(&prop, dev));
  printf("Device: %s  (SM %d.%d, %d SMs)\n", prop.name, prop.major, prop.minor,
         prop.multiProcessorCount);

  const int blocks = 1024, threads = 256;
  uint64_t perThread = (argc > 1) ? strtoull(argv[1], nullptr, 10) : 4000ULL;
  uint64_t totalGpu = (uint64_t)blocks * threads * perThread;
  uint64_t seed = 0xC0FFEEULL;

  unsigned long long* dInside;
  CK(cudaMalloc(&dInside, sizeof(unsigned long long)));

  // Warm up (kernel JIT/load) so timing reflects steady state.
  CK(cudaMemset(dInside, 0, sizeof(unsigned long long)));
  mcPiKernel<<<blocks, threads>>>(seed, 1, dInside);
  CK(cudaDeviceSynchronize());

  CK(cudaMemset(dInside, 0, sizeof(unsigned long long)));
  auto t0 = std::chrono::high_resolution_clock::now();
  mcPiKernel<<<blocks, threads>>>(seed, perThread, dInside);
  CK(cudaDeviceSynchronize());
  auto t1 = std::chrono::high_resolution_clock::now();

  unsigned long long hInside = 0;
  CK(cudaMemcpy(&hInside, dInside, sizeof(unsigned long long), cudaMemcpyDeviceToHost));
  double piGpu = 4.0 * (double)hInside / (double)totalGpu;
  double gpuMs = std::chrono::duration<double, std::milli>(t1 - t0).count();

  // CPU baseline (cap so it stays a few hundred ms; rates are the fair metric).
  uint64_t totalCpu = totalGpu < 200000000ULL ? totalGpu : 200000000ULL;
  auto c0 = std::chrono::high_resolution_clock::now();
  double piCpu = cpuPi(seed, totalCpu);
  auto c1 = std::chrono::high_resolution_clock::now();
  double cpuMs = std::chrono::duration<double, std::milli>(c1 - c0).count();

  double gpuRate = totalGpu / (gpuMs / 1000.0);
  double cpuRate = totalCpu / (cpuMs / 1000.0);

  printf("GPU: pi=%.6f  err=%.2e  samples=%.3e  time=%.2f ms  rate=%.3e /s\n",
         piGpu, fabs(piGpu - M_PI), (double)totalGpu, gpuMs, gpuRate);
  printf("CPU: pi=%.6f  err=%.2e  samples=%.3e  time=%.2f ms  rate=%.3e /s\n",
         piCpu, fabs(piCpu - M_PI), (double)totalCpu, cpuMs, cpuRate);
  printf("Throughput speedup (samples/s, GPU/CPU, 1 CPU thread): %.1fx\n", gpuRate / cpuRate);

  bool ok = fabs(piGpu - M_PI) < 1e-3;
  printf("VALIDATION: %s (|pi_gpu - pi| < 1e-3)\n", ok ? "PASS" : "FAIL");
  cudaFree(dInside);
  return ok ? 0 : 2;
}
