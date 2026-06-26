// cuPythia kernel 05 — counter-based per-event REPRODUCIBLE RNG.
//
// A capability stock Pythia lacks. Pythia's RNG (RANMAR, a sequential lagged-
// Fibonacci state machine) means event N is only reachable by advancing the
// generator N times or by serializing/restoring state. With a COUNTER-BASED RNG
// each event's randomness is a pure function of (globalSeed, eventId): you can
// regenerate ANY single event independently, on ANY node, in O(1) — no replay,
// no checkpoint. That is exactly what large-scale GRID/distributed production and
// debugging want (re-run only the one failed event of a million-event job;
// bit-identical results across heterogeneous nodes).
//
// Each "event" here gets an independent substream seeded by (globalSeed, eventId)
// and draws a VARIABLE number of randoms (mimicking variable event multiplicity),
// producing one observable. We then prove:
//   (1) regenerating a SHUFFLED subset of events out of order is BIT-IDENTICAL;
//   (2) partitioning events across S "nodes" reproduces the full ensemble exactly.
//
// Build: nvcc -O3 -arch=sm_120 -o repro_events repro_events.cu
// Run:   ./repro_events [nEvents=2000000]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include "../common/rng.cuh"

__host__ __device__ inline double pow2(double x) { return x * x; }
__host__ __device__ inline double gg2gg_sigma(double sH, double tH, double uH, double alpS) {
  double sH2 = sH * sH, tH2 = tH * tH, uH2 = uH * uH;
  double a = (9. / 4.) * (tH2 / sH2 + 2. * tH / sH + 3. + 2. * sH / tH + sH2 / tH2);
  double b = (9. / 4.) * (uH2 / sH2 + 2. * uH / sH + 3. + 2. * sH / uH + sH2 / uH2);
  double c = (9. / 4.) * (tH2 / uH2 + 2. * tH / uH + 3. + 2. * uH / tH + uH2 / tH2);
  return (M_PI / sH2) * pow2(alpS) * 0.5 * (a + b + c);
}

// The whole point: observable(seed, e) depends ONLY on (seed, e) — no sequential
// state, no dependence on which other events ran, or in what order, or on which
// node. event e owns an independent substream seeded by splitmix64(seed ^ mix(e)).
__host__ __device__ inline double eventObservable(uint64_t seed, uint64_t e) {
  uint64_t c = splitmix64(seed ^ (e * 0x9E3779B97F4A7C15ULL)); // per-event substream
  double cosT = 2.0 * u01(splitmix64(c++)) - 1.0;
  double s = 10000.0, alphaS = 0.118;            // sHat=(100 GeV)^2
  double t = -0.5 * s * (1.0 - 0.98 * cosT);
  double acc = gg2gg_sigma(s, t, -s - t, alphaS);
  int nExtra = 1 + (int)(splitmix64(c++) % 16);  // VARIABLE draw count per event
  for (int i = 0; i < nExtra; ++i) acc += u01(splitmix64(c++));
  return acc;
}

// out[i] = observable of event eventIds[i]. Same kernel for canonical, shuffled,
// and per-shard runs -> the computation is identical; only the id list differs.
__global__ void eventKernel(uint64_t seed, const uint64_t* eventIds, double* out, int n) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) out[i] = eventObservable(seed, eventIds[i]);
}

#define CK(call) do { cudaError_t e_=(call); if(e_!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); return 1; } } while(0)

static int runIds(uint64_t seed, const std::vector<uint64_t>& ids, std::vector<double>& out) {
  int n = (int)ids.size(); out.resize(n);
  uint64_t *dIds; double *dOut;
  CK(cudaMalloc(&dIds, n * sizeof(uint64_t))); CK(cudaMalloc(&dOut, n * sizeof(double)));
  CK(cudaMemcpy(dIds, ids.data(), n * sizeof(uint64_t), cudaMemcpyHostToDevice));
  int th = 256, bl = (n + th - 1) / th;
  eventKernel<<<bl, th>>>(seed, dIds, dOut, n);
  CK(cudaDeviceSynchronize());
  CK(cudaMemcpy(out.data(), dOut, n * sizeof(double), cudaMemcpyDeviceToHost));
  cudaFree(dIds); cudaFree(dOut);
  return 0;
}

int main(int argc, char** argv) {
  int N = (argc > 1) ? atoi(argv[1]) : 2000000;
  uint64_t seed = 0x5EED1234ULL;

  // Canonical run: events [0, N) in order.
  std::vector<uint64_t> ids(N);
  for (int e = 0; e < N; ++e) ids[e] = e;
  std::vector<double> canon;
  if (runIds(seed, ids, canon)) return 1;

  // (1) Regenerate a SHUFFLED subset out of order, independently.
  const int K = 100000;
  std::vector<uint64_t> sub(K);
  uint64_t r = 0xC0FFEEULL;
  for (int k = 0; k < K; ++k) { r = splitmix64(r); sub[k] = r % (uint64_t)N; }
  std::vector<double> regen;
  if (runIds(seed, sub, regen)) return 1;
  double maxDiff1 = 0.0;
  for (int k = 0; k < K; ++k) maxDiff1 = fmax(maxDiff1, fabs(regen[k] - canon[sub[k]]));

  // (2) Partition events across S "nodes" (strided) and reassemble.
  const int S = 8;
  double maxDiff2 = 0.0;
  for (int s = 0; s < S; ++s) {
    std::vector<uint64_t> shard;
    for (int e = s; e < N; e += S) shard.push_back(e);
    std::vector<double> sout;
    if (runIds(seed, shard, sout)) return 1;
    for (size_t j = 0; j < shard.size(); ++j)
      maxDiff2 = fmax(maxDiff2, fabs(sout[j] - canon[shard[j]]));
  }

  printf("Counter-based per-event reproducible RNG\n");
  printf("  events                         = %d\n", N);
  printf("  (1) %d events regenerated out-of-order, max |diff| vs canonical = %.1e\n", K, maxDiff1);
  printf("  (2) %d-way node partition reassembled, max |diff| vs canonical  = %.1e\n", S, maxDiff2);
  bool ok = (maxDiff1 == 0.0) && (maxDiff2 == 0.0);
  printf("VALIDATION: %s (any event regenerates bit-identically, independent of order/partition)\n",
         ok ? "PASS" : "FAIL");
  printf("  -> stock Pythia (RANMAR) cannot do this in O(1): event N needs N advances or a state restore.\n");
  return ok ? 0 : 2;
}
