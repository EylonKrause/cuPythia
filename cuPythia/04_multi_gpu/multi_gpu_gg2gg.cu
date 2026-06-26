// cuPythia kernel 04 — MULTI-GPU / MULTI-NODE Monte-Carlo (g g -> g g).
//
// Why this scales (honestly): MC event generation is embarrassingly parallel —
// independent RNG substreams + one tiny reduction at the end. This is the ONE
// place a cluster gives ~linear speedup; it does NOT speed up the sequential
// within-event chain (shower/hadronization).
//
// Decomposition: a fixed global thread grid (totalBlocks x threads) is split into
// S disjoint SHARDS by contiguous block ranges. Every thread keys its RNG off its
// GLOBAL thread id, so the union of all shards == one deterministic run: the total
// sample COUNT is identical for any S (only the FP reduction order differs). That
// is what lets a single GPU validate the cluster reduction — run S=1 vs S=16 and
// the count is identical and sigma is unchanged.
//
// Mapping:
//   - shards map round-robin to the node's GPUs (cudaGetDeviceCount), one host
//     thread per shard; distinct GPUs run concurrently -> throughput scales ~#GPU.
//   - with -DUSE_MPI, each rank owns a disjoint shard range on its local GPUs and
//     MPI_Allreduce combines the partial sums + counts across nodes.
//
// Build (single node, validatable here):
//   nvcc -O3 -arch=sm_120 -o multi_gpu_gg2gg multi_gpu_gg2gg.cu
//   ./multi_gpu_gg2gg [trialsPerThread=2000] [--shards S]
// Build (cluster, one rank per node/GPU):
//   nvcc -O3 -arch=sm_120 -DUSE_MPI -ccbin mpicxx -o multi_gpu_gg2gg \
//        multi_gpu_gg2gg.cu -lmpi
//   mpirun -np <ranks> ./multi_gpu_gg2gg 2000        # or srun under Slurm

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>
#include <vector>
#include <thread>
#include <cuda_runtime.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#include "../common/rng.cuh"

__host__ __device__ inline double pow2(double x) { return x * x; }

// verbatim Pythia 8.317 Sigma2gg2gg::sigmaKin (src/SigmaQCD.cc:115-129) = dσ/dt̂
__host__ __device__ inline double gg2gg_sigma(double sH, double tH, double uH, double alpS) {
  double is=1.0/sH, it=1.0/tH, iu=1.0/uH;  // reciprocal precompute: 3 FP64 div instead of 13
  double rts=tH*is, rst=sH*it, rus=uH*is, rsu=sH*iu, rtu=tH*iu, rut=uH*it;
  double sigTS = (9. / 4.) * (rts*rts + 2.*rts + 3. + 2.*rst + rst*rst);
  double sigUS = (9. / 4.) * (rus*rus + 2.*rus + 3. + 2.*rsu + rsu*rsu);
  double sigTU = (9. / 4.) * (rtu*rtu + 2.*rtu + 3. + 2.*rut + rut*rut);
  return (M_PI * is*is) * pow2(alpS) * 0.5 * (sigTS + sigUS + sigTU);
}

// one shard = a contiguous range of global blocks; RNG keyed off GLOBAL thread id.
__global__ void shardKernel(uint64_t seed, uint64_t blockOffset, uint64_t kPer,
                            double s, double alpS, double cMax, double* partial) {
  uint64_t gBlock = blockOffset + blockIdx.x;
  uint64_t gTid = gBlock * (uint64_t)blockDim.x + threadIdx.x;
  uint64_t ctr = seed + gTid * 0x100000001B3ULL;
  double local = 0.0;
  for (uint64_t i = 0; i < kPer; ++i) {
    double c = (2.0 * u01(splitmix64(ctr++)) - 1.0) * cMax;
    double t = -0.5 * s * (1.0 - c);
    local += gg2gg_sigma(s, t, -s - t, alpS);
  }
  atomicAdd(partial, local);
}

static double simpson(int N, double s, double alpS, double cMax) {
  double a = -cMax, h = (2.0 * cMax) / N, sum = 0.0;
  for (int i = 0; i <= N; ++i) {
    double c = a + i * h, t = -0.5 * s * (1.0 - c);
    double w = (i == 0 || i == N) ? 1.0 : (i % 2 ? 4.0 : 2.0);
    sum += w * gg2gg_sigma(s, t, -s - t, alpS);
  }
  return (sum * h / 3.0) * (s / 2.0);
}

using Clock = std::chrono::steady_clock;
static double msec(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

struct ShardOut { double sum = 0.0; double ms = 0.0; uint64_t blocks = 0; int device = 0; };

static void runShard(int device, uint64_t blockOffset, uint64_t blocksThis, int threads,
                     uint64_t seed, uint64_t kPer, double s, double alpS, double cMax,
                     ShardOut* out) {
  cudaSetDevice(device);
  double* d; cudaMalloc(&d, sizeof(double)); cudaMemset(d, 0, sizeof(double));
  shardKernel<<<1, threads>>>(seed, blockOffset, 1, s, alpS, cMax, d); // warmup
  cudaDeviceSynchronize(); cudaMemset(d, 0, sizeof(double));
  auto t0 = Clock::now();
  shardKernel<<<(unsigned)blocksThis, threads>>>(seed, blockOffset, kPer, s, alpS, cMax, d);
  cudaDeviceSynchronize();
  auto t1 = Clock::now();
  double h = 0.0; cudaMemcpy(&h, d, sizeof(double), cudaMemcpyDeviceToHost); cudaFree(d);
  out->sum = h; out->ms = msec(t0, t1); out->blocks = blocksThis; out->device = device;
}

int main(int argc, char** argv) {
  int rank = 0, nranks = 1;
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);
#endif

  uint64_t kPer = 2000; int shardsOverride = 0;
  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i], "--shards") && i + 1 < argc) shardsOverride = atoi(argv[++i]);
    else kPer = strtoull(argv[i], nullptr, 10);
  }

  int nDev = 0; cudaGetDeviceCount(&nDev); if (nDev < 1) nDev = 1;
  const int threads = 256;
  const uint64_t totalBlocks = 4096;          // fixed global grid (1,048,576 threads)
  // total shards: explicit override, else one per (rank, local GPU).
  int totalShards = shardsOverride ? shardsOverride : nranks * nDev;
  uint64_t bps = (totalBlocks + totalShards - 1) / totalShards; // ceil

  // this rank owns a contiguous slice of the shard list.
  int shardsPerRank = (totalShards + nranks - 1) / nranks;
  int myFirst = rank * shardsPerRank;
  int myLast = (myFirst + shardsPerRank < totalShards) ? myFirst + shardsPerRank : totalShards;

  double s = 100.0 * 100.0, alphaS = 0.118, cMax = 0.9, conv = 0.3893793721; // GeV^-2 -> mb
  uint64_t seed = 0xABCDEF01ULL;

  std::vector<ShardOut> outs(myLast > myFirst ? myLast - myFirst : 0);
  std::vector<std::thread> ths;
  for (int g = myFirst; g < myLast; ++g) {
    uint64_t off = (uint64_t)g * bps;
    if (off >= totalBlocks) break;
    uint64_t blocksThis = (off + bps <= totalBlocks) ? bps : totalBlocks - off;
    int dev = (g - myFirst) % nDev;                 // local GPU, round-robin
    ths.emplace_back(runShard, dev, off, blocksThis, threads, seed, kPer,
                     s, alphaS, cMax, &outs[g - myFirst]);
  }
  auto w0 = Clock::now();
  for (auto& t : ths) t.join();
  auto w1 = Clock::now();

  double nodeSum = 0.0; uint64_t nodeBlocks = 0; double nodeMaxMs = 0.0;
  for (auto& o : outs) { nodeSum += o.sum; nodeBlocks += o.blocks; if (o.ms > nodeMaxMs) nodeMaxMs = o.ms; }

  double globalSum = nodeSum; uint64_t globalBlocks = nodeBlocks; double wallMs = msec(w0, w1);
#ifdef USE_MPI
  MPI_Allreduce(MPI_IN_PLACE, &globalSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  unsigned long long gb = nodeBlocks, gbAll = 0;
  MPI_Allreduce(&gb, &gbAll, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  globalBlocks = gbAll;
  double wmax = wallMs; MPI_Allreduce(MPI_IN_PLACE, &wmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  wallMs = wmax;
#endif

  if (rank == 0) {
    uint64_t totalSamples = globalBlocks * (uint64_t)threads * kPer;
    double sigma = (2.0 * cMax) * (s / 2.0) * (globalSum / (double)totalSamples) * conv;
    double ref = simpson(2000000, s, alphaS, cMax) * conv;
    double rel = fabs(sigma - ref) / ref;
    double rate = totalSamples / (wallMs / 1000.0);
    printf("MULTI-GPU g g -> g g   ranks=%d  GPUs/rank=%d  shards=%d  (sHat=(100 GeV)^2)\n",
           nranks, nDev, totalShards);
    printf("  global blocks  = %llu  (covers the fixed %llu-block grid exactly)\n",
           (unsigned long long)globalBlocks, (unsigned long long)totalBlocks);
    printf("  total samples  = %.3e\n", (double)totalSamples);
    printf("  Simpson ref    = %.6e mb\n", ref);
    printf("  multi-GPU MC   = %.6e mb   relerr=%.2e\n", sigma, rel);
    printf("  wall (max shard/rank) = %.1f ms   aggregate rate = %.3e samples/s\n", wallMs, rate);
    bool ok = (globalBlocks == totalBlocks) && (rel < 2e-3);
    printf("VALIDATION: %s (full grid covered exactly AND sigma matches quadrature)\n",
           ok ? "PASS" : "FAIL");
    if (nranks * nDev == 1)
      printf("  NOTE: 1 physical GPU -> extra shards validate the reduction, not speed.\n"
             "        On N distinct GPUs the shards run concurrently and rate scales ~N.\n");
#ifdef USE_MPI
    MPI_Finalize();
#endif
    return ok ? 0 : 2;
  }
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
