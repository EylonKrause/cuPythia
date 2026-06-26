# 04 — Multi-GPU / multi-node Monte-Carlo (gg→gg)

Scales the gg→gg MC integral across many GPUs and nodes. MC generation is
**embarrassingly parallel** (independent RNG substreams + one final reduction),
so this scales ~linearly with the number of GPUs — the one place a cluster gives
near-ideal speedup. It does **not** speed up the sequential within-event chain
(shower/hadronization); that stays per-event.

## How it works
- A fixed global thread grid is split into **S disjoint shards** (contiguous block
  ranges). Each thread keys its RNG off its **global** thread id, so the union of
  shards == one deterministic run: the total sample **count is identical for any S**.
- Shards map round-robin to the node's GPUs (`cudaGetDeviceCount`), one host thread
  per shard. Distinct GPUs run concurrently → throughput scales ~#GPU.
- `-DUSE_MPI`: each rank owns a disjoint shard range on its local GPUs;
  `MPI_Allreduce` combines the partial sums + counts across nodes (one rank per
  node, or per GPU).

## Validated here (1× RTX 5050)
Run at **1, 4, 16 shards** — each covers the 4096-block grid exactly, same
2.097×10⁹ samples, **σ = 1.323804×10⁻⁴ mb (identical), relerr 2.7e-5 vs Simpson,
PASS**. Wall time is ~constant (one physical GPU serializes the shards), which is
the point: it validates the **decomposition + reduction**, not speed. On N distinct
GPUs the shards run concurrently and the aggregate rate scales ~N×.

## What is NOT validated here (honest)
This box has 1 GPU and no MPI/sudo, so real multi-GPU concurrency and the MPI
multi-node path are **not run here**. The MPI layer uses only stable MPI core
(`Init`/`Comm_rank`/`Comm_size`/`Allreduce`/`Finalize`) and reuses the exact
validated shard+reduce logic — only the final cross-rank sum changes. Validate
the speedup on real multi-GPU / multi-node hardware.

## Build / run
Single node (uses all local GPUs):
```bash
nvcc -O3 -arch=sm_120 -Xcompiler -pthread -o multi_gpu_gg2gg multi_gpu_gg2gg.cu
./multi_gpu_gg2gg 2000 [--shards S]
```
Cluster (one rank per node or per GPU; needs MPI + CUDA toolchain):
```bash
nvcc -O3 -arch=sm_120 -Xcompiler -pthread -DUSE_MPI -ccbin mpicxx \
     -o multi_gpu_gg2gg multi_gpu_gg2gg.cu -lmpi      # or: make mpi
mpirun -np <ranks> ./multi_gpu_gg2gg 2000             # or: srun ... under Slurm
```
