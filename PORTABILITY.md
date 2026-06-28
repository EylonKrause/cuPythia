# cuPythia GPU portability — Pascal (2016) through Blackwell (2025)

cuPythia's CUDA kernels are written in **plain FP64 + integer math**: no architecture-specific
intrinsics (no warp-shuffle, cooperative groups, tensor cores, `__half`/bf16, or `__CUDA_ARCH__`
branching), and the only hardware feature they rely on — double-precision `atomicAdd` — is native
on **Pascal (sm_60) and up**, with a CAS fallback (`common/rng.cuh`) compiled only for older arches.
So the code runs unchanged on every NVIDIA GPU from **Pascal to Blackwell**; portability is almost
entirely a matter of **build flags** and **which CUDA toolkit you use**.

## The one thing to know: CUDA 13 removed Pascal/Volta

NVIDIA drops old architectures from each new toolkit. To *emit* code for an architecture you need a
toolkit that still supports it:

| GPU family | Arch (`sm_`) | Example cards | Min CUDA | **Max CUDA that can target it** |
|---|---|---|---|---|
| Pascal   | 60, 61 | Tesla P100, GTX 1080/1070, Titan X | 8.0  | **12.9** (removed in 13.0) |
| Volta    | 70     | Tesla V100, Titan V                | 9.0  | **12.9** (removed in 13.0) |
| Turing   | 75     | RTX 2080, T4, GTX 16xx             | 10.0 | 13.x (still supported) |
| Ampere   | 80, 86 | A100, RTX 3090/3080               | 11.0 (sm_86: 11.1) | 13.x |
| Ada      | 89     | RTX 4090                          | 11.8 | 13.x |
| Hopper   | 90     | H100                              | 11.8 | 13.x |
| Blackwell| 120    | RTX 50-series                     | 12.8 | 13.x |

So: **for Pascal/Volta, build with a CUDA 11.x or 12.x toolkit** (12.9 is the last). Ampere and newer
build with any toolkit ≥ their min. cuPythia requires **C++17**, i.e. **CUDA ≥ 11.0** (the CUB
event-compaction kernel needs it); this still covers all Pascal+ GPUs via CUDA 11/12.

> **Pascal is compile-verified, not just asserted.** The repo's Linux/WSL toolkit is CUDA 13.x, whose
> `nvcc --list-gpu-arch` minimum is `sm_75` (`nvcc -arch=sm_60` → *"Unsupported gpu architecture"*).
> But the same machine's **CUDA 12.9** toolkit *can* emit Pascal, and the **entire pipeline compiles
> for `sm_60`** there — base kernels (build_events, shower_fsr, hadronize, hadronize_mr, pdf_xsec,
> region_test, decay_test, zlund_test, baryon_test) **and** every opt-in physics flag
> (`-DME_FIRST -DGLUON_SPLIT -DZFLAV` shower; `-DDECAYS -DHFDECAY -DDALITZ_ME` decay chain). The
> `atomicAdd(double)` CAS shim path also compiles for `sm_50`. (CUDA 13.x is also covered: the same
> sources build for `sm_75`, the oldest 13.x emits, and `make check` stays byte-identical.)

## Building for a specific / older GPU

**Easiest — zero config (`./run.sh` at the repo root):** auto-detects your GPU's compute capability,
picks the newest installed CUDA toolkit that can compile for it (so Pascal/Volta automatically use a
CUDA ≤ 12.9 toolkit if present), builds for that arch once (cached), and runs — e.g. `./run.sh`,
`./run.sh shower_fsr 200000`, `./run.sh --arch sm_61 ...`, `./run.sh --cuda /usr/local/cuda-12.9 ...`.
On Windows, `cuPythia\build.ps1` does the same auto-detect + compatible-CUDA selection. The manual
recipes below are for fine-grained control.

**Multi-GPU & mixed-architecture clusters (automatic).** `run.sh` enumerates *every* GPU. If they have
**different** microarchitectures it builds **one fatbinary covering all of them** (e.g. an A100 sm_80 +
an RTX 4090 sm_89 → `SMS="80 89"`), and picks a CUDA toolkit that can emit them all (erroring with
guidance if none can — e.g. a Pascal+Blackwell mix needs a toolkit supporting both). A **generation run
with an output file is sharded across all GPUs**: each GPU gets a disjoint slice of the counter-RNG
event stream: the kernel adds an **event-index offset** `shard·M` to its local index `e`
(`CUPYTHIA_SHARD`/`CUPYTHIA_SHARD_N`), so shard *i* computes exactly global events `[i·M, i·M+cnt)` and
every seed is a pure function of the *global* index. Each GPU runs its own native SASS from the
fatbinary; the per-GPU dumps are merged. Because the slices are structurally disjoint and contiguous,
**the merged result is bit-identical to one GPU running all N events** — verified: a 2-shard merged
`hadronize_mr` dump is **byte-for-byte identical (same sha256)** to the single-GPU run of all N events.
With no `CUPYTHIA_SHARD` set the offset is 0 → single-GPU runs stay byte-identical. `--gpus 0,2` picks
devices; `--single` forces one. To pool GPUs across **multiple machines on a LAN** (any mix of
architectures), `cluster.sh` does the same global-index sharding over SSH and merges the per-host
dumps — see [`CLUSTER.md`](CLUSTER.md).

**Pipeline** (`cuPythia/pipeline/`, GNU make):
```
make ARCH=sm_60                 # Pascal P100   (needs CUDA <= 12.8)
make ARCH=sm_61 hadronize_mr_hf # GTX 10-series, one stage
make                            # default sm_120 (this box); override NVCC=/usr/local/cuda-12.8/bin/nvcc
```
**Portable multi-arch fatbinary** — one binary that runs on many GPUs (and JIT-runs on future ones
via the embedded PTX of the highest arch):
```
make SMS="60 61 70 75 80 86 89 90 120"      # with CUDA <= 12.8
make SMS="75 80 86 89 90 120"               # with CUDA 13.x (Turing+)
```
Point `make` at an older toolkit with `NVCC=/usr/local/cuda-12.8/bin/nvcc` (and matching
`PY8`/`HEPMC3` if you use the optional tools).

**Demo kernels + cross-platform** (`cuPythia/`, CMake ≥ 3.20):
```
cmake -S . -B build                                  # default "native": auto-detects your GPU
cmake -S . -B build -DCMAKE_CUDA_ARCHITECTURES=60    # Pascal
cmake -S . -B build -DCMAKE_CUDA_ARCHITECTURES=all-major   # portable, every major arch the toolkit has
cmake --build build
```
**Windows** (`build.ps1`):
```
powershell -ExecutionPolicy Bypass -File build.ps1                       # native arch, newest CUDA
powershell -ExecutionPolicy Bypass -File build.ps1 -Arch 61 -CudaPath "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.8"
```

## FP64 performance caveat (correct everywhere, fast only on data-center parts)

cuPythia is FP64-heavy. Double-precision throughput varies enormously by GPU:
- **GP100 / V100 / A100 / H100** (data-center): FP64 at 1/2 of FP32 — fast.
- **Consumer Pascal/Turing/Ampere/Ada/Blackwell** (GeForce): FP64 at **1/32–1/64** of FP32 — *correct
  but slow* (e.g. the shower runs at ~0.2 M evt/s on an RTX 5050). Results are bit-for-bit the same;
  only wall-clock differs. A future FP32/mixed-precision path (kernel 06 shows the approach) is the
  lever for consumer-card speed; it does not affect portability.

## What is guaranteed portable

- No arch-gated intrinsics; FP64 `atomicAdd` native on sm_60+, CAS fallback below (`rng.cuh`).
- The counter-based RNG (`splitmix64`) is integer math → bit-identical host↔device on every arch.
- File-scope physics constants (MZ, EBEAM, the Lund `H_*` params, …) are `constexpr`, so they are
  usable in device code on **every** CUDA toolkit. (Older toolkits — e.g. CUDA 12.9 — reject a
  file-scope `static const double` referenced from `__device__` code with *"identifier undefined in
  device code"*; CUDA 13.x silently accepts it. `constexpr` is the portable form and compiles on both.)
- `M_PI` fallback for MSVC; C++17 throughout (CUDA ≥ 11.0). Builds on Linux (g++) and Windows (MSVC).
