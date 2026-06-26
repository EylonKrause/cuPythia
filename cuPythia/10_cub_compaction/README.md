# 10 — On-GPU event compaction with CUB

A concrete example of a **CUDA toolkit library helping cuPythia**. Kernel 07
filtered accepted (unweighted) events on the **host**; here we generate candidates
on the GPU, flag the accepted ones, and use **`cub::DeviceSelect::Flagged`** to
compact them into a dense array **entirely on the GPU** — no host round-trip. This
is the scalable unweighted-event I/O pattern modern GPU generators (madgraph4gpu)
rely on. CUB ships header-only with the CUDA toolkit, so it builds on Windows and
Linux with no extra dependency.

## Result — RTX 5050
- 16,777,216 candidates → **1,680,778 accepted, CUB-compacted on-GPU**
- CUB count **== independent atomic count** (exact) ⇒ compaction is correct
- unweighting efficiency 10.0%, σ from the compacted set matches Simpson (relerr 2.4e-4)
- `VALIDATION: PASS`

## Build / run
```bash
nvcc -O3 -arch=sm_120 -o cub_compaction cub_compaction.cu
./cub_compaction [nCandidates=16777216]
```
