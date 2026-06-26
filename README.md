# cuPy — CUDA Pythia

GPU-acceleration experiments on top of [Pythia 8](https://pythia.org), the
Monte-Carlo event generator. The name is short for **cu**da **Py**thia.

> Note: this is unrelated to the [CuPy](https://cupy.dev) NumPy-on-GPU library —
> the name collision is incidental.

## What this is

A study + prototyping repo that:

1. Vendors an unmodified copy of **Pythia 8.317** (`pythia8317/`) as the baseline.
2. Studies which parts of event generation are genuinely data-parallel.
3. Implements those parts as CUDA kernels under `cupy/`, each one **validated**
   against the CPU/analytic result and **benchmarked** against the Pythia baseline.

## Honest scope (read this first)

Pythia event generation is a **sequential causal chain** — hard process →
parton shower → hadronization → decays — where each stage consumes the previous
stage's output. GPUs accelerate *data-parallel* work, not causal chains, so
**no "exponential" or whole-generator speedup is physically possible.** By
Amdahl's law the end-to-end gain is bounded by the sequential fraction.

What *is* real and worth doing:

- **10–100× on individual data-parallel kernels** (Monte-Carlo phase-space /
  cross-section integration, matrix-element evaluation, O(N²) hadronic
  rescattering) measured in isolation.
- A **meaningful end-to-end factor** by batching across *events* (not within an
  event) and keeping data **GPU-resident across stages** to avoid PCIe traffic —
  the thing that capped prior GPU-rescattering work at ~3×.

Every speedup claim in this repo comes with a reproducible benchmark and a
correctness check. If a number looks too good, it's a measurement bug.

## Baseline status

- Pythia 8.317 downloaded from pythia.org (sha256 `1ae551d1…745adf`), configured
  with `g++ 13.3 / -O2 / -std=c++11`, built to `libpythia8.{a,so}`.
- Verified: `examples/main101` runs, producing a charged-multiplicity
  distribution (mean ≈ 184 charged/event, 100 pp events @ default settings).
- Toolchain: WSL2 Ubuntu 24.04, CUDA 13.3, NVIDIA RTX 5050 Laptop GPU.

## Layout

```
pythia8317/      vendored upstream Pythia 8.317 (build artifacts gitignored)
cupy/            GPU kernels + benchmarks (added incrementally)
README.md
```

## Build the baseline

```bash
cd pythia8317 && ./configure && make -j"$(nproc)"
cd examples && make main101 && ./main101
```
