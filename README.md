# cuPythia — CUDA Pythia

GPU-acceleration experiments on top of [Pythia 8](https://pythia.org), the
Monte-Carlo event generator. The name is short for **cu**da **Pythia**.
(Renamed from "cuPy" to avoid confusion with the unrelated
[CuPy](https://cupy.dev) NumPy-on-GPU array library.)

## What this is

1. Vendors an unmodified copy of **Pythia 8.317** (`pythia8317/`) as the baseline.
2. Studies which parts of event generation are genuinely data-parallel (`AUDIT.md`).
3. Implements those parts as CUDA kernels under `cuPythia/`, each **validated**
   against a CPU/analytic result and **benchmarked** on an RTX 5050.

## Honest scope (read first)

Pythia event generation is a **sequential causal chain** (hard process → shower →
hadronization → decays), so **no "exponential" or whole-generator speedup is
physically possible** — Amdahl's law bounds the end-to-end gain. What *is* real:

- **10–100×** on individual data-parallel kernels in isolation;
- a meaningful **end-to-end factor** by batching across *events* and keeping data
  **GPU-resident** (avoiding PCIe traffic).

Every speedup here ships with a reproducible benchmark and a correctness check.

## Kernels (validated on RTX 5050, SM 12.0, CUDA 13.3)

| # | kernel | validated against | speedup |
|---|---|---|---|
| 00 | Monte-Carlo π | π (known) | ~21× |
| 01 | σ(e⁺e⁻→μ⁺μ⁻) | 4πα²/3s | ~17.7× |
| 02 | QCD gg→gg ME (Pythia `Sigma2gg2gg`) | CPU port + textbook | 4.5× kern / 1.3× e2e |
| 03 | fused resident gg→gg MC | Simpson quadrature | ~6.8× |

See `cuPythia/README.md`. The 02→03 jump (1.3×→6.8×) is the project's core lesson:
keep data GPU-resident; the remaining ceiling is consumer-GPU FP64 throughput.

## Audit

A 32-agent study of Pythia 8.317 (`AUDIT.md`) surfaced **17 adversarially-verified
issues**. The two real correctness bugs are fixed here (`SigmaProcess.cc:1171`
factor-scale copy-paste; `Basics.cc` `Rndm::pick` OOB); the rest are documented as
upstream-PR candidates.

## Layout

```
pythia8317/      vendored Pythia 8.317 (build artifacts gitignored)
cuPythia/        GPU kernels + benchmarks (00..03, common/rng.cuh, Makefile)
AUDIT.md         verified Pythia findings
```

## Build

```bash
cd pythia8317 && ./configure && make -j"$(nproc)"   # baseline library
cd ../cuPythia && make check                         # build + validate all kernels
```
