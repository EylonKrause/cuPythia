# cuPythia — CUDA Pythia

GPU-acceleration experiments on top of [Pythia 8](https://pythia.org), the
Monte-Carlo event generator. The name is short for **cu**da **Pythia**.
(Not affiliated with the [CuPy](https://cupy.dev) GPU array library.)

## What this is

1. Vendors an unmodified copy of **Pythia 8.317** (`pythia8317/`) as the baseline.
2. Studies which parts of event generation are genuinely data-parallel (`AUDIT.md`).
3. Implements those parts as CUDA kernels under `cuPythia/`, each **validated**
   against a CPU/analytic result, **benchmarked** on an RTX 5050, and scalable
   across **multiple GPUs and nodes**.

## Honest scope (read first)

Pythia event generation is a **sequential causal chain** (hard process → shower →
hadronization → decays), so **no "exponential" or whole-generator speedup is
physically possible** — Amdahl's law bounds the end-to-end gain. What *is* real:

- **10–100×** on individual data-parallel kernels in isolation;
- a meaningful **end-to-end factor** by batching across *events* and keeping data
  **GPU-resident** (avoiding PCIe traffic);
- **near-linear** scaling of Monte-Carlo generation across a **cluster of GPUs**
  (independent RNG substreams + one reduction).

Every speedup here ships with a reproducible benchmark and a correctness check.

## Kernels (validated on RTX 5050, SM 12.0, CUDA 13.3)

| # | kernel | validated against | speedup |
|---|---|---|---|
| 00 | Monte-Carlo π | π (known) | ~21× |
| 01 | σ(e⁺e⁻→μ⁺μ⁻) | 4πα²/3s | ~17.7× |
| 02 | QCD gg→gg ME (Pythia `Sigma2gg2gg`) | CPU port + textbook | 4.5× kern / 1.3× e2e |
| 03 | fused resident gg→gg MC | Simpson quadrature | ~6.8× |
| 04 | multi-GPU / multi-node MC | exact grid coverage + quadrature | ~N× per GPU |
| 05 | reproducible per-event RNG | out-of-order regen + node partition | bit-identical |
| 06 | FP32 / mixed precision | FP64 + Simpson | **30× over FP64**, same accuracy |

See `cuPythia/README.md`. The 02→03 jump (1.3×→6.8×) is the core lesson (keep data
GPU-resident); 04 scales MC across a cluster (near-linear — the one place that holds).

## What this adds beyond stock Pythia

See `HEP_FEATURES.md` — capabilities the HL-LHC / generators-on-accelerators
community wants that stock Pythia lacks: GPU-accelerated ME/MC, cluster scaling,
and **counter-based per-event reproducible RNG** (regenerate any single event
independently on any node — useful for GRID production and debugging).

## Audit

A 32-agent study of Pythia 8.317 (`AUDIT.md`) surfaced **17 adversarially-verified
issues**. The two real correctness bugs are fixed here (`SigmaProcess.cc:1171`
factor-scale copy-paste; `Basics.cc` `Rndm::pick` OOB); the rest are documented as
upstream-PR candidates.

## Layout

```
pythia8317/      vendored Pythia 8.317 (build artifacts gitignored)
cuPythia/        GPU kernels 00..05 + common/rng.cuh + Makefile
AUDIT.md         verified Pythia findings
HEP_FEATURES.md  gap analysis vs community needs
```

## Build

```bash
cd pythia8317 && ./configure && make -j"$(nproc)"   # baseline library
cd ../cuPythia && make check                         # build + validate all kernels
cd ../cuPythia && make mpi                            # optional: multi-node build (needs mpicxx)
```
