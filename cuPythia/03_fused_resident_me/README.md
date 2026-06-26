# 03 — Fused, GPU-resident gg→gg (the payoff)

Each thread generates its own trials **on-device** (counter-based RNG, many per
thread, in registers), evaluates the verbatim Pythia gg→gg ME, and accumulates a
Monte-Carlo integral. Nothing crosses PCIe but one scalar — the opposite of
kernel 02. Physics validated against a deterministic CPU Simpson quadrature of
the identical integrand (RNG-independent).

## Result — RTX 5050, 1.05×10⁹ trials, |cosθ| < 0.9

| | σ (mb) | note |
|---|---|---|
| Simpson quadrature | 1.323768e-4 | deterministic reference |
| GPU fused MC | 1.323719e-4 | **relerr 3.7e-5** |
| CPU MC (same RNG) | 1.324063e-4 | |

`VALIDATION: PASS`. GPU 7.3×10⁸/s vs CPU 1.07×10⁸/s → **6.8×**.

## What this shows (the project's thesis, in numbers)
- **Resident/fusion removes the transfer wall:** end-to-end **1.3× (kernel 02) →
  6.8×**. The kernel is now compute-bound, not PCIe-bound.
- **The new ceiling is FP64 division:** the ME does ~15 double divisions/trial,
  and consumer Blackwell (RTX 5050) runs FP64 at ~1/64 of FP32. Kernels 00/01 hit
  17–21× precisely because their per-trial work is division-light.
- **Next levers:** FP32 / mixed precision where the physics tolerates it (~30–60×
  more FP throughput on this GPU), and reciprocal precompute to cut divisions.
  On a data-center GPU (full-rate FP64) this same kernel would already scale far
  better — the 6.8× is a laptop-FP64 artifact, not an algorithmic limit.

## Build / run
```bash
nvcc -O3 -arch=sm_120 -o fused_gg2gg fused_gg2gg.cu
./fused_gg2gg [trialsPerThread=4000]
```
