# 06 — FP32 / mixed-precision (breaking the FP64 ceiling)

Kernel 03 topped out ~7–9×, FP64-division-bound (consumer Blackwell runs FP64 at
~1/64 of FP32). This runs the **same** fused gg→gg MC in both `double` and `float`
with **identical RNG samples**, so the comparison is apples-to-apples.

## Result — RTX 5050, 1.05×10⁹ trials
| precision | σ (mb) | relerr vs Simpson | time | rate |
|---|---|---|---|---|
| FP64 | 1.323719e-4 | 3.7e-5 | 1440 ms | 7.3e8/s |
| FP32 | 1.323712e-4 | 4.2e-5 | **48 ms** | **2.2e10/s** |

**FP32/FP64 speedup = 30×**, `PASS`. The key point: **FP32 gives the same physics
accuracy** here — the Monte-Carlo statistical error (~1/√N ≈ 1e-4) dwarfs FP32
rounding, so both match the quadrature to ~4e-5. For MC integration on a consumer
GPU, single precision is the right default; reserve FP64 for the few spots that
genuinely need it (mixed precision).

## Caveat
FP32 is appropriate where the integrand is well-conditioned and MC error dominates
(true here). Near the t̂,û→0 collinear poles, or for observables sensitive to
catastrophic cancellation, validate FP32 against FP64 before trusting it — which is
exactly what this kernel does.

## Build / run
```bash
nvcc -O3 -arch=sm_120 -o fp32_gg2gg fp32_gg2gg.cu
./fp32_gg2gg [trialsPerThread=4000]
```
