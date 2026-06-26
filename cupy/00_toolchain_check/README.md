# 00 — Toolchain check (GPU Monte Carlo π)

Proves the CUDA toolchain + GPU work end-to-end, and establishes the
host/device-identical **SplitMix64** counter-based RNG that later cuPy kernels
reuse (no cuRAND dependency). π has a known answer, so correctness is
unambiguous.

## Result — RTX 5050 Laptop (SM 12.0, 20 SMs), CUDA 13.3

| | π estimate | error | samples | time | rate (samples/s) |
|---|---|---|---|---|---|
| GPU | 3.141699 | 1.07e-4 | 1.05e9 | 84 ms | 1.25e10 |
| CPU (1 thread) | 3.141696 | 1.04e-4 | 2.0e8 | 337 ms | 5.94e8 |

**Throughput speedup ≈ 21× (GPU vs 1 CPU thread).** Both agree with π to 4
digits using the *same* RNG — correctness and the RNG are validated together.

This is the honest "10–100× on an embarrassingly-parallel Monte-Carlo kernel"
regime — **not** an end-to-end generator speedup (Amdahl's law caps that).

## Build / run

```bash
nvcc -O3 -arch=sm_120 -o mc_pi mc_pi.cu   # sm_120 = Blackwell / RTX 5050
./mc_pi [samplesPerThread]                # default 4000  -> ~1.05e9 samples
```
