# cuPythia kernels

GPU-accelerated, **validated** building blocks toward CUDA Pythia. Every kernel
is checked against an analytic or CPU reference and benchmarked on an RTX 5050
(SM 12.0, CUDA 13.3). Speedups are per-kernel constant factors — Amdahl-bounded
end-to-end, **not** whole-generator claims.

| # | kernel | validated against | result | speedup* |
|---|---|---|---|---|
| 00 | Monte-Carlo π | π (known constant) | 3.141699, err 1.1e-4 — PASS | ~21× |
| 01 | σ(e⁺e⁻→μ⁺μ⁻) | 4πα²/3s (closed form) | 0.868544 nb, relerr 8.5e-7 — PASS | ~17.7× |
| 02 | QCD gg→gg ME (Pythia `Sigma2gg2gg`) | CPU port + textbook analytic | relerr 3e-16 / 8e-16 — PASS | 4.5× kern / **1.3× e2e** |
| 03 | Fused resident gg→gg MC | Simpson quadrature | relerr 3.7e-5 — PASS | **6.8×** |

\* throughput vs 1 CPU thread. 00/01 do many trials per thread in registers
(compute-bound → big speedup); 02 transfers SoA arrays + 1 eval/thread
(transfer/FP64-bound → small e2e speedup); 03 fuses generation+ME on-device,
recovering 1.3×→6.8× (now FP64-division-bound, not transfer-bound). The central
lesson: **keep data GPU-resident; fuse stages; mind FP64 on consumer GPUs.**

Shared: `common/rng.cuh` — host/device-identical SplitMix64 (counter-based, no
cuRAND), which is what makes the exact GPU-vs-CPU validation possible.

## Build / validate
```bash
make            # build all kernels   (override ARCH=... / NVCC=...)
make check      # build + run every kernel's self-validation
```

## Roadmap (sequenced by GPU-friendliness)
- [x] 00 — toolchain + RNG harness (MC π)
- [x] 01 — matrix-element MC (fixed-angle σ)
- [ ] 02 — 2→2 phase-space generation (full final-state kinematics)
- [ ] 03 — O(N²) hadronic rescattering (all-pairs, heavy-ion)
- [ ] 04 — batched-across-events parton shower (the hard one)

Targets/ordering are refined by the subsystem study of Pythia 8.317.
