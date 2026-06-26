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
| 04 | Multi-GPU / multi-node MC | exact grid coverage + Simpson | identical σ @ 1/4/16 shards — PASS | **~N× per GPU** |
| 05 | Reproducible per-event RNG | out-of-order regen + node partition | max\|diff\|=0 — PASS | (reproducibility, not speed) |
| 06 | FP32 / mixed precision | FP64 + Simpson (same samples) | relerr 4e-5 — PASS | **30× over FP64** |
| 07 | Unweighting + LHE output | η==⟨w⟩/w_max + σ vs Simpson | η=10%, relerr 3.5e-5 — PASS | (production metric + I/O) |
| 08 | QCD 2→2 process library | Pythia verbatim vs textbook (5 processes) | all PASS <1e-12 | (physics coverage) |
| 09 | Neutrino DIS (parton model) | flat vs (1−y)² + σ(νq)/σ(νq̄)=3 | PASS | (EW / ν sector) |
| 10 | CUB on-GPU event compaction | CUB count == atomic count + σ | PASS | (scalable I/O, CUDA lib) |
| 11 | VEGAS importance sampling | η 10%→76% + integral vs Simpson | PASS | **7.6× efficiency** |

\* throughput vs 1 CPU thread. 00/01 do many trials per thread in registers
(compute-bound → big speedup); 02 transfers SoA arrays + 1 eval/thread
(transfer/FP64-bound → small e2e); 03 fuses generation+ME on-device, recovering
1.3×→6.8× (now FP64-division-bound); 04 shards the MC across GPUs/nodes
(embarrassingly parallel → near-linear, the one place a cluster scales).

Shared: `common/rng.cuh` — host/device-identical SplitMix64 (counter-based, no
cuRAND), which is what makes the exact GPU-vs-CPU validation and independent
per-shard substreams possible.

## Build / validate
```bash
make            # build all kernels   (override ARCH=... / NVCC=...)
make check      # build + run every kernel's self-validation
make mpi        # optional multi-node build of kernel 04 (needs mpicxx)
```

## Roadmap
- [x] 00 — toolchain + counter-based RNG harness (MC π)
- [x] 01 — matrix-element MC (e⁺e⁻→μ⁺μ⁻ vs closed form)
- [x] 02 — batched QCD 2→2 ME (the transfer-bound lesson)
- [x] 03 — fused resident ME (the resident win; FP64 ceiling)
- [x] 04 — multi-GPU / multi-node scaling (cluster)
- [x] 05 — counter-based per-event reproducible RNG (GRID production / debugging)
- [x] 06 — FP32 / mixed precision (30× over FP64 at the same MC accuracy)
- [x] 07 — GPU unweighting efficiency + Les Houches (.lhe) output
- [x] 08 — QCD 2→2 process library (gg→gg, qg→qg, qq'→qq', qqbar→gg, gg→qqbar)
- [x] 09 — neutrino DIS parton-model cross sections (EW sector, (1−y)² signature)
- [x] 10 — on-GPU event compaction with CUB (scalable unweighted-event I/O)
- [x] 11 — VEGAS adaptive importance sampling (η 10%→76%)
- [ ] full 2→2 phase-space generation (final-state kinematics)
- [ ] O(N²) hadronic rescattering (all-pairs, heavy-ion)
- [ ] batched-across-events parton shower (the hard one)

Targets/ordering are refined by the Pythia 8.317 subsystem study (`../AUDIT.md`).
