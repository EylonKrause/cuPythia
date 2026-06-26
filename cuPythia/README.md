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
| 03 | Fused resident gg→gg MC | Simpson quadrature | relerr 3.7e-5 — PASS | **12.2×** (recip-opt) |
| 04 | Multi-GPU / multi-node MC | exact grid coverage + Simpson | identical σ @ 1/4/16 shards — PASS | **~N× per GPU** |
| 05 | Reproducible per-event RNG | out-of-order regen + node partition | max\|diff\|=0 — PASS | (reproducibility, not speed) |
| 06 | FP32 / mixed precision | FP64 + Simpson (same samples) | relerr 4e-5 — PASS | FP32 **~10×** over (recip-opt) FP64 |
| 07 | Unweighting + LHE output | η==⟨w⟩/w_max + σ vs Simpson | η=10%, relerr 3.5e-5 — PASS | (production metric + I/O) |
| 08 | QCD 2→2 process library | Pythia verbatim vs textbook (5 processes) | all PASS <1e-12 | (physics coverage) |
| 09 | Neutrino DIS (parton model) | flat vs (1−y)² + σ(νq)/σ(νq̄)=3 | PASS | (EW / ν sector) |
| 10 | CUB on-GPU event compaction | CUB count == atomic count + σ | PASS | (scalable I/O, CUDA lib) |
| 11 | VEGAS importance sampling | η 10%→76% + integral vs Simpson | PASS | **7.6× efficiency** |
| 12 | 2→2 phase-space generation | 4-mom conservation + on-shell + σ | PASS | (event kinematics) |
| 13 | O(N²) hadronic rescattering | GPU all-pairs == CPU (exact) | PASS | (heavy-ion) |
| 14 | Batched parton shower (Sudakov) | no-emission == exp(−CL) + mean mult | PASS | (sequential→batched) |

\* throughput vs 1 CPU thread. 00/01 do many trials per thread in registers
(compute-bound → big speedup); 02 transfers SoA arrays + 1 eval/thread
(transfer/FP64-bound → small e2e); 03 fuses generation+ME on-device, recovering
1.3×→**12.2×** (after the reciprocal-precompute ME optimization, see below); 04
shards the MC across GPUs/nodes (near-linear — the one place a cluster scales).

Shared: `common/rng.cuh` — host/device-identical SplitMix64 (counter-based, no
cuRAND), which is what makes the exact GPU-vs-CPU validation and independent
per-shard substreams possible.

## Optimizations applied
- **Reciprocal-precompute in the gg→gg ME** (kernels 02–12): the ME did 13 FP64
  divisions; precomputing 1/sH,1/tH,1/uH cuts it to 3. FP64 division is ~1/64 rate
  on consumer Blackwell and nvcc won't rewrite `a/b → a·(1/b)` without fast-math,
  so this is **~2.9× on the compute-bound ME** (kernel 03 6.8→12.2×; kernel 06 FP64
  7.3e8→2.1e9/s), result preserved to ~2e-15. Memory/PCIe-bound kernel 02 is
  unchanged (its bottleneck isn't division).
- **`sincos`** in the phase-space kernel (12) — one call instead of 2×sin+2×cos.
- **Benchmarking honesty:** clean numbers need an **idle** GPU; under contention
  timings vary ~10× (a first reciprocal attempt under load looked like a *regression*
  and was wrongly reverted before clean min-of-many timing settled it). Correctness
  is always validated by `make check`.
- Reviewed and **deliberately not applied**: FP32-narrowing in the FP64 MC kernels
  (erodes the tolerance margin), loop-invariant 1/s hoist (marginal after the
  reciprocal change; invasive), and several "optimizations" `-O3` already performs.

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
- [x] 12 — 2→2 phase-space generation (full final-state four-momenta)
- [x] 13 — O(N²) hadronic rescattering (heavy-ion all-pairs)
- [x] 14 — batched parton shower (Sudakov veto vs analytic)
- [ ] learned sampler (normalizing flow) past VEGAS; full veto shower w/ recoil+colour

Targets/ordering are refined by the Pythia 8.317 subsystem study (`../AUDIT.md`).
