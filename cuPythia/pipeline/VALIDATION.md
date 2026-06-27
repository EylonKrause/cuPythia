# cuPythia pipeline — validation record

Every stage of the device-resident generator is gated by an explicit, reproducible
check (run `make check`). This file records *what* is validated and the *numbers*
achieved, so claims are auditable. The guiding rule (the project's standing stance):
**no stage is called "validated" without a real check, and no result is over-claimed.**

The architectural lever that makes all of this possible is the **counter-based
SplitMix64 RNG** (`../common/rng.cuh`): every random draw is a pure function of a
per-event counter, so (a) any event is O(1)-reproducible on any node, and (b) the GPU
and an *identical CPU port* consume the same draws — enabling exact GPU-vs-CPU checks.

## Stage-by-stage

| stage | file | what is checked | result |
|-------|------|-----------------|--------|
| 0 hard process | `build_events.cu` | 4-momentum conservation; record integrity; σ vs Simpson quadrature | imbalance **0**; σ relerr 5.7e-3 (MC) |
| 1 PDF convolution | `pdf_xsec.cu` | interp fidelity at the σ level; GPU-vs-CPU determinism (identical samples) | fidelity **1.7e-4**; determinism **1.4e-11** |
| 2 reweighting | `reweight.cu` | N scale weights in one pass == N independent pinned re-runs | **bit-identical (max\|diff\|=0)** |
| orchestrator | `generate.cu` | build→reweight→unweight→CUB-compact on one record, no host round-trip; σ; CUB count | σ relerr 4.8e-5; CUB count **exact** |
| 3 FSR shower | `shower_fsr.cu` | momentum cons.; on-shell; GPU re-run reproducibility; GPU-vs-CPU control flow; **thrust vs Pythia** | see below |
| 4 hadronization | `hadronize.cu` | momentum cons.; on-shell; GPU≡CPU; reproducible; **multiplicity vs Pythia** | see below |
| 5 LHE I/O | `lhe_writer.cu` | spec-valid LHEF read back by Pythia (frameType=4) | **10⁴/10⁴** read+showered, imbalance **2.3e-5** |
| 5 HepMC3 I/O | `hepmc3_writer.cc` | spec-valid HepMC3 Asciiv3, read back by HepMC3's reader | **5000/5000**, round-trip **2.9e-7** |
| bridge | `bridge.cu`+`bridge_pythia.cc` | GPU shower kinked chains hadronized by Pythia | **5000/5000**, ⟨n_ch⟩ **20.75 ≈ LEP 21** |
| 4b multi-region hadronization | `hadronize_mr.cu` | all-GPU gluon-kinked Lund: conservation, on-shell, mult vs Pythia | cons **5e-11**, mult **21.0 vs 22.57** (~7%), **3.4% drop** |
| (unit) Lund f(z) | `zlund_test.cu` | sampled z vs analytic f(z), 3 regimes | χ²/ndf **0.98/1.05/0.92/0.91** |
| (unit) multi-region table | `region_test.cu` | gluon-kinked region basis: lightlike + orthonormal + project/pHad inverse | **5.7e-16** (q-g-q̄), **8.9e-16** (q-g-g-q̄) |
| (tool) FFT luminosity | `pdf_lumi_fft.cu` | cuFFT convolution vs direct O(N²) | agree **1.7e-7**, **94×/464×** faster |

## Stage 3 (the headline) — FSR dipole shower, in detail

Setup: e⁺e⁻→Z→qq̄ at √s = 91.1876 GeV, final-state radiation only, one event per GPU
thread (the GAPS decomposition), Pythia `SimpleTimeShower` splitting kernels
$(1+z^2)/2$ and $(1+z^3)/2$, running-α_s trial generation with **flavour-threshold
matching** (n_f = 5,4,3 across m_b,m_c), and the exact Pythia local-dipole recoil.

Correctness:
- **4-momentum conservation** max|Δ| = 1.8e-9 GeV
- **on-shellness** max|p²| = 1.5e-12 GeV²
- **reproducibility** GPU re-run diffs = 0 (counter-RNG)
- **GPU vs CPU port** control flow **100% bit-identical** (every accept/veto decision);
  summed momenta agree to 1.4e-12 (GPU/CPU IEEE transcendental limit — never flips a
  branching decision)

Physics (vs Pythia's own `SimpleTimeShower`, identical setup — `thrust_pythia.cc`):
- thrust agrees to **4.0%** (mean |ratio−1| over all bins) vs Pythia with ME-corrections
  **off** — the apples-to-apples leading-log comparison — ⟨1−T⟩ **0.0690 vs 0.0698**.
- vs *default* Pythia (ME-corrections on) the spread is 23.8%; that gap **is** Pythia's
  ME corrections, a higher-order effect no LL shower carries (GAPS included).

Honest residuals / scope: FSR-only, massless, q→qg & g→gg. The ~4% LL residual is g→qq̄
(kept on in the Pythia reference because `TimeShower:nGluonToQuark=0` *hangs* Pythia
8.317), the exact-Λ α_s normalisation, and MC statistics. ME corrections (LL→NLL) are
the next step. An early far-tail (1−T>0.35) excess was traced to **the observable** (a
4-seed thrust axis missing the true axis on rare multi-jet events), not the shower — it
collapses when the axis is seeded from every particle.

## Stage 4 — Lund string hadronization, in detail

Setup: one straight u-ū colour-singlet string at √s = 91.1876 GeV, one string per GPU
thread. Faithful port of Pythia 8.317 string fragmentation: the **zLund f(z) sampler**
(validated standalone, χ²/ndf≈1 across z_max 0.09→0.93), StringFlav meson selection (the
always-drawn spin → pseudoscalars **and** vectors, η/η′ suppression, uds mixing), StringPT
pT (enhancedFraction draw kept), constituent-mass stop test, light-cone longitudinal
kinematics, and an exact 2-body finalTwo with refragment-on-failure.

Correctness: **4-momentum conservation 5.7e-14**, **on-shellness 1.6e-12**, GPU re-runs
bit-identical, **100% GPU≡CPU** (same hadron count per event over 20k), **0 failed events**,
2.7 M strings/s.

Physics (vs Pythia `forceHadronLevel` on the identical single u-ū string,
`multiplicity_pythia.cc`, decays+baryons off both sides): primary multiplicity
**12.64 vs 12.15 (4.0%)**, charged **6.98 vs 6.68 (4.5%)**. The GPU being slightly high is
consistent with the documented residuals (pole masses vs Pythia's Breit-Wigner vectors;
simplified finalTwo). Honest scope: pseudoscalar+vector mesons, uds, pole masses, single
straight string (no gluon kinks yet), no baryons/decays. Not aware of a prior *algorithmic*
(non-ML) GPU Lund port.

## Reproducing
```bash
cd cuPythia/pipeline && make check                 # all GPU stages, gated
# Pythia thrust reference (needs the vendored Pythia built + PYTHIA8DATA):
export PYTHIA8DATA=$PWD/../../pythia8317/share/Pythia8/xmldoc
g++ -O2 thrust_pythia.cc -o thrust_pythia $(../../pythia8317/bin/pythia8-config --cxxflags --libs)
./shower_fsr 200000 && ./thrust_pythia 100000 mecoff && bash compare_thrust.sh thrust_pythia_mecoff.dat
```
