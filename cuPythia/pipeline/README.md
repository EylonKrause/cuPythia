# cuPythia/pipeline — the device-resident generator (in progress)

Kernels 00–14 are isolated, validated *demonstrations*. **This** directory is the
real thing being built: a fused, device-resident event generator where an event
**never leaves the GPU between stages**.

## Data plane — `event.cuh`
A Structure-of-Arrays `DeviceEvents` record: per-particle four-momenta + mass, PDG
id, status, colour/anticolour tags, mother indices; per-event seed/weight/scale.
Every stage reads and writes particles here, in device memory. Per-event RNG is
counter-based, so any single event is **O(1)-reproducible** on any node.

## Stages (built incrementally, each validated; correctness gated by `make`/tests)
- [x] **stage 0 — `build_events.cu`**: populate gg→gg hard-process events into the
  record on-GPU. Validated: exact 4-momentum conservation (**0 imbalance**), record
  integrity (all events well-formed), cross-section sanity vs quadrature.
- [x] **stage 1 — `pdf_xsec.cu`**: device PDF evaluator (`pdf.cuh`: log-x/log-Q²
  grid, **log(xf) bilinear interpolation** with edge freezing) convolved to a real
  **hadronic** gg→gg σ (13 TeV pp, pT-hat>50 GeV). Validated: interp fidelity at the
  σ level **1.7e-4**, GPU-vs-CPU on identical samples **1.4e-11** (determinism). A
  real LHAPDF `.dat` grid plugs into the same arrays + interpolator with no kernel change.
- [x] **stage 2 — `reweight.cu`**: N scale-variation weights per event in one pass,
  **bit-identical (max|diff|=0) to N independent pinned re-runs** (the counter-RNG
  advantage); physical ±25% LO scale band. μ_F/PDF variations await stage 1.
- [x] **orchestrator — `generate.cu`**: the device-resident parton-level generator —
  build → reweight → unweight → CUB-compact, **all on one record, no host round-trip**
  (the gap Pepper/madgraph4gpu concede). σ vs quadrature, scale band, CUB count exact.
- [x] **stage 3 — `shower_fsr.cu`**: a physical final-state (timelike) **dipole shower**,
  one event per GPU thread (GAPS pattern), with Pythia `SimpleTimeShower` splitting kernels
  $(1+z^2)/2$, $(1+z^3)/2$, running-α_s trial generation, z-sampling, and exact local-dipole
  **recoil kinematics**. Validated on e⁺e⁻→Z→qq̄: 4-momentum conservation **1.8e-9**,
  on-shellness **1.5e-12**, GPU re-runs bit-identical, and the **control flow is 100%
  bit-identical to an independent CPU port** (momenta agree to 1.4e-12, GPU/CPU IEEE limit).
  **Physics validated vs Pythia's own SimpleTimeShower** (`thrust_pythia.cc`, identical
  setup; `compare_thrust.sh`): the **thrust** distribution agrees to **4.0% (mean over all
  bins)** vs Pythia with ME-corrections off — the apples-to-apples LL comparison — with
  ⟨1−T⟩ 0.0690 vs 0.0698. The 23.8% spread vs *default* Pythia is precisely its ME
  corrections, a higher-order effect no LL shower (GAPS included) carries. Now includes
  **flavour-threshold α_s** (n_f=5,4,3 across m_b,m_c; tightened the fit 5.2%→4.0%). Scope:
  FSR-only, massless, q→qg & g→gg; the ~4% residual is g→qq̄ + exact-Λ α_s + MC statistics,
  and ME corrections are the LL→NLL next step. (NB: `TimeShower:nGluonToQuark=0` hangs
  Pythia 8.317, so the reference keeps g→qq̄ on.)
- [ ] stage 4 — hadronization on device (feasibility-gated) + decays
- [ ] stage 5 — standard I/O (spec-valid LHE / HepMC3 / Rivet smoke test)

Design grounded in a web-research pass on GPU-generator architecture
(Pepper, madgraph4gpu, GPU showers, device PDFs) → see `../../ARCHITECTURE.md`.

## Build
```bash
nvcc -O3 -std=c++17 -arch=sm_120 -o build_events build_events.cu && ./build_events
```
