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
- [ ] stage 1 — PDF convolution + process selection (a real hadronic σ, not partonic)
- [ ] stage 2 — multi-weight systematic reweighting (scale/αs/PDF, one resident pass)
- [ ] stage 3 — physical parton shower (recoil + colour + variable P(z)) into the record
- [ ] stage 4 — hadronization on device (feasibility-gated) + decays
- [ ] stage 5 — standard I/O (spec-valid LHE / HepMC3 / Rivet smoke test)

Design grounded in a web-research pass on GPU-generator architecture
(Pepper, madgraph4gpu, GPU showers, device PDFs) → see `../../ARCHITECTURE.md`.

## Build
```bash
nvcc -O3 -std=c++17 -arch=sm_120 -o build_events build_events.cu && ./build_events
```
