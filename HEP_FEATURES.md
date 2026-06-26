# What cuPythia adds beyond stock Pythia — mapped to community needs

Framed around the priorities of the **HEP Software Foundation (HSF) Physics Event
Generator working group** and the **HEP-CCE "generators on accelerators"** effort.
cuPythia is a prototype; this is an honest map of its design to real needs, with
the gaps stated plainly.

## Why this matters

The HL-LHC era needs roughly an order of magnitude more simulated events than
today, and Monte-Carlo event generation is a recognized cost driver of the WLCG
compute budget. The community wants generators that:

- **run on GPUs** — the hardware new HPC centers actually field (Perlmutter,
  Frontier, Aurora, LUMI), where most FLOPs now live;
- **scale across clusters** — not just one node;
- **produce bit-reproducible events** for distributed (GRID) production and
  debugging.

Stock Pythia 8 is a serial C++ generator whose RNG (RANMAR, a sequential
lagged-Fibonacci-with-carry state machine) is inherently single-stream and stateful.

## What cuPythia demonstrates that stock Pythia lacks

| # | capability | stock Pythia | cuPythia |
|---|---|---|---|
| 1 | GPU matrix-element / MC evaluation | CPU only | kernels 01–03 (validated) |
| 2 | Multi-GPU + multi-node scaling | independent processes only (`PythiaParallel`) | kernel 04 (shards + MPI), near-linear |
| 3 | O(1) per-event reproducible RNG | sequential RANMAR (replay or state-restore) | kernel 05 (counter-based), bit-exact |

**(1) GPU acceleration.** Stock Pythia has no device backend. cuPythia ports the
data-parallel pieces (a real Pythia matrix element, `Sigma2gg2gg`, and the MC
integration around it) and measures them honestly — Amdahl caps end-to-end, and
on a consumer GPU the FP64 division rate caps the ME kernel; both are documented.

**(2) Cluster scaling.** Stock Pythia parallelizes only by running independent
instances. cuPythia's kernel 04 shards the MC across GPUs and nodes (MPI
`Allreduce`) with disjoint RNG substreams — embarrassingly parallel, so throughput
scales ~linearly with the number of GPUs. This is the one place a cluster gives
near-ideal speedup (it does **not** speed up a single event's sequential chain).

**(3) Reproducibility (kernel 05).** With a counter-based RNG, event *N*'s
randomness is a pure function of `(seed, N)`. You can regenerate **any single
event independently, on any node, in O(1)** — no replay, no checkpoint. RANMAR
cannot: reaching event *N* requires *N* advances or serializing/restoring state.
Demonstrated bit-exact: 100k events regenerated out of order and an 8-way node
partition both reassemble with zero difference. This is what GRID production and
debugging want — re-run only the one failed event of a million-event job, and get
bit-identical results across heterogeneous nodes. (This is the established HPC
pattern — cf. Random123 / Philox counter-based RNGs.)

## Honest gaps (what a real production-grade tool still needs)

- **Standard I/O:** HepMC3 output, LHE (Les Houches Event) read/write, and
  Rivet + ROOT interfaces — need external libraries; not added here.
- **Full physics on device:** only `gg→gg` ME plus toy MC are ported; the parton
  shower, hadronization, decays, PDFs, and the other `Sigma` processes remain
  CPU/sequential. The `AUDIT.md` study ranks which are worth porting next.
- **Portability layer:** real deployment wants one source compiling to CUDA + HIP
  + SYCL (Kokkos / Alpaka / SYCL) for Frontier (AMD) and Aurora (Intel), not raw
  CUDA. The kernels are written to make that port mechanical, but it isn't done.
- **Physics validation:** against Rivet analyses / data, not only analytic cross
  sections.

## References (direction, not endorsements)

- HSF Physics Event Generator WG — generator computing costs toward HL-LHC.
- HEP-CCE / generators-on-accelerators — `madgraph4gpu`, MadFlow.
- Random123 (Salmon et al.) — counter-based RNG, the reproducibility pattern used
  in kernel 05.
