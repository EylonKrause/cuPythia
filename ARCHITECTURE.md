# cuPythia ARCHITECTURE — building the device-resident GPU event generator

The master design for turning the 15 validated demo kernels into a real, fused,
**device-resident** event generator. Grounded in a web-research pass on the state
of the art (sources in §8). Honest about what is achievable on one consumer GPU
and what is genuinely research-scale.

## 1. Goal & honest scope
A generator where **an event never leaves the GPU between stages**: hard process →
PDF/luminosity → multi-weight → parton shower → string construction →
hadronization → decays → event I/O, all reading/writing one resident event record.
- **Achievable now:** a fully device-resident **parton-level** generator (hard
  process + PDFs + reweighting + a real FSR shower) with standard I/O. This already
  fills a gap Pepper and madgraph4gpu concede — they stop at parton level and hand
  events to **CPU** Pythia/Sherpa to shower.
- **Research-scale (scoped as such, not promised in a v1):** GPU hadronization (no
  published implementation exists) and NLL/multi-jet-merged accuracy.

## 2. Architecture decision: STAGED kernels over one resident record (NOT a megakernel)
Three sourced facts drive this:
1. **GAPS** (Sule & Seymour, arXiv:2403.08692 / 2511.19633) — the *only* published,
   validated GPU parton shower — keeps the **sequential Sudakov veto unchanged** and
   parallelizes **one event per thread**, splitting each emission cycle into separate
   lock-step kernels (`SelectWinner`/`CheckCutoff`/`AcceptOrVeto`/`GenerateEmission`)
   with a sync per cycle. We adopt this pattern verbatim for the shower.
2. **madgraph4gpu / CUDACPP** (arXiv:2510.05392, Oct 2025) **split** their monolithic
   ME kernel into many small kernels to cut register pressure and raise occupancy.
   On a small consumer card (RTX 5050, ~20 SMs, weak FP64) register pressure/occupancy
   are the binding constraints — a 15-stage megakernel would spill and serialize.
   **Split kernels win.**
3. **Residency is the actual differentiator.** `pipeline/event.cuh` already encodes it:
   one `DeviceEvents` SoA buffer; device→host copies happen ONLY at final I/O and tiny
   scalar reductions (cross-section sums, validation residuals).

**Design:** a host-side **orchestrator** issues a sequence of kernel launches over one
persistent `DeviceEvents`. Iterative stages (shower, hadronization) use the GAPS
split-kernel + per-cycle-sync + **CUB stream-compaction** sub-loop (we already have
`cub::DeviceSelect` validated in kernel 10). Raw CUDA first; portability layer
(Kokkos/Alpaka) deferred. One CUDA stream for the prototype.

## 3. The data plane — `pipeline/event.cuh`
Extend the proven SoA record (flat index `e*maxPart + p`; build_events.cu already gets
**0 four-momentum imbalance**). Current: `px,py,pz,e,m` (FP64); `pdg,status,col,acol,
mo1,mo2` (int); per-event `nPart, seed, weight, scale`. **To add:**
- per-particle: `d1,d2` (daughter links — a decay tree needs both mother and daughter).
- per-event: **weight VECTOR** `weight[nEvents*NVAR]` (slot 0 = nominal, MCnet naming
  convention arXiv:2203.08230; FP32 unless precision needs FP64 — weights×events is the
  memory ceiling); `x1,x2` (parton momentum fractions — mandatory for μ_F/PDF reweight
  and ISR); `flavA,flavB` (incoming flavours); `active`/`endShower` flag (GAPS idle flag;
  the field stream-compaction reads); a small per-event string-endpoint table.

## 4. Stage pipeline (R=reads, W=writes the record)
- **Stage 0 — hard process / build** *(DONE: build_events.cu)*. One thread/event:
  sample, fill the gg→gg partons. W: momenta, status, nPart, seed, scale (+ x1,x2,flav).
- **Stage 1 — PDF convolution + luminosity** *(NEW, biggest missing piece)*. Device
  xfx(x,Q²) grid (constant/texture, bilinear→bicubic, all 11 flavours/call, GAPS-v2
  ~70× vs 1 CPU core). **Low-x/Q² freezing from day one** (GAPS's open bug without it).
- **Stage 2 — multi-weight reweight** *(μ_R DONE: reweight.cu; bit-identical to N pinned
  re-runs)*. Add μ_F/PDF-member variations using the Stage-1 xfx cache.
- **Stage 3 — parton shower** *(toy DONE: kernel 14; real version is the headline build)*.
  GAPS split-kernel loop: SelectWinner→CheckCutoff→AcceptOrVeto→GenerateEmission, recoil
  + colour reassignment + variable P(z) + running αs, appending partons. Validate
  event-by-event bit-identity vs a CPU port on the SAME counter-RNG, THEN Rivet
  histograms (Durham jet rates, thrust, heavy-jet-mass) vs Pythia SimpleTimeShower.
- **Stage 4 — string construction** from col/acol (prerequisite for hadronization).
- **Stage 5 — hadronization** *(research-scale)*. Track A: batch-over-strings native
  Lund + action-based per-step kernel + CUB compaction (exploiting causally-disconnected
  breaks); start pion-restricted. Track B: an MLHad-style normalizing-flow surrogate.
  Validate both vs CPU Pythia 8.317 multiplicity / z / pT spectra. **Not promised in v1.**
- **Stage 6 — decays** (table-driven; easy after hadronization).
- **Stage 7 — unweighting + I/O** (CUB compaction → spec-valid LHE / HepMC3 weight vector).

## 5. Build order (validatable, incremental — each step shippable)
1. **Lock the data plane:** extend `event.cuh` (weight vector, x1,x2, flav, active, d1,d2,
   string table); re-run build_events, assert still 0 imbalance + same σ.
2. **Orchestrator:** one driver, `DeviceEvents` allocated once, Stage 0→2→7 with **no host
   round-trip** — the Pepper-gap parton-level residency win, achievable now.
3. **Device PDF evaluator** (Stage 1): one LHAPDF grid → constant/texture, bilinear then
   bicubic, low-x/Q² freezing; validate xfx vs CPU LHAPDF <1e-4 and a hadronic σ vs Pythia.
4. **Full (μ_R,μ_F,PDF) reweight** on the xfx cache; HepMC3 named weights + Rivet smoke test.
5. **Real FSR dipole shower** (GAPS pattern): bit-identity vs a CPU port first, then Rivet
   observables vs Pythia. **The headline novel result.**
6. **Stream-compaction in the shower** + benchmark the break-even batch size on the RTX 5050
   (GAPS saw it HURT below ~10× the FP64 core count); then massive partons, then ISR.
7. **String construction** from colour tags; validate colour-consistent + momentum-complete.
8. **Hadronization** (two tracks, honest research scope) — pion-restricted Lund first.
9. **Decays + I/O hardening** (wMax estimation, surplus-event tracking, efficiency report).
10. **Counter-RNG discipline throughout** — every draw a pure function of
    `(event/string id, stream id, sub-counter)` with disjoint key namespaces, preserving
    the bit-identical CPU-equivalence test as the cornerstone of validation.

## 6. Hardest parts (ranked) & honest risks
1. **Hadronization** — unsolved in the literature; research-scale (multi-month+), ship as
   a prototype reproducing Pythia distributions, NOT a production fragmenter.
2. **Real shower** — divergence is *mitigated, not solved*; small-card occupancy may make
   GAPS-style compaction net-negative at feasible batch sizes — must be measured.
3. **Device PDFs with correct extrapolation** — on the critical path; build freezing in from
   day one.
4. **Unweighting efficiency** — needs reliable wMax + surplus bookkeeping, no silent drops.
5. **Warp-divergence tax** — cross-cutting; benchmark compaction break-even per stage.
6. **SoA refactor / orchestration glue** — without breaking the 0-imbalance / bit-identical
   guarantees.

Other risks: memory is the binding constraint (FP32 weights, size maxPart from the 99.9th
percentile, not worst case); CPU-equivalence validation breaks if any stage adds hidden
mutable state or sampling-steering non-deterministic reductions; accuracy ceiling is
leading-log dipole (GAPS already shows worse subleading-jet agreement); **avoid scope creep**
into a megakernel or premature portability layer.

## 7. Validation philosophy
Two complementary gates at every stage: (1) **bit-identity** vs a CPU reference consuming the
same counter-RNG stream (the counter-based RNG makes this trivial — GAPS's decisive check);
(2) **physical distributions** vs analytic results and, where external libs allow, Pythia /
Rivet. Reductions that only *report* (σ sums) are fine; reductions that *steer sampling* must
stay deterministic.

## 8. Sources
GAPS arXiv:2403.08692, 2511.19633 · Pepper arXiv:2311.06198 · madgraph4gpu/CUDACPP
arXiv:2106.12631, 2510.05392 · MadtRex arXiv:2510.05100 · MCnet weight-naming arXiv:2203.08230 ·
MLHad (ML hadronization) line · HSF/HEP-CCE generator-on-accelerator reviews.
