# cuPythia vs Pythia — precision, and corrections toward theory

## Who is more precise?

**Regular Pythia 8 is more physically precise, by a wide margin.** It is NLO-matched,
matrix-element–corrected, full-flavour (baryons, excited mesons, decays), uses 2-loop +
CMW α_s, and is tuned to LEP/LHC data. cuPythia is a from-scratch GPU port with
*deliberately documented* simplifications:

| aspect | Pythia 8 | cuPythia (this pipeline) |
|---|---|---|
| shower | NLL-ish, ME-corrected, g→qqbar, CMW α_s | **leading-log**, no ME corr, q→qg & g→gg, 1-loop α_s (threshold-matched) |
| hadronization | full Lund: BW masses, all multiplets, baryons, decays | pseudoscalar+vector, **pole masses**, uds, no baryons/decays |
| hard process | NLO/merged, LHAPDF | LO gg→gg, toy PDF (real LHAPDF drops in) |
| tuning | decades of data tunes | untuned (Pythia defaults) |

So on *matching theory and data*, Pythia wins. cuPythia's measured agreements with
Pythia are honest LL/LO-level: shower thrust **4.0%** (vs Pythia-LL), hadron multiplicity
**~7%**, end-to-end LEP n_ch≈21.

### Where cuPythia is *more* precise — methodology, not physics accuracy
- **Exact reproducibility & bit-identical reweighting.** The counter-based RNG makes any
  event O(1)-reproducible on any node, and N scale-variation weights are **bit-identical
  (max|diff|=0)** to N independent pinned re-runs (`reweight.cu`). Pythia's sequential
  RANMAR cannot do per-event O(1) reproducibility; its on-the-fly reweighting is
  statistically equivalent, not bit-exact.
- **Determinism of the GPU port.** Control flow is bit-identical GPU-vs-CPU; momenta agree
  to the IEEE transcendental limit (~1e-12).
- **Throughput → statistical precision** on the accelerated stages (Amdahl- and, on this
  consumer GPU, FP64-bounded).

These are real precision advantages, but of *method* (reproducibility, statistics), not of
*physics accuracy*.

## Corrections to bring cuPythia closer to theory (prioritized)

1. **First-emission ME correction — DONE & VALIDATED (highest impact).** The 23.8% thrust
   gap vs *default* Pythia was the ME correction. Two routes were tried; the honest record:
   - **Naive reweight R=(x1²+x2²)/2 — FAILED (reverted):** it assumes the *old mass-ordered*
     shower's eikonal density; on this pT-ordered dipole shower (which already applies the
     (1+z²)/2 kernel) it *double-suppresses* (⟨1−T⟩ 0.069→0.0517).
   - **Direct ME generation (POWHEG-style) — WORKS** (`-DME_FIRST`, `sampleFirstEmission` in
     `shower_inc.cuh`): the hardest emission is generated *directly* from the exact Dalitz
     density (x1²+x2²)/((1-x1)(1-x2)) via a **pT-ordered Sudakov veto** (t=pT²/Q², y=½ln((1-x1)/(1-x2)),
     rate ∝ α_s(t)(x1²+x2²) dt/t dy), then the LL shower runs below pT_first. No (1+z²) kernel
     on the first emission → no double-count. (One debug: the Sudakov `Iover` must carry the
     C_F/2π prefactor — without it the rate was ~4.7× too high, ⟨1−T⟩→0.133.) **Result:** thrust
     vs *default* Pythia (ME corrections ON) **23.8% → 7.9%**, ⟨1−T⟩ **0.0641 vs 0.0635**;
     conservation exact, 100% GPU≡CPU. The first emission now reproduces the exact O(α_s) ME.
2. **CMW (Catani–Marchesini–Webber) α_s — DONE & VALIDATED.** Soft-gluon NLL coherence:
   α_s → α_s(1 + α_s K/2π), K = C_A(67/18 − π²/6) − 5n_f/9. Added behind `-DUSE_CMW`
   (default off, so the LL results above are unchanged). It correctly increases soft
   radiation: ⟨1−T⟩ 0.069→**0.0768**, matching Pythia with `alphaSuseCMW=on` (0.070→0.0791)
   to **5.4%** — confirming the NLL rescaling is implemented correctly. (It is a shower-soft
   refinement; a full-event tune would re-fit α_s(M_Z) since CMW raises the multiplicity.)
3. **g→qq̄ in the shower** — restores the missing splitting (correct flavour/multiplicity);
   needs multi-string bookkeeping (a gluon→qqbar splits the colour chain).
4. **Breit–Wigner resonance masses** in hadronization (ρ, K*, ω, φ, …) — correct mass
   spectrum instead of pole masses.
5. **Real PDF — DONE & VALIDATED.** `genpdf.cc` fills the `pdf.cuh` log-x/log-Q² grid from
   Pythia 8.317's **real proton gluon PDF** (the default NNPDF set) and the device interpolator
   reproduces it to **3.24e-3** in the σ-support region (xf_g(0.01,100)=7.9801 vs 7.9796) —
   the real-PDF analog of the toy validation. The grid (`real_pdf.grid`) drops into the same
   arrays for a **physical** hadronic σ instead of the illustrative toy.
6. **2-loop α_s** and **NLO matching** (POWHEG/MC@NLO) — the precision frontier.

Items 2, 4, 5 are clean and tractable; item 1 (done correctly) is the biggest physics win;
items 3, 6 are research-scale.

## Bottom line
Pythia is the more physically precise generator. cuPythia's contribution is a
**reproducible, GPU-native** generator validated at the LL/LO level (~4–7% vs Pythia), with
a clear, prioritized path to higher theoretical precision — and an honest record that the
naive ME correction over-corrects a dipole shower and must be derived for *this* shower.
