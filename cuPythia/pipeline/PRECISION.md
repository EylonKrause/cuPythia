# cuPythia vs Pythia — precision, and corrections toward theory

## Who is more precise?

**Regular Pythia 8 is more physically precise, by a wide margin.** It is NLO-matched,
matrix-element–corrected, full-flavour (baryons, excited mesons, decays), uses 2-loop +
CMW α_s, and is tuned to LEP/LHC data. cuPythia is a from-scratch GPU port with
*deliberately documented* simplifications:

| aspect | Pythia 8 | cuPythia (this pipeline) |
|---|---|---|
| shower | NLL-ish, ME-corrected, g→qqbar, CMW α_s | **leading-log**, q→qg & g→gg (+ g→qqbar `-DGLUON_SPLIT`), 1-loop α_s (threshold-matched) |
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
3. **g→qq̄ shower splitting — DONE & VALIDATED** (`-DGLUON_SPLIT`, `shower_inc.cuh` +
   `hadronize_mr.cu`). Adds the missing final-state gluon splitting as a second trial channel for
   gluon radiators: the plain DGLAP kernel **T_R[z²+(1−z)²]** (T_R=½, flat-z), summed over flavours
   with Pythia's quark masses (u,d=0.33, s=0.5, c=1.5, b=4.8) and the **THRESHM2=4.004** pair
   threshold + β_Q. Three design points, each verified:
   - **RNG phase-lock:** a fixed 5 draws/trial-iteration for *every* dipole end (the channel- and
     flavour-draws are taken even for quark radiators, where they go unused) so control flow never
     depends on parton flavour → **100% GPU≡CPU**. Flag OFF draws 3 → **byte-identical** to the
     committed LL shower (mult 12.719, ⟨1−T⟩ 0.0690, and N_gqq=**0** exactly: the clean A/B gate).
   - **No double-count:** a gluon sits in two dipoles and each end generates the trial, so the
     single g→qqbar conversion is **shared ½ per end** (coefficient T_R·½ = Pythia's `wtPSqqbar`).
     Without the ½ the rate came out 2.14× high.
   - **Colour-string fork:** g→qqbar cuts the chain into two strings (q…q̄′)(q′…q̄); `findStrings`
     slices at each q̄-then-q boundary and the hadronizer fragments each sub-string independently.
   **Validation** (200k, vs Pythia `weightGluonToQuark=1` — the same plain kernel — MEcorr off,
   reference counter `thrust_pythia g2q1`): secondary-pair rate **N_gqq = 0.578 vs 0.566 (+2.2%)**,
   flavour-resolved uds +2.2%, **c +0.3%**, b +6.8%; thrust non-regression (0.0702 vs 0.0700, <1%);
   hadron-level **4-momentum conservation exact (1.43e-10)** across the forked multi-string events,
   reproducible, drop rate 3.0% (lower — forks simplify the region geometry). Hadron multiplicity
   shifts −0.4% (20.92 vs 21.01): a g→qqbar replaces a multiplicity-enhancing gluon kink with a
   quark string-break, so the small decrease is physical, not a regression. **Honest caveat:** the
   GPU shower is massless and uses Pythia **option 1**; Pythia's *default* is option 4 (massive
   zCosThe reshape + pow3(1−m²/m²_dip) damping → a lower 0.5233), which needs the massive-recoil
   `doKin` path (future). So this validates the LL kernel, not the option-4 phenomenology.
4. **Breit–Wigner vector masses — DONE & VALIDATED** (`bw_inc.cuh`, `-DUSE_BW`). Vector
   mesons (ρ, K*, ω, φ) get a relativistic Breit-Wigner mass (Lorentzian in s, truncated ±3Γ,
   exact inverse-CDF, one RNG draw, on-shell check uses the *sampled* mass). On the straight-
   string hadronizer the ρ spectrum is correctly broadened: **mean 0.779 (pole 0.775), RMS
   0.132 ≈ Γ=0.149**, with 4-momentum conservation still exact (5.7e-14) and multiplicity
   unchanged. Applied to **both** hadronizers (straight-string and multi-region; the
   multi-region also conserves exactly 4.95e-11 with BW). Default (pole) build is unaffected.
5. **Real PDF — DONE & VALIDATED.** `genpdf.cc` fills the `pdf.cuh` log-x/log-Q² grid from
   Pythia 8.317's **real proton gluon PDF** (the default NNPDF set) and the device interpolator
   reproduces it to **3.24e-3** in the σ-support region (xf_g(0.01,100)=7.9801 vs 7.9796) —
   the real-PDF analog of the toy validation. The grid (`real_pdf.grid`) drops into the same
   arrays for a **physical** hadronic σ instead of the illustrative toy.
6. **2-loop α_s** and **NLO matching** (POWHEG/MC@NLO) — the precision frontier.

Items 1–5 are DONE & VALIDATED (six corrections total, all opt-in via `-D` flags so the
committed LL/LO results stay byte-stable); item 6 is the remaining research frontier. The g→qqbar
addition (item 3) was the last *correctness* gap (flavour/string topology), not just a precision
knob — it took a dedicated research+design pass to get the sharing, masses and colour-fork right.

## Bottom line
Pythia is the more physically precise generator. cuPythia's contribution is a
**reproducible, GPU-native** generator validated at the LL/LO level (~2–7% vs Pythia across six
independent observables: thrust, multiplicity, g→qqbar rate, …), with a clear, prioritized path to
higher theoretical precision — and an honest record (the naive ME correction over-corrects a dipole
shower; the g→qqbar rate matches Pythia's *plain-kernel* option, not its damped default).
