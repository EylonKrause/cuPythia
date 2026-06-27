# cuPythia vs Pythia — precision, and corrections toward theory

## Who is more precise?

**Regular Pythia 8 is more physically precise, by a wide margin.** It is NLO-matched,
matrix-element–corrected, full-flavour (baryons, excited mesons, decays), uses 2-loop +
CMW α_s, and is tuned to LEP/LHC data. cuPythia is a from-scratch GPU port with
*deliberately documented* simplifications:

| aspect | Pythia 8 | cuPythia (this pipeline) |
|---|---|---|
| shower | NLL-ish, ME-corrected, g→qqbar, CMW α_s | **leading-log**, q→qg & g→gg (+ g→qqbar `-DGLUON_SPLIT`), 1-loop α_s (threshold-matched) |
| hadronization | full Lund: BW masses, all multiplets, baryons, decays | pseudoscalar+vector, pole masses (+BW `-DUSE_BW`), uds, **no baryons** (+ hadron decays `-DDECAYS`) |
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
   GPU shower is massless and uses Pythia **option 1** by default. Pythia's *default* weight is
   **option 4** (zCosThe reshape + (1+m²Rat)(1−m²Rat)² high-mass damping); `-DG2QQ_WEIGHT4` applies
   that full weight (with a veto oversample, since zCosThe makes it exceed 1) and reproduces the
   **default rate 0.5233 to +1.9%** — so cuPythia matches *both* of Pythia's kernel choices. The
   produced kinematics stay massless (the zCosThe massive-recoil construction is a future `doKin`
   path), so `-DG2QQ_WEIGHT4` is a **rate-level** match to the default, not its exact z-distribution.
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
6. **2-loop α_s — DONE & VALIDATED, with an instructive negative result** (`-DAS_2LOOP`,
   `shower_inc.cuh`; validator `as2loop_validate.cc`). Implements Pythia's `AlphaStrong` order=2
   (RPP 2006 eq. 9.5): α_s = 12π/(b₀L)·(1 − b₁ lnL/L), L=ln(μ²/Λ_nf²), with Λ_{5,4,3} iteratively
   matched to α_s(M_Z)=0.1365 (mc=1.5, mb=4.8). The running is **bit-identical to Pythia** (max
   relErr **0.0e+00** at every perturbative scale; Λ₅,₄,₃ = 0.5530, 0.7229, 0.7565 GeV; α_s(M_Z²)
   recovered exactly) including the SAFETYMARGIN2 low-scale freeze. **The honest finding:** plugged
   into this LL shower it **over-radiates badly** — parton multiplicity 12.7 → **24.3** (≈2×),
   ⟨1−T⟩ 0.069 → 0.078 — because 2-loop α_s is far larger near the cutoff (4.35 vs 1.62 at pT²ₘᵢₙ)
   and an LL shower has no NLL terms to compensate. So 2-loop α_s is *correctly implemented* but
   *physically inappropriate at LL* — it concretely demonstrates **why showers use 1-loop(+CMW)**.
   Kept opt-in and off by default.
7. **Hadron decays — DONE & VALIDATED** (`-DDECAYS`, `decay_inc.cuh`; validator `decay_test.cu` in
   `make check`). The real-ALEPH-data (Rivet) comparison identified hadron decays as the **#1 gap**;
   this closes most of it. GPU recursive decays of the primary unstable mesons (ρ/K*/ω/φ/η/η'/K⁰ →
   π/K/γ/K_S/K_L) into the ALEPH particle-level stable set: a baked decay table (BRs from Pythia
   8.317, renormalized), isotropic **2-body** + flat-Dalitz **3-body** kinematics, a recursion-free
   LIFO stack, a **fixed counter-RNG draw budget** (host≡CPU bit-identical) on a **separate stream**
   so the no-decay build stays byte-identical. **Validated:** per-decay 4-momentum closure **1.07e-14**,
   on-shell 1.45e-14, GPU≡CPU 100%, reproducible; in the pipeline charged multiplicity **11.3 → 18.69**
   (vs ALEPH 20.73, −9.8%) and the ALEPH **thrust χ²/ndf 39 → 10.3** (3.8× better), conservation still
   exact (1.23e-10). **Honest caveat:** the −9.8% residual is the floor from the *other* omissions
   (no baryons ~5–6%, heavy flavour only via g→qqbar, flat Dalitz, untuned Lund, no detector); π⁰/K_S
   kept stable per the Rivet particle-level convention. v1 = light-unstable-meson 2-/3-body flat phase
   space; D/B left stable.
8. **Baryon production — DONE & VALIDATED (flavour mechanism)** (`-DBARYONS`, in `hadronize_mr.cu`;
   validator `baryon_test.cu` in `make check`; decay rows in `decay_inc.cuh`). The Lund diquark
   mechanism, ported faithfully from Pythia 8.317 `FragmentationFlavZpT.cc` (`pick()` diquark-vs-quark
   + `combineId` SU(6) Clebsch–Gordan with octet/decuplet split and Λ/Σ disambiguation): a string
   break makes a diquark–antidiquark pair (prob `probQQtoQ`=0.081), forming p/n/Λ/Σ/Ξ/Ω/Δ/Σ*/Ξ* with
   the correct Pythia sign conventions; decuplet + Σ⁰ resonances decay. **Validated:** electric-charge
   and **baryon-number conservation 0 / 2,000,000 strings** (`baryon_test`); and in the **full chain**
   (ME+g→qqbar+decays+baryons, `hadronize_mr_full`) charge + baryon-number conserved **0 / 19334 events**,
   4-momentum exact (1.23e-10), with **p+p̄ 0.745, Λ+Λ̄ 0.214** per event (right order, ~60–70% of PDG).
   **The validator earned its keep — it caught two real bugs that exact 4-momentum conservation
   completely masked:** (a) a finalTwo *two-diquark* combine → baryon-number violation (fixed: finalTwo
   uses a single qq̄ pair), and (b) out-of-scope *charm/bottom baryons* from g→qqbar endpoints → charge
   violation (fixed: baryon breaks only off light uds endpoints; c/b stays D/B mesons). **Honest finding:**
   with *untuned* `probQQtoQ`/Lund params, baryons slightly *reduce* charged mult (18.69 → 17.30) and
   *worsen* the ALEPH thrust fit (χ²/ndf 10.3 → 13.1) — a conservation-correct mechanism can move the
   wrong way on data when its rate isn't tuned; closing to 20.73 needs a Lund tune + D/B decays, not
   more mechanisms. v1 is light (uds) baryons, octet+decuplet, no popcorn; charm/bottom baryons deferred.
   **Build note (key enabler):** the combined kernel was one monstrous inlined `__global__` that
   exhausted the compiler; marking the heavy device functions (`showerEvent`, `tryFragmentMR`, `kinHad`,
   `finalRegion`, `decayEvent`, `regionSetUp`, …) **`__noinline__`** splits it into tractable pieces
   with **bit-identical physics** (make check 12/12, byte-identical). The full chain then compiles
   (~10 min at -O2). The default -O3 build (flags off) is unaffected and byte-identical.
9. **NLO matching** (POWHEG/MC@NLO) — the genuine remaining frontier (a research program, not a
   flag): it requires the NLO virtual+real subtraction and a matching scheme, beyond this LL/LO
   port. Documented honestly as out of scope rather than approximated. (The first-emission ME
   correction, item 1, already supplies the O(α_s) real ME for the hardest emission.)

Items 1–8 are DONE & VALIDATED (all opt-in via `-D` flags so the committed LL/LO results stay
byte-stable); item 9 (NLO matching) is the genuine remaining frontier, documented honestly rather
than approximated. The g→qqbar addition (item 3) was the last *correctness* gap (flavour/string
topology), not just a precision knob — it took a dedicated research+design pass to get the sharing,
masses and colour-fork right. Item 6 (2-loop α_s) is an honest *negative* result: implemented
bit-identically to Pythia, it over-radiates at LL — showing exactly why the shower uses 1-loop.

## Bottom line
Pythia is the more physically precise generator. cuPythia's contribution is a
**reproducible, GPU-native** generator validated at the LL/LO level (~2–7% vs Pythia across six
independent observables: thrust, multiplicity, g→qqbar rate, …), with a clear, prioritized path to
higher theoretical precision — and an honest record (the naive ME correction over-corrects a dipole
shower; the g→qqbar rate matches Pythia's *plain-kernel* option, not its damped default).
