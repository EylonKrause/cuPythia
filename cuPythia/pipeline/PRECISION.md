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

1. **First-emission ME correction (highest impact).** The 23.8% thrust gap vs *default*
   Pythia is exactly the ME correction: the LL shower over-radiates the hard, wide-angle
   first emission; the exact O(α_s) γ*/Z→qqg matrix element ∝ (x1²+x2²)/((1-x1)(1-x2))
   tames it. **Attempted and reverted (honest finding):** the textbook Bengtsson–Sjöstrand
   reweight R=(x1²+x2²)/2 assumes the *old mass-ordered* Pythia shower's first-emission
   density. Applied to this **pT-ordered dipole** shower (whose local recoil already
   captures part of the ME), it *double-suppresses* — ⟨1−T⟩ overshot 0.069→0.0517 (target
   0.0635). A correct MEC must use **this shower's own** first-emission density to form the
   ME/PS ratio (a derivation, not a one-liner). This is the right #1 next step, done right.
2. **CMW (Catani–Marchesini–Webber) α_s** — soft-gluon NLL coherence:
   α_s → α_s(1 + α_s K/2π), K = C_A(67/18 − π²/6) − 5n_f/9. A clean, ~1-line NLL
   improvement (validate vs Pythia with `TimeShower:alphaSuseCMW=on`).
3. **g→qq̄ in the shower** — restores the missing splitting (correct flavour/multiplicity);
   needs multi-string bookkeeping (a gluon→qqbar splits the colour chain).
4. **Breit–Wigner resonance masses** in hadronization (ρ, K*, ω, φ, …) — correct mass
   spectrum instead of pole masses.
5. **Real LHAPDF grid** instead of the toy PDF (the `pdf.cuh` log-x/log-Q² interpolator
   already accepts a real `.dat` unchanged) — correct hard cross sections.
6. **2-loop α_s** and **NLO matching** (POWHEG/MC@NLO) — the precision frontier.

Items 2, 4, 5 are clean and tractable; item 1 (done correctly) is the biggest physics win;
items 3, 6 are research-scale.

## Bottom line
Pythia is the more physically precise generator. cuPythia's contribution is a
**reproducible, GPU-native** generator validated at the LL/LO level (~4–7% vs Pythia), with
a clear, prioritized path to higher theoretical precision — and an honest record that the
naive ME correction over-corrects a dipole shower and must be derived for *this* shower.
