# 07 — GPU unweighting efficiency + Les Houches (.lhe) output

Two community-wanted capabilities, grounded in the literature:

1. **Unweighting efficiency** — the headline metric of GPU parton-level generators
   (Pepper, arXiv:2311.06198; madgraph4gpu): `η = <w>/w_max` = the fraction of
   generated (weighted) events kept by acceptance-rejection. It sets the real cost:
   you generate `1/η` events for every unweighted event a detector sim consumes.
2. **Standard Les Houches event output** — the parton-level interchange format
   Pythia / Herwig / the experiments read — written with **no external library**.

## Result — RTX 5050, 1.05×10⁹ trials, gg→gg, |cosθ| < 0.9
- unweighting efficiency **η = 10.0%** (cross-check `<w>/w_max` = 10.0% — exact)
- σ (unweighted MC) matches Simpson quadrature, relerr 3.5e-5
- **1000 unweighted events** written to `events.lhe` (valid LHE 3.0: proton beams,
  `IDWTUP=3`, back-to-back gluons, colour flow, energy-momentum conserved)
- `VALIDATION: PASS`

η = 10% illustrates the generator-cost problem head-on: 90% of generated events
are thrown away. Raising η is exactly what importance sampling / MadNIS
(arXiv:2212.06172) targets — the motivation for the next kernel (VEGAS-on-GPU).

## Build / run
```bash
nvcc -O3 -arch=sm_120 -o unweight_lhe unweight_lhe.cu
./unweight_lhe [trialsPerThread=4000]   # writes events.lhe
```
