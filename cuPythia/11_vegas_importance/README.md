# 11 — VEGAS adaptive importance sampling

Unweighting efficiency `η = ⟨w⟩/w_max` is the generator cost driver: with **uniform**
sampling (kernel 07) gg→gg gives η ≈ 10% — 90% of generated events thrown away.
**VEGAS** (Lepage 1978) adapts a piecewise sampling grid to the integrand so the MC
weight `w = f/p` flattens, pushing η up. It is the classical precursor to neural
importance sampling (MadNIS, arXiv:2212.06172).

Each iteration: a GPU kernel samples from the current grid and accumulates the
integral, the max weight, and per-bin importance (shared-memory binning); the host
rebins the grid for equal importance. 128 bins, 12 iterations.

## Result — RTX 5050, gg→gg
- unweighting efficiency: **uniform 10.0% → VEGAS 76.0% (7.6× better)**
- integral matches Simpson quadrature to **relerr 1.2e-5**
- `VALIDATION: PASS`

7.6× fewer wasted events for the same physics — directly the HL-LHC generator-cost
lever. Next step on this axis: a learned sampler (normalizing flow), the MadNIS
direction.

## Build / run
```bash
nvcc -O3 -arch=sm_120 -o vegas vegas.cu
./vegas [samplesPerThread=2000]
```
