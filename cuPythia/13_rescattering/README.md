# 13 — O(N²) hadronic rescattering (heavy-ion)

The Pythia team flagged hadronic rescattering as the prime GPU target: in heavy-ion
events thousands of hadrons are produced and collision-finding costs grow as the
**square** of the multiplicity — an embarrassingly-parallel all-pairs screen. Here
N hadrons are generated at freeze-out and interacting pairs (closest approach within
the interaction radius `d = √(σ/π)`) are counted on the GPU and validated **exactly**
against a CPU reference.

## Result — RTX 5050, N=20,000 hadrons
- 2.0×10⁸ pairs screened, interaction radius 1.128 fm
- interacting pairs: **GPU 399,366 == CPU 399,366** (exact)
- GPU 1.0×10¹⁰ pairs/s, 7.8× over CPU (grows with N)
- `VALIDATION: PASS`

This is the geometric collision-**finding** core. Processing collisions in **time
order** is intrinsically sequential (each rewrites momenta) — the genuinely hard
research piece, and the next step here.

```bash
nvcc -O3 -std=c++17 -arch=sm_120 -o rescatter rescatter.cu && ./rescatter [N]
```
