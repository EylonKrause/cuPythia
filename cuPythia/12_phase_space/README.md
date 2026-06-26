# 12 — 2→2 phase-space generation (full final-state kinematics)

Kernels 02/03 sampled only the scattering angle. Here each GPU thread generates
the full final-state **four-momenta** of a 2→2 event (what you'd write to an LHE /
HepMC record), for a massless (gg→gg) and a massive (m=1.5 GeV, charm-like) final
state.

## Result — RTX 5050, √s=100 GeV
| | max \|p conservation\| | max \|on-shell\| | max \|ŝ−s\| |
|---|---|---|---|
| massless gg→gg | 3.6e-15 | 1.7e-12 | 0 |
| massive (m=1.5) | 3.6e-15 | 1.9e-12 | 0 |

4-momentum conservation to **machine precision**, on-shell masses exact, invariant
mass = s exactly; the massless cross section matches Simpson (relerr 3.7e-5).
`VALIDATION: PASS`.

```bash
nvcc -O3 -std=c++17 -arch=sm_120 -o phase_space phase_space.cu && ./phase_space
```
