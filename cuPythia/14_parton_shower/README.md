# 14 — Batched parton shower (Sudakov veto) — the hard one

The parton shower is the **sequential** heart of event generation: emissions happen
in order, each depending on the last, so it does not vectorise within one event. The
GPU pattern is to batch **across** events — one shower per thread, thousands in
flight. Each shower evolves a scale `t = pT²` downward from `t_max` via the Sudakov
form factor; for a 1/t kernel with constant integrated splitting C, the next scale is
the exact inversion `t → t · R^(1/C)`.

The clean validation is the **Sudakov no-emission probability** and the Poisson
emission multiplicity:

## Result — RTX 5050, 8×10⁶ showers (C=0.078, ln(t_max/t_min)=9.21)
- no-emission fraction: MC **0.48873** vs Sudakov `exp(−C·L)` **0.48863** (rel 2e-4)
- mean multiplicity: MC **0.71593** vs analytic `C·L` **0.71615** (rel 3e-4)
- `VALIDATION: PASS`

The sequential-but-batched pattern, validated against first-principles physics.
Next on this axis: variable splitting kernels with the full veto (accept/reject) and
recoil/colour — the step toward a real device shower.

```bash
nvcc -O3 -std=c++17 -arch=sm_120 -o shower shower.cu && ./shower [nShowers]
```
