# 01 — σ(e⁺e⁻ → μ⁺μ⁻) GPU Monte-Carlo integration

First **physics** kernel: tree-level QED cross section integrated by GPU Monte
Carlo and validated against the closed form.

- Differential: `dσ/dΩ = α²/(4s) · (1 + cos²θ)`
- Total (analytic): `σ = 4π α² / (3s)`
- MC: sample cosθ ∈ [−1,1] (uniform in solid angle with φ; measure 4π), average
  the integrand, ×4π. Same host/device SplitMix64 RNG as `00_toolchain_check`.

## Result — RTX 5050, 5.24×10⁹ samples

| √s | analytic σ | GPU MC σ | rel. err | GPU rate | CPU rate (1t) | speedup |
|---|---|---|---|---|---|---|
| 10 GeV | 0.868545 nb | 0.868544 nb | 8.5e-7 | 1.72e10/s | 9.71e8/s | 17.7× |
| 20 GeV | 0.217136 nb | 0.217136 nb | 8.5e-7 | 1.70e10/s | 9.72e8/s | 17.5× |

`VALIDATION: PASS` (relerr < 1e-3). The **1/s scaling is exact** — doubling √s
quarters σ — confirming the physics, not just the integral.

## Build / run
```bash
nvcc -O3 -arch=sm_120 -o xsec xsec_ee_mumu.cu
./xsec [sqrt_s_GeV=10] [samplesPerThread=20000]
```

## Note
This integrand is φ-independent and fixed-angle; the next step toward real Pythia
is a **2→2 phase-space generator** (sampling full final-state kinematics) feeding
a matrix element — the data-parallel core of `PhaseSpace.cc` / the `Sigma*`
family. That target comes from the subsystem study.
