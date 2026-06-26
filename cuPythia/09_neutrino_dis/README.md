# 09 — Neutrino-quark deep-inelastic scattering (parton model)

**Honest framing first:** neutrino-**nucleus** interactions (nuclear effects, Fermi
motion, final-state interactions) are the domain of **specialized generators —
GENIE, NuWro, GiBUU**, not Pythia. Pythia's role in neutrino physics is the DIS
final-state **hadronization**. cuPythia does **not** add nuclear physics, and air-
shower / ultra-high-energy forward physics (CORSIKA + SIBYLL/QGSJET/EPOS) is
likewise out of scope — I won't fake it.

What this **can** add honestly is the **partonic** neutrino DIS cross section — the
textbook electroweak result neutrino experiments rest on, and the signature that
revealed valence quarks:

- `dσ/dy(ν q)    ∝ 1`        (flat in inelasticity y)
- `dσ/dy(ν q̄) ∝ (1 − y)²`
- `dσ/dy = (G_F² s / π) · shape(y)`

## Result — RTX 5050, E_ν = 100 GeV
| channel | ⟨shape⟩ (analytic) | σ |
|---|---|---|
| ν q | 1.00000 (1) | 3.16×10⁻³⁶ cm² |
| ν q̄ | 0.33333 (1/3) | 1.05×10⁻³⁶ cm² |

**σ(νq)/σ(νq̄) = 3.0000** — the valence-quark signature — `VALIDATION: PASS`.
σ ≈ 10⁻³⁶ cm² is the correct order of magnitude for νN DIS. (Only the analytic
y-structure is here; the x-dependence needs real PDFs / LHAPDF.)

## Build / run
```bash
nvcc -O3 -arch=sm_120 -o nu_dis nu_dis.cu
./nu_dis [trialsPerThread=4000]
```
