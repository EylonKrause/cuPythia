# 02 — Batched QCD 2→2 matrix element (gg→gg)

Verbatim GPU port of Pythia 8.317 `Sigma2gg2gg::sigmaKin` (`src/SigmaQCD.cc:115-129`)
— the study's **#1 GPU target**: pure branchless double arithmetic on
(ŝ, t̂, û, α_s), one CUDA thread per trial over an SoA batch.

## Validation — RTX 5050, 10⁷ trials
- GPU vs CPU (same formula): **relerr 3.0e-16** — bit-perfect port.
- Pythia formula vs textbook `(9/4)(3 − t̂û/ŝ² − ŝû/t̂² − ŝt̂/û²)`: **relerr 7.7e-16**.
  - The coefficient is 9/4, not 9/2: Pythia folds in the identical-gluon ½
    (`SigmaQCD.cc:126`). The cross-check initially failed at relerr **exactly 1.0**,
    which flagged the convention; verified by hand that Pythia's rearranged bracket
    sum `B = 2·(3 − t̂û/ŝ² − …)`. Independent check doing its job.
- `VALIDATION: PASS`

## Performance — the honest Amdahl lesson
| metric | time (10⁷) | speedup |
|---|---|---|
| CPU loop | 64 ms | — |
| GPU kernel only | 14 ms | **4.5×** |
| GPU incl. transfer | 49 ms (H2D 27 + kern 14 + D2H 8) | **1.3×** |

Only **1.3× end-to-end** — and that is the point. This kernel has low arithmetic
intensity (4 loads + ~30 flops, **15 of them slow FP64 divisions** per element) and
pays a full PCIe round-trip per batch, so it is transfer- and FP64-div-bound —
exactly what the study predicted ("modest multiple, NOT proportional to raw kernel
speedup; never exponential").

Contrast kernels 00/01 (**17–21×**): they do *many* trials per thread in registers
with *no* per-trial transfer. The takeaway driving the next kernel: **fuse
phase-space generation + RNG + ME so data stays GPU-resident**, instead of
transferring pre-generated arrays across PCIe.

## Build / run
```bash
nvcc -O3 -arch=sm_120 -o qcd_2to2 qcd_2to2.cu
./qcd_2to2 [nTrials=10000000]
```
