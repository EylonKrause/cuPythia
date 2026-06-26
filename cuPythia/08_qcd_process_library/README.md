# 08 — QCD 2→2 process library on GPU

**Honest framing first:** Pythia does **not** lack these interactions — it has them
all (and far more: EW, Higgs, top, SUSY, a large BSM suite). What's genuinely at
the frontier is *precision* (NLO/NNLO, via external matching to POWHEG/MadGraph)
and *high-multiplicity* hard processes (via external ME generators) — not missing
fundamental interactions. So this kernel doesn't invent physics; it **broadens
cuPythia's GPU coverage** from one process to the complete tree-level QCD 2→2 set.

Each is a **verbatim port of Pythia 8.317** (`src/SigmaQCD.cc`) cross-checked **on
GPU** against the independent textbook (Ellis-Stirling-Webber / Combridge) analytic
form. Massless light quarks.

## Result — RTX 5050, 1.05×10⁹ trials/process, |cosθ|<0.8
| process | σ_cut [pb] | Pythia == textbook |
|---|---|---|
| gg→gg | 5.99e4 | PASS (0 mismatches) |
| qg→qg | 2.61e4 | PASS |
| qq'→qq' | 1.07e4 | PASS |
| qqbar→gg | 1.72e3 | PASS |
| gg→qqbar | 4.83e2 | PASS |

`VALIDATION: PASS` — all 5 agree to <1e-12. The cross-section ordering is physical
(gluon-rich channels dominate).

## Build / run
```bash
nvcc -O3 -arch=sm_120 -o qcd_library qcd_library.cu
./qcd_library [trialsPerThread=4000]
```
