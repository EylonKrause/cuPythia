# Research directions — where collision simulation is heading, and how cuPythia maps

Grounded in current literature (HSF Physics Event Generator WG, HEP-CCE, and the
GPU-generator groups). Honest about which directions cuPythia already touches and
which are open.

## The driver
HL-LHC needs roughly an order of magnitude more simulated events; Monte-Carlo
event generation is a recognized and growing fraction (~5–15% today) of the WLCG
CPU budget, and the community is reengineering generators for GPUs/accelerators
because that is where new HPC FLOPs live (Perlmutter, Frontier, Aurora, LUMI).
*(HSF Generator WG: arXiv:2004.13687, arXiv:2109.14938; "The Critical Importance
of Software for HEP" arXiv:2504.01050.)*

## Directions and cuPythia's status

| direction | key refs | cuPythia |
|---|---|---|
| GPU parton-level ME generation | Pepper arXiv:2311.06198; madgraph4gpu arXiv:2106.12631, 2510.05392 | gg→gg ME + MC ported & validated (01–03, 06) |
| **Unweighting efficiency** (headline metric) | Pepper; madgraph4gpu I/O | **kernel 07: η measured on GPU (≈10% for gg→gg) + LHE output** |
| Adaptive / neural importance sampling | MadNIS arXiv:2212.06172; flow matching arXiv:2506.18987 | uniform now; VEGAS-on-GPU = next kernel |
| Negative weights | Sherpa arXiv:2110.15211 | positive for gg→gg; relevant at NLO (open) |
| Generic event reweighting | MadtRex arXiv:2510.05100 | counter-based RNG (05) makes per-event recompute trivial; N-variation reweight = natural next kernel |
| Portability (CUDA+HIP+SYCL one-source) | HEP-CCE; arXiv:2203.09945 | raw CUDA, written to port; layer not yet added |
| Standard I/O (LHE / HepMC3 / Rivet) | — | **LHE added (07)**; HepMC3/Rivet need external libs (open) |
| Reproducibility for GRID production | Random123 / counter-based RNG | **kernel 05** |

## Who might use cuPythia, and what they'd need

- **madgraph4gpu / Pepper teams** — they have GPU matrix elements; the *downstream*
  half (shower, hadronization, unweighting, event I/O) is what cuPythia targets.
- **ATLAS/CMS production** — unweighting efficiency, reweighting, LHE/HepMC3, and
  GRID-reproducible events; cuPythia now has η + LHE + per-event reproducibility.
- **ML-generator groups (MadNIS, flow-matching)** — a fast, reproducible classical
  MC target to benchmark learned samplers against.
- **HEP-CCE** — a Pythia-flavoured testbed for portability-layer evaluation.
- **Heavy-ion (ALICE, sPHENIX) and EIC/ePIC** — the O(N²) hadronic rescattering
  kernel (roadmap) is most relevant there.

## Next kernels this motivates (in priority order)
1. **VEGAS / adaptive importance sampling** on GPU — directly raises the kernel-07
   unweighting efficiency (the biggest cost lever).
2. **On-the-fly reweighting** — N systematic-variation weights per event in one pass.
3. **Portability layer** (Alpaka or Kokkos) — one source → NVIDIA + AMD + Intel.
4. **O(N²) hadronic rescattering** — heavy-ion / EIC.

## References
arXiv:2004.13687 · 2109.14938 · 2504.01050 · 2311.06198 · 2106.12631 · 2510.05392 ·
2212.06172 · 2506.18987 · 2110.15211 · 2510.05100 · 2203.09945.
