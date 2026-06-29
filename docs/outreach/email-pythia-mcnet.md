# Draft email — the PYTHIA authors / MCnet (and the GPU-generator groups)

Look up current addresses yourself (pythia.org / the authors' institutional pages / MCnet). Keep it
short and respectful; you are sharing, not claiming endorsement.

---

Subject: A from-scratch GPU (CUDA) reimplementation of PYTHIA 8 event generation — feedback welcome

Dear PYTHIA authors,

I'm an independent developer and I've written **cuPythia**, an open-source (GPL-2) GPU port of parts of
PYTHIA 8.317 — the LO hard process, an FSR dipole parton shower, Lund string fragmentation (single- and
multi-region), and hadron decays — running device-resident, one event per thread, with a counter-based
RNG. I want to be upfront that it is a research / proof-of-concept reimplementation, **not** affiliated
with or endorsed by the PYTHIA Collaboration, and it credits and ports from PYTHIA 8.317 under GPL-2.

I thought it might interest you because:
- It validates against PYTHIA at the LL/LO level and against **real ALEPH LEP1 data via Rivet** (charged
  multiplicity 18.99, thrust chi2/ndf 5.22), with exact per-event 4-momentum/charge/baryon-number
  conservation and bit-identical host/device reproducibility.
- The corrections are opt-in and individually validated (first-emission ME, CMW & 2-loop alpha_s,
  g->qqbar, BW masses, baryons, Z-flavour init, D/B decays, Dalitz ME shapes), and I documented honestly
  where a correction did *not* improve the data fit (e.g. untuned baryons), and that NLO is out of scope.
- It runs Pascal->Blackwell and pools GPUs across devices and LAN hosts for one run.

I'd be very grateful for any feedback — especially on correctness/scope framing, and on whether a
GPU-native shower/hadronizer is useful to the wider effort. Repo: github.com/EylonKrause/cuPythia
(DOI: <zenodo>). Thank you for PYTHIA, and for your time.

Best regards,
Eylon Krause
eylonk@advstg.com

---

## Variant for the GPU-generator groups (madgraph4gpu / Pepper / MadFlow / HEP-CCE)

Same opening, then: "Your work (madgraph4gpu / Pepper / MadFlow / ...) is cited in cuPythia's
RESEARCH_DIRECTIONS.md as the map I built against. cuPythia covers the shower+hadronization side
device-resident; I'd value your view on where a GPU-native Lund hadronizer + shower fits alongside the
GPU ME efforts, and on the reproducibility approach (counter-RNG, bit-identical reweighting)."
