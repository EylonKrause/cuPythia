# Conference abstract — CHEP / ACAT (also reusable for NVIDIA GTC)

Title:
**cuPythia: a reproducible, GPU-native reimplementation of PYTHIA 8 event generation**

Author: Eylon Krause (independent)

Abstract (~200 words):

Monte Carlo event generation is a major and growing fraction of the HL-LHC computing budget, yet the
parton shower and hadronization remain largely CPU-bound. We present **cuPythia**, an open-source
(GPL-2) from-scratch GPU (CUDA) reimplementation of parts of the PYTHIA 8 chain: the LO hard process, a
final-state dipole parton shower, Lund string fragmentation (single- and multi-region), and hadron
decays, executed entirely device-resident, one event per thread. A counter-based (SplitMix64) RNG makes
every event O(1)-reproducible and bit-identical between host and device, and makes N-weight
scale-variation reweighting bit-identical to N independent pinned re-runs — a reproducibility guarantee
the sequential CPU generators do not offer. Opt-in precision corrections (first-emission matrix element,
CMW and 2-loop alpha_s, g->qqbar splitting, Breit-Wigner masses, baryons, Z-flavour initialisation, D/B
decays, Dalitz matrix-element shapes) are validated against PYTHIA and against real ALEPH LEP1
event-shape data via Rivet (charged multiplicity 18.99, thrust chi2/ndf 5.22), with an honest account of
where corrections help and where they do not without tuning. cuPythia builds across NVIDIA Pascal
through Blackwell, auto-detects a compatible CUDA toolkit, and pools GPUs across heterogeneous devices
and LAN hosts for a single run. We discuss the design, the validation methodology, the LL/LO scope, and
the path toward NLO matching and tuning. cuPythia is a research/proof-of-concept port, not the official
PYTHIA.

Keywords: event generation, parton shower, Lund hadronization, GPU, CUDA, reproducibility, HL-LHC.
