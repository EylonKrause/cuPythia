# HSF Generators WG — forum post + lightning-talk abstract

## Forum / mailing-list post (paste, add the repo + DOI links)

Subject: cuPythia — a from-scratch GPU (CUDA) reimplementation of PYTHIA 8 event generation (open source)

Hi all,

I'd like to share **cuPythia**, an independent, open-source (GPL-2) GPU port of parts of the PYTHIA 8
event-generation chain, in case it's useful to the generators / GPU-offloading discussion. It is a
research / proof-of-concept effort and is **not** the official PYTHIA or affiliated with the
collaboration — I've tried to be careful about scope and validation.

What's implemented, device-resident and one-event-per-thread:
- LO hard process, an FSR dipole parton shower, Lund string fragmentation (single- and multi-region,
  gluon-kinked), hadron decays, and supporting machinery.
- A counter-based (SplitMix64) RNG: every event is O(1)-reproducible and **host/device bit-identical**,
  and N-weight scale-variation reweighting is **bit-identical** to N pinned re-runs.
- Opt-in precision corrections (first-emission ME, CMW & 2-loop alpha_s, g->qqbar splitting, BW masses,
  baryons, Z-flavour init, D/B decays, Dalitz ME shapes), each validated against PYTHIA and, where
  possible, against **real ALEPH LEP1 data via Rivet** (charged multiplicity 18.99, thrust chi2/ndf 5.22).
- Portability Pascal (sm_60) -> Blackwell (sm_120); auto-detects GPU + a compatible CUDA toolkit;
  pools GPUs across multiple devices and across LAN hosts (any mix of architectures) for one run, with
  the merged output bit-identical to a single-GPU run (same-arch).

Honest scope: LL/LO shower, untuned Lund parameters, NLO matching explicitly out of scope; the repo
documents every approximation (PRECISION.md, RIVET.md). For production physics, use the official PYTHIA.

Repo: <github.com/EylonKrause/cuPythia>   DOI: <zenodo>   Would value any feedback, and happy to give a
short WG talk if useful.

Thanks,
Eylon Krause

## Lightning-talk abstract (5-10 min WG slot)

**cuPythia: a reproducible, GPU-native reimplementation of PYTHIA 8 event generation.**
cuPythia is an open-source (GPL-2) CUDA port of the PYTHIA 8 chain — LO hard process, FSR dipole shower,
Lund string fragmentation, and decays — run entirely on the GPU, one event per thread, with a
counter-based RNG that makes events host/device bit-identical and O(1)-reproducible and makes N-weight
reweighting bit-identical to N re-runs. Opt-in precision corrections are validated against PYTHIA and
against real ALEPH LEP1 event-shape data via Rivet. It builds across NVIDIA Pascal->Blackwell, auto-
selects a compatible CUDA toolkit, and pools GPUs across heterogeneous devices and LAN hosts for a
single run. I'll cover the design (device-resident SoA, counter-RNG sharding), the validation (and an
honest account of where corrections did and did not help vs data), and the LL/LO scope and open
directions (NLO matching, tuning). Not the official PYTHIA; a proof-of-concept exploring GPU-native
generation for HL-LHC-scale computing.
