# Short announcements (Show HN / Reddit / Bluesky / Mastodon / X)

Add the repo link + the Zenodo DOI. Keep the honest framing.

## One-liner
From-scratch CUDA reimplementation of PYTHIA 8 hadronization — validated against real ALEPH LEP data,
runs Pascal→Blackwell, and pools GPUs across a LAN (any mix of architectures). GPL-2, not the official PYTHIA.

## Show HN
**Show HN: cuPythia — a GPU (CUDA) reimplementation of PYTHIA 8 particle-physics event generation**

I rebuilt parts of PYTHIA 8 (the Monte Carlo generator used across the LHC experiments) to run entirely
on the GPU: the hard process, the parton shower, Lund string fragmentation, and decays, one event per
thread. A counter-based RNG makes every event bit-identical between CPU and GPU and O(1)-reproducible,
and makes N-weight reweighting bit-identical to N re-runs. The opt-in physics corrections are validated
against PYTHIA and against real ALEPH LEP1 data via Rivet. It builds for every NVIDIA GPU from Pascal to
Blackwell, auto-detects a compatible CUDA toolkit, and can pool GPUs across several machines on a LAN
(e.g. a 5050 laptop + a Jetson, or a 4090 + an A100) for one run. It's a research/proof-of-concept port
(GPL-2), not the official PYTHIA — and I've tried to document the scope and every approximation honestly.

## Reddit (r/CUDA, r/HPC, r/Physics)
Title: cuPythia — GPU (CUDA) port of PYTHIA 8 event generation: bit-identical reproducibility,
Pascal→Blackwell, multi-GPU + multi-host
Body: as Show HN above, trimmed; emphasize the relevant angle (CUDA design / HPC deployment / HEP physics).

## Bluesky / Mastodon (fediscience) / X
Built cuPythia: a from-scratch CUDA reimplementation of @PYTHIA 8 event generation 🚀 device-resident
shower + Lund hadronization, bit-identical CPU↔GPU reproducibility, validated vs real ALEPH LEP data,
Pascal→Blackwell, and it pools GPUs across a LAN (any arch mix). GPL-2, proof-of-concept (not official
PYTHIA). Repo + DOI 👇 #HEP #CUDA #GPU #HPC
