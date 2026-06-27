# cuPythia → Rivet: real LEP data comparison (ALEPH_2004_S5765862)

This is the first comparison of cuPythia against **real experimental data** (not just Pythia):
the GPU events are written as HepMC3 and run through **Rivet**'s `ALEPH_2004_S5765862` LEP1
event-shape analysis, which carries the actual ALEPH reference histograms.

## The toolchain (built from source, NO sudo / NO pip → `~/.local`)

A genuine dependency slog — documented honestly so it is reproducible (`~/rivetbuild/build.sh`):
- **fastjet 3.4.0** (autotools) + **FastJet Contrib 1.051** (Rivet 3.1.x requires the contrib
  headers + the `fastjetcontribfragile` shared lib).
- **zlib 1.3.1** from source (the box has `libz.so.1` runtime but no `zlib.h` dev header; YODA
  and Rivet's gzipped reference data need it).
- **YODA 1.9.11** with `--disable-pyext` (no `Python.h`/`python3-dev`, no pip — pyext avoided).
- **HepMC3 3.3.1 rebuilt with `-DHEPMC3_ENABLE_SEARCH=ON`** — the original install omitted the
  Search module, so `HepMC3/Relatives.h` (which Rivet's `RivetHepMC.hh` needs) was missing.
- **Rivet 3.1.11** with `--disable-pyext` → the C++ driver **`rivet-nopy`** (the Python `rivet`
  wrapper needs pyext). `make install` skips only the Python doc-index; components installed
  directly. `rivet-nopy` is the libtool `.libs/` binary (the bin wrapper isn't installed).

Run command (literal paths; the WSL shell strips `$VAR`, so no variables):
```
LD_LIBRARY_PATH=~/.local/lib RIVET_ANALYSIS_PATH=~/.local/lib/Rivet \
RIVET_DATA_PATH=~/.local/share/Rivet RIVET_REF_PATH=~/.local/share/Rivet \
~/.local/bin/rivet-nopy cupythia_best.hepmc3 ALEPH_2004_S5765862
```

## The interface
- `hepmc3_writer.cc` emits an **e+e- beam pair** (status 4, ids ±11 at √s/2) so Rivet identifies
  the LEP1 beams (without them: "CANNOT FIND ANY BEAMS"). Validated by HepMC3 readback (1.4e-7).
- `hadronize_mr.cu` gains an optional per-event **hadron dump** (`./hadronize_mr N out.dat`,
  gated on argv — default `make check` behaviour byte-identical). `hadronize_mr_best` uses the
  full physics (`-DME_FIRST -DGLUON_SPLIT`). `compare_rivet.py` parses MC vs REF YODA blocks.

## Results — HONEST (vs real ALEPH data, 19.4k hadron-level events, best physics ME+g→qqbar)

| observable | cuPythia | ALEPH | gap | dominant cause |
|---|---|---|---|---|
| observable | cuPythia (no decay) | cuPythia **+decays** | ALEPH | residual w/ decays |
|---|---|---|---|---|
| **charged multiplicity** (√s=91.2) | ~11.3 | **18.69** | **20.73 ± 0.21** | **−9.8%** (no baryons/heavy-flavour floor) |
| thrust **distribution** (d54) χ²/ndf | 39 | **10.3** | — | a **3.8× better** fit (T→1 spike softens) |

**The real-data comparison drove the next correction — and it worked.** The no-decay Rivet run cleanly
identified hadron decays as the #1 gap; adding GPU recursive decays (`-DDECAYS`, see `decay_inc.cuh` /
PRECISION.md) closes most of it: charged multiplicity **11.3 → 18.69** (toward ALEPH 20.73) and the
thrust fit improves **χ²/ndf 39 → 10.3** (the T→1 pencil-2-jet excess softens as multiplicity rises),
with 4-momentum **conservation still exact** (1.23e-10) and the decay module bit-identical GPU≡CPU.
The residual −9.8% in charged multiplicity is the **documented floor**: no baryons (~5–6% of LEP1
tracks), heavy flavour only from g→qqbar (not a realistic Z→bb̄/cc̄ rate, D/B left stable), flat-Dalitz
3-body (no ME shape), untuned Lund params, and no detector simulation. This is the honest landing —
decays were necessary and sufficient to close most of the gap, the rest is structurally out of scope.
The Rivet toolchain + HepMC3 interface remain the reusable infrastructure deliverable; the conclusion
(decays first) is one only a real-data comparison could reveal, not a Pythia-vs-Pythia one.

## Scope note
The toolchain lives in `~/.local` + `~/rivetbuild/` (not vendored into the repo — it is large and
machine-specific). `build.sh` records the exact recipe. The HepMC3 interface, the hadron dump, and
`compare_rivet.py` are in-repo and reproducible given the toolchain.
