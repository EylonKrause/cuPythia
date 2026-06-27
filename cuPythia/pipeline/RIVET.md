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
| ⟨1−T⟩ (thrust mean) | ~0.063 | ~0.07 | small | — |
| thrust **distribution** (d54) | — | — | **~38% bin-avg, χ²/ndf 39** | **no hadron decays** (→ T→1 excess) |
| **charged multiplicity** (d01, √s=91.2) | **~11.3** | **20.73 ± 0.21** | **−45%** | **no hadron decays** (≈×2) |

**The honest bottom line:** the Rivet toolchain + HepMC3 interface *work* end-to-end — this is the
deliverable. The physics gap to **real data** is the sum of cuPythia's documented simplifications,
and the comparison cleanly identifies the **#1 missing piece: hadron decays.** With ρ/K*/ω/… stable,
cuPythia produces ≈half the charged multiplicity (11.3 vs 20.7) and too many pencil-like 2-jet events
(the T→1 thrust excess). ME corrections (`-DME_FIRST`) measurably help the hard 3-jet tail (χ²/ndf
47→39), but cannot close the decay gap. So the priority for *data* agreement is **decays**, ahead of
further shower/hadronization refinement — a conclusion only a real-data (Rivet) comparison reveals,
not a Pythia-vs-Pythia one. (Other caveats: no baryons; uds + g→qqbar c/b; pole/BW masses.)

## Scope note
The toolchain lives in `~/.local` + `~/rivetbuild/` (not vendored into the repo — it is large and
machine-specific). `build.sh` records the exact recipe. The HepMC3 interface, the hadron dump, and
`compare_rivet.py` are in-repo and reproducible given the toolchain.
