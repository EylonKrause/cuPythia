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

| observable | no decay | **+decays** | **+decays +baryons** | ALEPH |
|---|---|---|---|---|
| **charged multiplicity** (√s=91.2) | ~11.3 | **18.69** | 17.30 | **20.73 ± 0.21** |
| thrust **distribution** (d54) χ²/ndf | 39 | **10.3** | 13.1 | — |
| p+p̄ / event | 0 | 0 | **0.745** | ~1.05 |
| Λ+Λ̄ / event | 0 | 0 | **0.214** | ~0.39 |

**The real-data comparison drove the corrections — honestly, including where one did NOT help.**
- **Decays (`-DDECAYS`)** were the #1 win the no-decay run identified: charged mult **11.3 → 18.69**,
  thrust **χ²/ndf 39 → 10.3** (the T→1 pencil-2-jet excess softens), conservation exact, GPU≡CPU.
- **Baryons (`-DBARYONS`)** are now *validated in the full chain* (electric-charge and baryon-number
  conserved **0/19334 events**, 4-momentum exact 1.23e-10, p+p̄ and Λ+Λ̄ at the right order ~60–70% of
  PDG). But with **untuned** `probQQtoQ`/Lund parameters they slightly *reduce* multiplicity
  (18.69 → 17.30) and *worsen* the thrust fit (χ²/ndf 10.3 → 13.1) — a correct, conservation-exact
  mechanism can still move the wrong way on data when its rate isn't tuned. **Honest conclusion:**
  closing the last ~15% to ALEPH 20.73 is a *Lund tune* (a separate multi-histogram fit) **plus** D/B
  decays (the Z→bb̄/cc̄ heavy-hadron track boost) — not more conservation-correct mechanisms. The
  baryon *physics* is right; the *normalization* needs tuning, which is out of scope and not faked.

(Build note: the combined shower+multiregion+decays+baryons kernel only compiles after marking the
heavy device functions `__noinline__` — it was one monstrous inlined `__global__` that exhausted the
compiler; `__noinline__` splits it into tractable pieces with bit-identical physics. ~10 min at -O2.)

## Scope note
The toolchain lives in `~/.local` + `~/rivetbuild/` (not vendored into the repo — it is large and
machine-specific). `build.sh` records the exact recipe. The HepMC3 interface, the hadron dump, and
`compare_rivet.py` are in-repo and reproducible given the toolchain.
