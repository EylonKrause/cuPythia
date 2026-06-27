# Pythia 8.317 — upstream patch package

Eleven **behavior-preserving** optimizations found by the cuPythia subsystem audit
(`../AUDIT.md`, 32-agent study, each finding adversarially re-verified at `file:line`;
7 of 25 raw findings were rejected by that pass). They are read-only / fewer-copies /
O(1)-lookup changes — no physics changes — packaged here as individually-reviewable
upstream merge requests. (The two *correctness* bugs B1/B2 are already applied in this
repo, commit `da7db7f`.)

Pythia 8 is developed on its own GitLab (gitlab.com/Pythia8/releases). These are ready
to submit there as separate small MRs; they are **not** auto-pushed from here.

## Verified

| id | file | change | benefit | status |
|----|------|--------|---------|--------|
| **OPT1** | `src/SigmaProcess.cc` + `SigmaProcess.h` | precompute the id→beam-slot index in `initFlux`; replace the per-channel `inBeamA/inBeamB` linear scan in `sigmaPDF` with O(1) indexed reads | removes ~O(sizePair·sizeBeam) int-compares **per phase-space trial** (the hottest per-event path) | **VERIFIED: `examples/main101` (HardQCD, 8 TeV) bit-identical (sha256) before/after a full rebuild** → `OPT1-sigmaPDF-slot-lookup.patch` |

OPT1 is the flagship: it is on the hot Monte-Carlo path *and* main101 exercises it, so
bit-identity is a real regression. The patch in this directory is the verified diff.

## Documented (on the default-QCD path — bit-identity verifiable the same way)

| id | file | change |
|----|------|--------|
| OPT2 | `src/SimpleTimeShower.cc:2598`, `SimpleSpaceShower.cc:855` (+headers) | external `pTnext` takes `Event` **by value** inside `noEmissionProbability`'s 10000-trial loop → `Event&`/`const Event&` + `const vector<…DipoleEnd>&` (a deep `Event` copy per trial). SpaceShower can take `const Event&`; TimeShower needs non-const `Event&` (it forwards to non-const `pT2nextQCD`/`finiteCorrection`) — minimal blast radius. |
| OPT4 | `src/StringFragmentation.cc:744` | `rapPairs = colConfig.rapPairs` deep-copies a `vector<vector<pair<double,double>>>` every `fragment()`, read-only → `const&`. |
| OPT5 | `src/StringFragmentation.cc:1143,1203,1224,1255,1307…` | `StringRegion X = system.region(...)` copies an 8-`Vec4` struct; `region()` returns a reference → `const StringRegion&` per site. **Caveat (found by attempting it): the locals call `StringRegion` methods (`pHad`, `particle`, …) that are not marked `const`, so the MR must first const-ify those accessors (they are read-only) — otherwise `const&` fails to compile. A real, slightly larger MR, not a one-liner.** |

## Documented (off the default path — needs a targeted regression, not main101)

| id | file | change |
|----|------|--------|
| OPT3 | `src/StringFragmentation.cc:3100` | `kappaEffModifier(StringEnd, vector<int>)` by value (StringEnd embeds ~100-double StringFlav), deep-copied per hadron when ClosePacking on → `const&` both. |
| OPT6 | `src/PhaseSpace.cc:4274` | `Rambo::genPoint` `vector<double> mIn` by value + per-call heap churn → `const&` + `reserve`/reuse. |
| OPT7 | `src/SimpleTimeShower.cc:2301` | `globalRecoilMode==2` rescans the whole event O(N) per dipole-end to count `nFinal`; dipole-independent → hoist once per `pTnext`. |
| OPT8 | `src/HadronWidths.cc:297` | `getResonances` returns `set<int>` by value per call, from the O(N²) rescattering screen → return `const&`/out-param. |
| OPT9 | `src/HadronLevel.cc:660,802` | screened `sigmaTotal` computed in `queueDecResc` then recomputed in `rescatter` → cache on the queued node. |
| OPT10 | `src/SigmaLowEnergy.cc:446,634` | `sigmaPartial`/`pickProcess` allocate fresh vectors per call then linear-scan → reuse member buffers / index directly. |
| OPT11 | `src/PhaseSpace.cc:618,772,930` | the set{1,2,3}Kin+sigmaPDF+weights block is triplicated verbatim → extract one inlined helper (maintainability). |

Full per-finding crafted diffs + PR bodies + targeted regressions were produced by the
audit (see `../AUDIT.md` for the verified `file:line` and exact fix for each). Each MR
should be submitted separately with its own targeted before/after comparison; the
off-path ones (OPT3/6–10) need a config that exercises their path (ClosePacking,
junctions, low-energy rescattering), not a blanket main101 run.

## How to submit
1. `git apply upstream/OPT1-sigmaPDF-slot-lookup.patch` onto a clean Pythia 8.317 checkout.
2. Rebuild, run the regression for that finding (for OPT1: `examples/main101`, diff stdout).
3. Open one MR per finding on gitlab.com/Pythia8 with the PR title/body and the regression result.
