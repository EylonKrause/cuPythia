# Pythia 8.317 audit — cuPythia subsystem study

Findings from a 6-subsystem, 32-agent study of Pythia 8.317, **each
adversarially re-verified at file:line** (17 confirmed, 7 rejected, 1 uncertain).
Two correctness bugs are fixed in this repo; the rest are documented here as
individually-reviewable upstream-PR candidates. They are micro-optimizations on
real code paths, but `examples/main101` does not exercise most of them
(hadronization ClosePacking, junctions, merging, low-energy rescattering), so
each needs a *targeted* regression, not a blanket one — hence documented, not
blind-applied.

## Fixed in this repo (correctness, verified, `main101` bit-identical)

| # | location | bug | fix | commit |
|---|---|---|---|---|
| B1 | `src/SigmaProcess.cc:1171` | `Sigma3Process::store3Kin` isSChannel branch overrides `Q2RenSave` (renorm) instead of `Q2FacSave` (factor) for fixed factorization scale — copy-paste; clobbers the renorm scale set just above. Latent (no shipped 2→3 s-channel process). | `Q2RenSave` → `Q2FacSave`, matching `store1Kin:806`. | da7db7f |
| B2 | `src/Basics.cc:109-110` | `Rndm::pick` do/while dereferences `prob[++index]` before the bound test, so an FP residue (`flat()` ~1 ULP < 1) can read `prob[size]` and return `index==size` — OOB for callers (`SigmaEW.cc:1417`, `SigmaLowEnergy.cc:641/674`). ~2⁻⁵³ trigger. | bound → `size()-1`; last element still selectable, never OOB. | da7db7f |

## Do NOT "fix"

- **`src/FragmentationFlavZpT.cc:1124-1131` (`StringZ::initShape`, bug, confirmed).**
  When `useOldAExtra` is set, the `(1-z)^a` exponent uses the OLD endpoint's
  strange/diquark a-shift. It is **wrong on purpose** — kept behind the flag for
  tune backwards-compatibility (self-documented in code). Touching it changes
  physics tunes. Leave it.
- **`src/SimpleSpaceShower.cc:1099-1213` (ISR PDF refresh, confirmed inefficiency).**
  Re-evaluates ~2·nQuarkIn+1 `xfISR` interpolations per dipole refresh — but this
  is inherent/intentional (in-code comment). It is a **GPU batching target**, not
  a correctness fix.

## Needs physics-reference confirmation before ANY edit (risky)

- **`src/SigmaLowEnergy.cc:232-240` (`PelaezpiK32ElData`, possible data bug, confirmed).**
  The `LinearInterpolator` for the πK isospin-3/2 *elastic* cross section begins
  with `0.64527, 1.800` — which are exactly the interpolator's own left/right mass
  bounds passed on the same line. They look like the constructor arguments
  accidentally pasted into the data array. **If** real, it distorts the low-energy
  πK elastic σ. But correcting it requires the original Pelaez parametrization
  values — do not guess physics data. Flag upstream with the reference in hand.

## PR-ready inefficiencies (confirmed read-only / behavior-preserving)

Highest impact first.

1. **`src/SigmaProcess.cc:454-466` (hot path).** `sigmaPDF` linear-scans
   `inBeamA/inBeamB` for every channel in `inPair` every trial — O(sizePair·sizeBeam)
   (~1000 int-compares for a qq flux), though the id→beam-slot map is fixed after
   `initFlux`. → precompute the slot index in `initFlux`, store per `InPair`,
   O(1) lookup. *On the per-event Monte-Carlo path — the one with real throughput
   impact.*
2. **`src/SimpleTimeShower.cc:2598` / `src/SimpleSpaceShower.cc:855`.** External
   `pTnext` overload takes `Event` **by value**; called inside
   `noEmissionProbability`'s 10000-trial loop → a full `Event` deep-copy per trial.
   → `const Event&` (+ `const vector<TimeDipoleEnd>&`). *Virtual-signature change;
   match the base declaration.*
3. **`src/StringFragmentation.cc:3100-3103`.** `kappaEffModifier(StringEnd end,
   vector<int> partonList)` by value; `StringEnd` embeds a ~100-double `StringFlav`
   block, deep-copied per produced hadron when ClosePacking is on. → `const&` both
   (decl `include/Pythia8/StringFragmentation.h:257-259`).
4. **`src/StringFragmentation.cc:744`.** `rapPairs = colConfig.rapPairs;` deep-copies
   a `vector<vector<pair<double,double>>>` every `fragment()` call, read-only. →
   `const&` (and the two callee params).
5. **`src/StringFragmentation.cc:1143,1203,1224,1255,1307,…`.** `StringRegion X =
   system.region(...)` copies an 8-`Vec4` struct by value; `region()` returns a
   reference and X is read-only. → `const StringRegion&` (verify each site).
6. **`src/PhaseSpace.cc:4274-4326` (`Rambo::genPoint`).** `vector<double> mIn` by
   value (read-only) + `energies/mXi/energiesXi` heap-churned per call. →
   `const vector<double>&` + `reserve()`/reuse buffers (decl `PhaseSpace.h:670`).
7. **`src/SimpleTimeShower.cc:2301-2304`.** `globalRecoilMode==2` rescans the whole
   event (O(N)) per dipole-end to count `nFinal`; dipole-independent. → hoist above
   the dipole loop (once per `pTnext`).
8. **`src/HadronWidths.cc:297-330` (`getResonances`).** Returns `set<int>` by value
   (heap red-black tree) per call, reached from the O(N²) rescattering screen
   (`SigmaLowEnergy.cc:876`). → return `const&` / fill an out-param.
9. **`src/HadronLevel.cc:802` & `:660`.** For a surviving hadron pair, `sigmaTotal`
   is computed in `queueDecResc`, then the full chain re-runs in `rescatter` via
   `pickProcess`. → cache the screened σ on the queued node.
10. **`src/SigmaLowEnergy.cc:446-449,634-642`.** `sigmaPartial`/`pickProcess`
    allocate fresh `vector<int>`/`vector<double>` per call to fetch one partial σ,
    then linear-scan. → reuse member buffers / index directly.
11. **`src/PhaseSpace.cc:618-653,772-807,930-967` (maintainability).** The
    set{1,2,3}Kin + `sigmaPDF` + weights block is triplicated verbatim across the
    three setup scans — divergence hazard, not runtime. → extract one inlined helper.

## Rejected by verification (NOT real — examples of the filter working)

`PhaseSpace.cc:704-711` "loop-invariant hoist" (the cited `log`/divisions actually
depend on the loop-varying `ratio34`/`unity34`), plus 6 others. The adversarial
pass rejected 7/25 raw findings — the reason these are worth trusting.
