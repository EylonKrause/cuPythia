# Outreach kit — getting cuPythia in front of the HEP-generator community

Ready-to-send drafts. **You** send these (they're you, representing yourself) — nothing here is posted
automatically. Lead everywhere with the *validated, honest* results and the "not the official Pythia"
framing; that is the credibility.

## Channels & order (highest leverage first)

1. **HSF (HEP Software Foundation) — Generators working group.** → `hsf-generators-wg.md`
   - Post to the HSF forum / generators mailing list; ask for a short slot at a WG meeting.
   - HSF coordinates exactly this audience and amplifies via its newsletter.
2. **MCnet / the PYTHIA authors (Lund — Sjöstrand, Bierlich, Prestel, …).** → `email-pythia-mcnet.md`
   - A from-scratch GPU port of their code, validated vs PYTHIA + ALEPH, is directly interesting to them.
   - Best case: feedback / a link / collaboration. Be explicit it is not official PYTHIA.
3. **The GPU event-generator groups you already cite.** (madgraph4gpu, Pepper, MadFlow, HEP-CCE.)
   - Short email + the `RESEARCH_DIRECTIONS.md` mapping; ask where a GPU shower/hadronizer fits.
4. **Conferences:** CHEP and ACAT (talk + proceedings); CERN LPCC generator meetings (Indico). → `conference-abstract.md`
5. **NVIDIA:** GTC talk + NVIDIA Developer Blog (they feature HEP-on-GPU). Use the same abstract.
6. **Broad:** Show HN, r/CUDA, r/HPC, Bluesky/Mastodon(fediscience)/X. → `announcement.md`

## Before sending — make sure these are live
- [ ] GPL-2 `LICENSE` + `NOTICE` pushed (done).
- [ ] `v0.1.0` tag + GitHub release (done).
- [ ] Zenodo connected to the repo (manual, one click at zenodo.org → GitHub) and the release archived → **DOI**. Put the DOI badge in the README and the abstract.
- [ ] Repo topics set (done) and the README headline results visible.
- [ ] `make check` green; `INSTALL.md` works.

## Standing honesty rules (your own)
- It is a research / proof-of-concept GPU port; **not** the official PYTHIA, **not** affiliated with the
  PYTHIA Collaboration or MCnet. Say so.
- Quote only validated numbers (Rivet/ALEPH, conservation, reproducibility) and state the scope/caveats
  (LL/LO shower, untuned Lund, NLO is out of scope). No over-claiming.
