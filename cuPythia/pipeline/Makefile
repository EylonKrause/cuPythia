# cuPythia device-resident generator pipeline (separate from the demo kernels 00-14).
# Each stage reads/writes the SoA event record (event.cuh) on the device.
NVCC      ?= /usr/local/cuda/bin/nvcc
# --- GPU architecture (portable: Pascal sm_60 through Blackwell sm_120) -------------------------
# Default = one arch, fast. Target a different GPU:                 make ARCH=sm_60   (Pascal P100)
#                                                                   make ARCH=sm_61   (GTX 10-series)
# Portable multi-arch fatbinary (one binary runs on many GPUs; embeds PTX of the highest for
# forward JIT on future GPUs):     make SMS="60 61 70 75 80 86 89 90 120"   (any order is fine)
# IMPORTANT: CUDA 13.x REMOVED Pascal/Volta -- to emit sm_60/61/70 you need CUDA <= 12.9 (the last
# 12.x). The code itself is portable to sm_60 (plain FP64/int, no arch intrinsics); see PORTABILITY.md.
# NB the heavy hadronize_mr_baryon/max/hf targets are SINGLE-arch only (use ARCH=, not SMS=) -- a
# multi-arch fatbinary of that monster kernel would exhaust the compiler (see the __noinline__ note).
ARCH      ?= sm_120
SMS       ?=
ifeq ($(strip $(SMS)),)
  GENCODE := -arch=$(ARCH)
else
  # numeric-max arch via `sort -n` so SMS order does NOT matter (its PTX is embedded for forward JIT)
  SMS_HI  := $(shell echo $(SMS) | tr ' ' '\n' | sort -n | tail -1)
  GENCODE := $(foreach s,$(SMS),-gencode arch=compute_$(s),code=sm_$(s)) -gencode arch=compute_$(SMS_HI),code=compute_$(SMS_HI)
endif
NVCCFLAGS ?= -O3 -std=c++17 $(GENCODE)

STAGES = build_events reweight generate pdf_xsec pdf_lumi_fft shower_fsr zlund_test hadronize lhe_writer bridge region_test hadronize_mr decay_test baryon_test dalitz_test

all: $(STAGES)

build_events: build_events.cu event.cuh physics.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

reweight: reweight.cu physics.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

generate: generate.cu event.cuh physics.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

pdf_xsec: pdf_xsec.cu pdf.cuh physics.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

pdf_lumi_fft: pdf_lumi_fft.cu
	$(NVCC) $(NVCCFLAGS) -o $@ $< -lcufft

shower_fsr: shower_fsr.cu shower_inc.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

# Optional shower precision variants (compile flags; NOT in `make check`, so the committed LL
# results stay byte-stable -- each flag is proven byte-identical when OFF). Each adds physics:
#   ME_FIRST    = direct-ME (POWHEG-style) first/hardest emission
#   USE_CMW     = CMW soft-coherence alpha_s rescaling
#   GLUON_SPLIT = final-state g->qqbar splitting (secondary qqbar pairs; forks the colour string).
#     Validated vs Pythia weightGluonToQuark=1: N_gqq +2.2%; reference counter in thrust_pythia (g2q1).
shower_fsr_me: shower_fsr.cu shower_inc.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -DME_FIRST -o $@ $<
shower_fsr_cmw: shower_fsr.cu shower_inc.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -DUSE_CMW -o $@ $<
shower_fsr_gsplit: shower_fsr.cu shower_inc.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -DGLUON_SPLIT -o $@ $<
# g->qqbar with Pythia's DEFAULT option-4 weight (zCosThe reshape + high-mass damping) -> matches
# the default rate 0.5233 (+1.9%) instead of the plain-kernel 0.5656; massless kinematics (rate match).
shower_fsr_g2q4: shower_fsr.cu shower_inc.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -DGLUON_SPLIT -DG2QQ_WEIGHT4 -o $@ $<
# 2-loop alpha_s (Pythia AlphaStrong order=2, bit-identical). NB it ~doubles the LL parton
# multiplicity (over-radiates) -- a leading-log shower wants 1-loop; offered for completeness.
shower_fsr_2loop: shower_fsr.cu shower_inc.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -DAS_2LOOP -o $@ $<
# ZFLAV = initial qqbar flavour from the Z hadronic branching fractions (default build: always d).
# Makes b/c primaries -> realistic heavy-hadron content (D/B mesons); pair with -DHFDECAY for the
# Z->bb/cc track boost. Validated: charge conserved per event, D~0.39/B~0.37 per event (right rate).
shower_fsr_zflav: shower_fsr.cu shower_inc.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -DZFLAV -o $@ $<

bridge: bridge.cu shower_inc.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

zlund_test: zlund_test.cu ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

hadronize: hadronize.cu zlund_inc.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

region_test: region_test.cu region_inc.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

hadronize_mr: hadronize_mr.cu region_inc.cuh zlund_inc.cuh shower_inc.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

decay_test: decay_test.cu decay_inc.cuh bw_inc.cuh shower_inc.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

# Dalitz ME-shape validator: decays omega/eta at rest, checks the P-wave edge-suppression and the
# eta linear slope <y> appear ONLY with -DDALITZ_ME (the flat build is the control). In `make check`
# we build+run the FLAT control (default, byte-stable); the ME build is exercised by `make checkhf`.
dalitz_test: dalitz_test.cu decay_inc.cuh bw_inc.cuh shower_inc.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -DDECAYS -o $@ $<

# Standalone validator of the Lund baryon flavour mechanism (copies the flavour functions from
# hadronize_mr.cu; tiny kernel -> fast compile, unlike the full baryon hadronizer build).
baryon_test: baryon_test.cu ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

# Optional hadronization precision variants (compile flags; NOT in `make check`):
#   USE_BW      = Breit-Wigner vector-meson masses
#   GLUON_SPLIT = consume the forked multi-string events from g->qqbar (the shower also splits);
#     findStrings slices the chain at each q->qqbar boundary -> independent strings, conservation
#     exact (1.4e-10). A non-forked event is one string -> identical to the default hadronizer.
hadronize_bw: hadronize.cu zlund_inc.cuh bw_inc.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -DUSE_BW -o $@ $<
hadronize_mr_bw: hadronize_mr.cu region_inc.cuh zlund_inc.cuh shower_inc.cuh bw_inc.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -DUSE_BW -o $@ $<
hadronize_mr_gsplit: hadronize_mr.cu region_inc.cuh zlund_inc.cuh shower_inc.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -DGLUON_SPLIT -o $@ $<
# DECAYS = decay the primary unstable hadrons (rho/K*/omega/phi/eta/eta'/K0 -> pi/K) into the ALEPH
# particle-level stable set. Best-physics build (ME+g->qqbar+decays): charged mult 11.3 -> 18.7 (vs
# ALEPH 20.73), thrust chi2/ndf vs ALEPH 39 -> 10.3. decay_inc.cuh; validated by decay_test (in check).
hadronize_mr_decay: hadronize_mr.cu region_inc.cuh zlund_inc.cuh shower_inc.cuh decay_inc.cuh bw_inc.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -DME_FIRST -DGLUON_SPLIT -DDECAYS -o $@ $<
# BARYONS = Lund diquark baryon production (p/Lambda/Sigma/...). The flavour mechanism is validated
# standalone by baryon_test (charge + baryon number conserved over 2M strings); this full build adds
# it to the shower+hadronizer+decays chain. The combined kernel only compiles because the heavy device
# functions are __noinline__ (one monstrous inlined __global__ otherwise exhausts cicc/ptxas) -> now
# -O2 is fine (~10 min). Decuplet/Sigma0 resonances decay (decay_inc.cuh, -DBARYONS).
hadronize_mr_baryon: hadronize_mr.cu region_inc.cuh zlund_inc.cuh shower_inc.cuh decay_inc.cuh bw_inc.cuh ../common/rng.cuh
	$(NVCC) -O2 -std=c++17 -arch=$(ARCH) -DME_FIRST -DGLUON_SPLIT -DDECAYS -DBARYONS -o $@ $<
# Full-chain "best realistic physics": ME + g->qqbar + decays + baryons + Z-flavour init + heavy-flavour
# (D/B) decays + Dalitz ME shapes. This is the build used for the real-data (ALEPH) comparison. Heavy
# (~10-12 min at -O2; __noinline__ makes it tractable). HFDECAY needs ZFLAV/gsplit to make the c/b.
hadronize_mr_max: hadronize_mr.cu region_inc.cuh zlund_inc.cuh shower_inc.cuh decay_inc.cuh bw_inc.cuh ../common/rng.cuh
	$(NVCC) -O2 -std=c++17 -arch=$(ARCH) -DME_FIRST -DGLUON_SPLIT -DDECAYS -DBARYONS -DZFLAV -DHFDECAY -DDALITZ_ME -o $@ $<
# Same realistic flavour + heavy decays but WITHOUT baryons (baryons are untuned -> they pull mult the
# wrong way; this is the cleaner multiplicity/thrust build for ALEPH). ~8 min.
hadronize_mr_hf: hadronize_mr.cu region_inc.cuh zlund_inc.cuh shower_inc.cuh decay_inc.cuh bw_inc.cuh ../common/rng.cuh
	$(NVCC) -O2 -std=c++17 -arch=$(ARCH) -DME_FIRST -DGLUON_SPLIT -DDECAYS -DZFLAV -DHFDECAY -DDALITZ_ME -o $@ $<

lhe_writer: lhe_writer.cu ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

# Build + validate every pipeline stage.
check: all
	./build_events 500000
	./reweight 1000000
	./generate 1000000
	./pdf_xsec 4000
	./pdf_lumi_fft 15
	./shower_fsr 200000
	./zlund_test 2000000
	./hadronize 200000
	./lhe_writer 10000 events.lhe
	./bridge 5000 shower_partons.dat
	./region_test
	./hadronize_mr 20000
	./decay_test 500000
	./baryon_test 2000000
	./dalitz_test 1000000

# Exercise the optional heavy-flavour / Dalitz precision flags (separate from `check` so the committed
# byte-stable results are untouched). Validates D/B decays (conservation + GPU==CPU) and the Dalitz ME
# shape (omega edge-suppression + eta <y><0 appear only with -DDALITZ_ME).
checkhf: decay_test.cu dalitz_test.cu decay_inc.cuh bw_inc.cuh shower_inc.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -DDECAYS -DHFDECAY -DDALITZ_ME -o decay_test_hf decay_test.cu
	./decay_test_hf 500000
	$(NVCC) $(NVCCFLAGS) -DDECAYS -o dalitz_test_flat dalitz_test.cu
	$(NVCC) $(NVCCFLAGS) -DDECAYS -DDALITZ_ME -o dalitz_test_me dalitz_test.cu
	@echo "--- Dalitz: flat control vs ME (ME must edge-suppress omega and give eta <y> < 0) ---"
	./dalitz_test_flat 2000000
	./dalitz_test_me 2000000

# Optional real-PDF tool — fills the device PDF grid from Pythia's real gluon PDF and
# validates the interpolator vs it (needs Pythia + PYTHIA8DATA; not in `make check`).
PY8 ?= ../../pythia8317
genpdf: genpdf.cc
	g++ -O2 $< -o $@ $$($(PY8)/bin/pythia8-config --cxxflags --libs)

# Optional: validate the 2-loop alpha_s (the -DAS_2LOOP path) bit-identically against Pythia's
# own AlphaStrong order=2, and print the Lambda_{3,4,5} constants (needs Pythia + PYTHIA8DATA).
as2loop_validate: as2loop_validate.cc
	g++ -O2 -std=c++17 $< -o $@ $$($(PY8)/bin/pythia8-config --cxxflags --libs)

# Optional HepMC3 output — needs HepMC3 installed in $(HEPMC3) (built from source to
# ~/.local; not in `make check` since it is a non-standard dependency). Validates by readback.
HEPMC3 ?= $(HOME)/.local
hepmc3_writer: hepmc3_writer.cc
	g++ -O2 $< -o $@ -I$(HEPMC3)/include -L$(HEPMC3)/lib -lHepMC3 -Wl,-rpath,$(HEPMC3)/lib

clean:
	rm -f $(STAGES) hepmc3_writer decay_test_hf dalitz_test_flat dalitz_test_me \
	      hadronize_mr_full hadronize_mr_max hadronize_mr_hf hadronize_mr_baryon hadronize_mr_zf hadronize_mr_zfhf

.PHONY: all check checkhf clean
