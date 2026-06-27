# cuPythia device-resident generator pipeline (separate from the demo kernels 00-14).
# Each stage reads/writes the SoA event record (event.cuh) on the device.
NVCC      ?= /usr/local/cuda/bin/nvcc
ARCH      ?= sm_120
NVCCFLAGS ?= -O3 -std=c++17 -arch=$(ARCH)

STAGES = build_events reweight generate pdf_xsec pdf_lumi_fft shower_fsr zlund_test hadronize lhe_writer bridge region_test hadronize_mr

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
	rm -f $(STAGES) hepmc3_writer

.PHONY: all check clean
