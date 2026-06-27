# cuPythia device-resident generator pipeline (separate from the demo kernels 00-14).
# Each stage reads/writes the SoA event record (event.cuh) on the device.
NVCC      ?= /usr/local/cuda/bin/nvcc
ARCH      ?= sm_120
NVCCFLAGS ?= -O3 -std=c++17 -arch=$(ARCH)

STAGES = build_events reweight generate pdf_xsec pdf_lumi_fft shower_fsr zlund_test

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

shower_fsr: shower_fsr.cu ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

zlund_test: zlund_test.cu ../common/rng.cuh
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

clean:
	rm -f $(STAGES)

.PHONY: all check clean
