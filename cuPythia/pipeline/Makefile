# cuPythia device-resident generator pipeline (separate from the demo kernels 00-14).
# Each stage reads/writes the SoA event record (event.cuh) on the device.
NVCC      ?= /usr/local/cuda/bin/nvcc
ARCH      ?= sm_120
NVCCFLAGS ?= -O3 -std=c++17 -arch=$(ARCH)

STAGES = build_events reweight

all: $(STAGES)

build_events: build_events.cu event.cuh physics.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

reweight: reweight.cu physics.cuh ../common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

# Build + validate every pipeline stage.
check: all
	./build_events 500000
	./reweight 1000000

clean:
	rm -f $(STAGES)

.PHONY: all check clean
