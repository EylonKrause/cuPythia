# cuPythia GPU kernels. Override NVCC / ARCH on the command line as needed:
#   make ARCH=sm_120            # Blackwell / RTX 5050 (default)
#   make NVCC=/path/to/nvcc
# Cross-platform CMake build (Windows + Linux) lives in CMakeLists.txt.
NVCC      ?= /usr/local/cuda/bin/nvcc
ARCH      ?= sm_120
NVCCFLAGS ?= -O3 -arch=$(ARCH)

BINS = 00_toolchain_check/mc_pi 01_xsec_ee_mumu/xsec 02_qcd_2to2_me/qcd_2to2 \
       03_fused_resident_me/fused_gg2gg 04_multi_gpu/multi_gpu_gg2gg \
       05_reproducible_events/repro_events 06_fp32_mixed/fp32_gg2gg \
       07_unweighting_lhe/unweight_lhe

all: $(BINS)

00_toolchain_check/mc_pi: 00_toolchain_check/mc_pi.cu common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

01_xsec_ee_mumu/xsec: 01_xsec_ee_mumu/xsec_ee_mumu.cu common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

02_qcd_2to2_me/qcd_2to2: 02_qcd_2to2_me/qcd_2to2.cu common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

03_fused_resident_me/fused_gg2gg: 03_fused_resident_me/fused_gg2gg.cu common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

04_multi_gpu/multi_gpu_gg2gg: 04_multi_gpu/multi_gpu_gg2gg.cu common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -Xcompiler -pthread -o $@ $<

05_reproducible_events/repro_events: 05_reproducible_events/repro_events.cu common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

06_fp32_mixed/fp32_gg2gg: 06_fp32_mixed/fp32_gg2gg.cu common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

07_unweighting_lhe/unweight_lhe: 07_unweighting_lhe/unweight_lhe.cu common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

# Build everything and run each kernel's self-validation.
check: all
	./00_toolchain_check/mc_pi
	./01_xsec_ee_mumu/xsec 10.0
	./02_qcd_2to2_me/qcd_2to2 2000000
	./03_fused_resident_me/fused_gg2gg 1000
	./04_multi_gpu/multi_gpu_gg2gg 1000 --shards 4
	./05_reproducible_events/repro_events 1000000
	./06_fp32_mixed/fp32_gg2gg 2000
	cd 07_unweighting_lhe && ./unweight_lhe 1000

# Multi-node build (needs mpicxx in PATH); not part of `all` / `check`.
mpi:
	$(NVCC) $(NVCCFLAGS) -Xcompiler -pthread -DUSE_MPI -ccbin mpicxx \
	  -o 04_multi_gpu/multi_gpu_gg2gg 04_multi_gpu/multi_gpu_gg2gg.cu -lmpi

clean:
	rm -f $(BINS)

.PHONY: all check clean mpi
