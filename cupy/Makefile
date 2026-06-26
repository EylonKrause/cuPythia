# cuPy GPU kernels. Override NVCC / ARCH on the command line as needed:
#   make ARCH=sm_120            # Blackwell / RTX 5050 (default)
#   make NVCC=/path/to/nvcc
NVCC      ?= /usr/local/cuda/bin/nvcc
ARCH      ?= sm_120
NVCCFLAGS ?= -O3 -arch=$(ARCH)

BINS = 00_toolchain_check/mc_pi 01_xsec_ee_mumu/xsec 02_qcd_2to2_me/qcd_2to2

all: $(BINS)

00_toolchain_check/mc_pi: 00_toolchain_check/mc_pi.cu common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

01_xsec_ee_mumu/xsec: 01_xsec_ee_mumu/xsec_ee_mumu.cu common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

02_qcd_2to2_me/qcd_2to2: 02_qcd_2to2_me/qcd_2to2.cu common/rng.cuh
	$(NVCC) $(NVCCFLAGS) -o $@ $<

# Build everything and run each kernel's self-validation.
check: all
	./00_toolchain_check/mc_pi
	./01_xsec_ee_mumu/xsec 10.0
	./02_qcd_2to2_me/qcd_2to2 2000000

clean:
	rm -f $(BINS)

.PHONY: all check clean
