# Installing cuPythia (including on a cluster)

cuPythia is plain CUDA + a Makefile, so it installs three ways. All of them auto-/explicitly target the
right GPU architecture (Pascal `sm_60` → Blackwell `sm_120`); **targeting Pascal/Volta needs CUDA ≤ 12.9**.

## 1. Zero-config (single machine / login node with a GPU)

```bash
git clone https://github.com/EylonKrause/cuPythia && cd cuPythia
./run.sh                                   # detect GPU + CUDA, build, validate
./run.sh hadronize_mr_hf 1000000 events.dat
```
`run.sh` auto-detects every GPU and a compatible CUDA toolkit and shards across them; `cluster.sh`
extends that across LAN hosts (see [CLUSTER.md](CLUSTER.md)).

## 2. Environment module (`make install` into a prefix)

Build for the cluster's GPU, then install the binaries + orchestration scripts into a prefix you can
turn into a module:

```bash
cd cuPythia/pipeline
make ARCH=sm_80 NVCC=/usr/local/cuda-12.9/bin/nvcc shower_fsr hadronize_mr hadronize_mr_hf
make install PREFIX=/opt/cupythia/0.1.0          # -> /opt/cupythia/0.1.0/bin/cupythia-*
```
Binaries are namespaced `cupythia-<stage>` (e.g. `cupythia-hadronize_mr_hf`); `run.sh`→`cupythia`,
`cluster.sh`→`cupythia-cluster`. Add `$PREFIX/bin` to a modulefile and `module load cupythia`.
(For a multi-arch fatbinary that runs on heterogeneous nodes, use `SMS="80 90 120"` instead of `ARCH=`.)

## 3. Spack (the HPC standard)

```bash
spack repo add /path/to/cuPythia/packaging/spack
spack install cupythia cuda_arch=80              # A100
spack install cupythia cuda_arch=80,90           # A100 + H100 fatbinary
spack install cupythia cuda_arch=70 ^cuda@12.9   # Volta (pin CUDA <= 12.9)
spack install cupythia targets=shower_fsr,hadronize_mr,hadronize_mr_hf
spack load cupythia                              # cupythia-* on PATH
```
`packaging/spack/packages/cupythia/package.py` is a `MakefilePackage`+`CudaPackage` recipe.

## 4. Containers (Apptainer / Docker — best for reproducible cluster deployment)

**Apptainer / Singularity** (ubiquitous on HPC; needs only the NVIDIA driver on the host):
```bash
apptainer build cupythia.sif packaging/cupythia.def
apptainer run --nv cupythia.sif cupythia-hadronize_mr 100000 events.dat
```
**Docker / Podman:**
```bash
docker build -t cupythia -f packaging/Containerfile .
docker run --gpus all -v "$PWD":/work cupythia cupythia-hadronize_mr 100000 /work/events.dat
```
Both build a **multi-arch fatbinary** (`SMS="75 80 86 89 90 120"` by default) so one image runs on any
GPU in that range. To include Volta (`70`)/Pascal, set a `<=12.9` base (`--build-arg CUDA_TAG=...` /
edit the `.def`). The heavy `hadronize_mr_hf` is not prebuilt multi-arch; build it inside for your GPU:
`make ARCH=sm_<cc> hadronize_mr_hf`.

## Requirements

- NVIDIA GPU, Pascal (`sm_60`) or newer; the NVIDIA driver.
- CUDA toolkit **≥ 11.0** (C++17). For Pascal/Volta (`sm_60/61/70`) use CUDA **≤ 12.9**.
- A C++17 host compiler (g++/clang on Linux, MSVC on Windows via `cuPythia/build.ps1`).
- Optional tools (Pythia consumer checks, Rivet, HepMC3) are separate; see `cuPythia/pipeline/RIVET.md`.

See [PORTABILITY.md](PORTABILITY.md) for the full GPU/CUDA matrix and the FP64 performance caveat.
