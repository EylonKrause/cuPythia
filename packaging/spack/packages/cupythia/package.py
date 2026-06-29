# Spack package for cuPythia — a GPU (CUDA) reimplementation of parts of PYTHIA 8.
#
# Install on an HPC cluster (then expose as a module):
#   spack repo add /path/to/cuPythia/packaging/spack
#   spack install cupythia cuda_arch=80          # one arch  (A100)
#   spack install cupythia cuda_arch=80,90       # fatbinary (A100 + H100)
#   spack load cupythia                          # puts cupythia-* on PATH
#
# NOTE: targeting Pascal/Volta (cuda_arch 60/61/70) requires CUDA <= 12.9 (CUDA 13 removed them);
# pin it with e.g.  spack install cupythia cuda_arch=70 ^cuda@12.9.
from spack.package import *


class Cupythia(MakefilePackage, CudaPackage):
    """An independent, from-scratch GPU (CUDA) reimplementation of parts of the PYTHIA 8 event
    generator (LO hard process, FSR dipole shower, Lund string fragmentation, hadron decays),
    device-resident with a counter-based RNG. Research / proof-of-concept port derived from
    PYTHIA 8.317 (GPL-2); not the official PYTHIA."""

    homepage = "https://github.com/EylonKrause/cuPythia"
    git = "https://github.com/EylonKrause/cuPythia.git"
    url = "https://github.com/EylonKrause/cuPythia/archive/refs/tags/v0.1.0.tar.gz"

    maintainers("EylonKrause")
    license("GPL-2.0-only")

    version("main", branch="master")
    version("0.1.0", tag="v0.1.0")

    # C++17 (the CUB compaction kernel) needs CUDA >= 11.0.
    depends_on("cuda@11:", type=("build", "run"))
    depends_on("gmake", type="build")

    # Build at least one GPU architecture.
    conflicts(
        "cuda_arch=none",
        msg="set a GPU arch, e.g. cuda_arch=80 (supported: 60,61,70,75,80,86,89,90,120)",
    )

    # Which pipeline generators to build (comma-separated). The heavy full-physics hadronizers
    # (hadronize_mr_hf / _max) are opt-in because they take ~10 min to compile per arch.
    variant(
        "targets",
        default="shower_fsr,hadronize_mr",
        description="comma-separated pipeline targets to build (e.g. add hadronize_mr_hf)",
        values=lambda x: True,
    )

    build_directory = join_path("cuPythia", "pipeline")

    def edit(self, spec, prefix):
        # MakefilePackage runs `make` in build_directory; we drive it via build_targets/args below.
        pass

    @property
    def common_make_args(self):
        spec = self.spec
        nvcc = join_path(spec["cuda"].prefix.bin, "nvcc")
        archs = list(spec.variants["cuda_arch"].value)
        args = ["NVCC=%s" % nvcc]
        if len(archs) == 1:
            args.append("ARCH=sm_%s" % archs[0])
        else:
            args.append("SMS=%s" % " ".join(archs))
        return args

    @property
    def build_targets(self):
        return self.common_make_args + self.spec.variants["targets"].value.split(",")

    @property
    def install_targets(self):
        return ["install", "PREFIX=%s" % self.prefix] + self.common_make_args
