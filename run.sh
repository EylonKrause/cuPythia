#!/usr/bin/env bash
# cuPythia one-command launcher — no manual GPU compilation required.
#
# On first run (or whenever it sees a different GPU) it AUTO-DETECTS your GPU's microarchitecture,
# picks a CUDA toolkit that can compile for it, builds the kernels for that arch ONCE, then runs.
# Every later run just runs the cached binary.
#
#   ./run.sh                       # detect + build everything + validate (make check)
#   ./run.sh shower_fsr 200000     # build if needed, then run a stage with its args
#   ./run.sh hadronize_mr_hf 50000 # the best-physics hadronizer (heavy: one-time ~10 min compile)
#   ./run.sh checkhf               # build + run the heavy-flavour / Dalitz validators
#
# Options (before the target):
#   --rebuild              force a clean rebuild
#   --arch sm_61           override the detected GPU arch
#   --cuda /usr/local/cuda-12.9   use a specific CUDA toolkit (e.g. for Pascal/Volta on a CUDA-13 box)
#
# Physicists: you should never need nvcc, gencode flags, or CUDA versions — just ./run.sh <thing>.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPE="$HERE/cuPythia/pipeline"
STAMP="$PIPE/.cupythia_build"
say(){ printf '\033[1;36m[cuPythia]\033[0m %s\n' "$*" >&2; }
die(){ printf '\033[1;31m[cuPythia] %s\033[0m\n' "$*" >&2; exit 1; }

# ---- optional flags, then: TARGET [args...] ----
FORCE=0; ARCH_OVR=""; CUDA_OVR=""
while [ $# -gt 0 ]; do case "$1" in
  --rebuild)  FORCE=1; shift;;
  --arch)     [ $# -ge 2 ] || die "--arch needs a value, e.g. --arch sm_61"; ARCH_OVR="${2#sm_}"; shift 2;;
  --cuda)     [ $# -ge 2 ] || die "--cuda needs a path, e.g. --cuda /usr/local/cuda-12.9"; CUDA_OVR="$2"; shift 2;;
  -h|--help)  sed -n '2,21p' "$0" | sed 's/^# \?//'; exit 0;;
  --) shift; break;;
  -*) die "unknown option $1 (see --help)";;
  *)  break;;
esac; done
TARGET="${1:-check}"; [ $# -gt 0 ] && shift || true

# ---- 1. detect compute capability -> build arch (12.0 -> 120, 6.1 -> 61) ----
detect_cc(){
  command -v nvidia-smi >/dev/null 2>&1 || return 1
  local cc name
  cc=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d ' .')
  if [ -n "$cc" ] && [ "$cc" -eq "$cc" ] 2>/dev/null; then echo "$cc"; return 0; fi
  # fallback for drivers older than --query-gpu=compute_cap: map the GPU name
  name=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)
  case "$name" in
    *P100*) echo 60;; *P40*|*P4*|*GTX\ 10*|*Titan\ Xp*) echo 61;;
    *V100*|*Titan\ V*) echo 70;; *T4*|*RTX\ 20*|*GTX\ 16*|*Quadro\ RTX*) echo 75;;
    *A100*|*A800*|*A30*) echo 80;; *RTX\ 30*|*RTX\ A2*|*RTX\ A4*|*RTX\ A5*|*RTX\ A6*|*A40*|*A10*|*A16*|*A2*) echo 86;;
    *RTX\ 40*|*L40*|*L4*|*Ada*) echo 89;; *H100*|*H200*|*GH200*) echo 90;;
    *RTX\ 50*|*B100*|*B200*|*GB200*) echo 120;; *) return 1;; esac
}

# ---- 2. newest installed CUDA toolkit that can EMIT that arch ----
pick_nvcc(){
  local cc="$1" c nv
  for c in $(ls -d /usr/local/cuda-* /usr/local/cuda 2>/dev/null | sort -Vr); do
    nv="$c/bin/nvcc"; [ -x "$nv" ] || continue
    "$nv" --list-gpu-arch 2>/dev/null | grep -qx "compute_$cc" && { echo "$nv"; return 0; }
  done
  if command -v nvcc >/dev/null 2>&1; then
    nvcc --list-gpu-arch 2>/dev/null | grep -qx "compute_$cc" && { command -v nvcc; return 0; }
  fi
  return 1
}

ARCH="${ARCH_OVR:-$(detect_cc || true)}"
[ -n "$ARCH" ] || die "could not detect a GPU (no nvidia-smi / no driver?). Override with: ./run.sh --arch sm_XX $TARGET"
say "GPU -> sm_$ARCH"

if [ -n "$CUDA_OVR" ]; then
  NVCC="$CUDA_OVR/bin/nvcc"; { [ -f "$NVCC" ] && [ -x "$NVCC" ]; } || NVCC="$CUDA_OVR"
  { [ -f "$NVCC" ] && [ -x "$NVCC" ]; } || die "no nvcc under --cuda $CUDA_OVR (expected $CUDA_OVR/bin/nvcc, or pass the nvcc binary directly)"
else
  NVCC="$(pick_nvcc "$ARCH" || true)"
fi
if [ -z "${NVCC:-}" ] || [ ! -x "$NVCC" ]; then
  m="no installed CUDA toolkit can compile for sm_$ARCH."
  case "$ARCH" in 60|61|70) m="$m Pascal/Volta were removed in CUDA 13.x -- install a CUDA <= 12.9 toolkit, then re-run (or pass --cuda /usr/local/cuda-12.9). See PORTABILITY.md.";; esac
  die "$m"
fi
say "CUDA -> $("$NVCC" --version | grep -oE 'release [0-9.]+' | head -1)  ($NVCC)"

# ---- 3. build once (rebuild only on first run, --rebuild, or a new GPU/toolkit) ----
KEY="sm_$ARCH|$NVCC"
if [ "$FORCE" = 1 ] || [ ! -f "$STAMP" ] || [ "$(cat "$STAMP" 2>/dev/null)" != "$KEY" ]; then
  say "first run / new GPU -> compiling for sm_$ARCH (one-time) ..."
  rm -f "$STAMP"                          # clear any stale stamp until THIS build actually succeeds
  make -C "$PIPE" clean >/dev/null 2>&1 || true
fi
say "make $TARGET (ARCH=sm_$ARCH)"
make -C "$PIPE" ARCH="sm_$ARCH" NVCC="$NVCC" "$TARGET" || die "build failed for '$TARGET'"
printf '%s' "$KEY" > "$STAMP"             # record success only after the build worked

# ---- 4. run (make-only targets like check/all/clean already did their work) ----
case "$TARGET" in
  check|checkhf|all|clean) say "done ($TARGET)"; exit 0;;
esac
[ -x "$PIPE/$TARGET" ] || { say "'$TARGET' built (no runnable binary of that name)"; exit 0; }
say "run: $TARGET $*"
exec "$PIPE/$TARGET" "$@"
