#!/usr/bin/env bash
# cuPythia one-command launcher — no manual GPU compilation; multi-GPU & mixed-arch aware.
#
# On first run (or whenever the GPU set changes) it AUTO-DETECTS every GPU in the machine, picks a
# CUDA toolkit that can compile for ALL their microarchitectures, builds (a fatbinary if the GPUs
# differ), then runs. A generation run is automatically SHARDED across all GPUs and the outputs are
# merged — including a box with different GPUs (e.g. an A100 + an RTX 4090) used together.
#
#   ./run.sh                            # detect + build everything + validate (make check)
#   ./run.sh shower_fsr 200000          # build if needed, then run a stage
#   ./run.sh hadronize_mr_hf 1000000 events.dat   # generate 1e6 events ACROSS ALL GPUs -> events.dat
#
# Options (before the target):
#   --rebuild            force a clean rebuild
#   --arch sm_61         override the detected arch (single-arch build)
#   --cuda <dir>         use a specific CUDA toolkit (e.g. for Pascal/Volta on a CUDA-13 box)
#   --gpus 0,2,3         use only these GPU indices (default: all). "0,0" runs two shards on GPU 0.
#   --single             use one GPU even if several are present
#
# Physicists: never touch nvcc, gencode, CUDA versions, or device IDs — just ./run.sh <thing>.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPE="$HERE/cuPythia/pipeline"
STAMP="$PIPE/.cupythia_build"
say(){ printf '\033[1;36m[cuPythia]\033[0m %s\n' "$*" >&2; }
die(){ printf '\033[1;31m[cuPythia] %s\033[0m\n' "$*" >&2; exit 1; }

FORCE=0; ARCH_OVR=""; CUDA_OVR=""; GPUS_OVR=""; SINGLE=0; BUILD_ONLY=0
while [ $# -gt 0 ]; do case "$1" in
  --rebuild)    FORCE=1; shift;;
  --build-only) BUILD_ONLY=1; shift;;   # build for this host's GPU(s) and exit (used by cluster.sh over SSH)
  --arch)     [ $# -ge 2 ] || die "--arch needs a value, e.g. --arch sm_61"; ARCH_OVR="${2#sm_}"; shift 2;;
  --cuda)     [ $# -ge 2 ] || die "--cuda needs a path, e.g. --cuda /usr/local/cuda-12.9"; CUDA_OVR="$2"; shift 2;;
  --gpus)     [ $# -ge 2 ] || die "--gpus needs a list, e.g. --gpus 0,1"; GPUS_OVR="$2"; shift 2;;
  --single)   SINGLE=1; shift;;
  -h|--help)  sed -n '2,25p' "$0" | sed 's/^# \?//'; exit 0;;
  --) shift; break;;
  -*) die "unknown option $1 (see --help)";;
  *)  break;;
esac; done
TARGET="${1:-check}"; [ $# -gt 0 ] && shift || true

# ---- compute capability of one GPU index (12.0 -> 120), with a name-based fallback ----
cap_of(){
  local idx="$1" cc name
  cc=$(nvidia-smi -i "$idx" --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d ' .')
  if [ -n "$cc" ] && [ "$cc" -eq "$cc" ] 2>/dev/null; then echo "$cc"; return 0; fi
  name=$(nvidia-smi -i "$idx" --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)
  case "$name" in
    *P100*) echo 60;; *P40*|*P4*|*GTX\ 10*|*Titan\ Xp*) echo 61;;
    *V100*|*Titan\ V*) echo 70;; *T4*|*RTX\ 20*|*GTX\ 16*|*Quadro\ RTX*) echo 75;;
    *A100*|*A800*|*A30*) echo 80;; *RTX\ 30*|*RTX\ A2*|*RTX\ A4*|*RTX\ A5*|*RTX\ A6*|*A40*|*A10*|*A16*|*A2*) echo 86;;
    *RTX\ 40*|*L40*|*L4*|*Ada*) echo 89;; *H100*|*H200*|*GH200*) echo 90;;
    *RTX\ 50*|*B100*|*B200*|*GB200*) echo 120;; *) return 1;; esac
}

# ---- build the list of GPU indices to use ----
GIDS=()
if [ -n "$ARCH_OVR" ]; then
  GIDS=(0)                                   # explicit arch: a single logical device
elif [ -n "$GPUS_OVR" ]; then
  IFS=',' read -r -a GIDS <<< "$GPUS_OVR"
else
  command -v nvidia-smi >/dev/null 2>&1 || die "no nvidia-smi / GPU driver. Override with: ./run.sh --arch sm_XX $TARGET"
  while IFS= read -r line; do [ -n "$line" ] && GIDS+=("$(echo "$line" | tr -d ' ')"); done \
    < <(nvidia-smi --query-gpu=index --format=csv,noheader 2>/dev/null)
fi
[ "${#GIDS[@]}" -ge 1 ] || die "no GPUs detected. Override with: ./run.sh --arch sm_XX $TARGET"
[ "$SINGLE" = 1 ] && GIDS=("${GIDS[0]}")

# ---- per-GPU arch + the DISTINCT arch set (for a single fatbinary covering them all) ----
CAPS=(); DISTINCT=""
for g in "${GIDS[@]}"; do
  if [ -n "$ARCH_OVR" ]; then c="$ARCH_OVR"; else c="$(cap_of "$g" || true)"; fi
  [ -n "$c" ] || die "could not detect arch of GPU $g. Override with --arch sm_XX."
  CAPS+=("$c")
  case " $DISTINCT " in *" $c "*) ;; *) DISTINCT="$DISTINCT $c";; esac
done
DISTINCT="$(echo $DISTINCT | tr ' ' '\n' | sort -n | tr '\n' ' ' | sed 's/ *$//')"
NDISTINCT=$(echo $DISTINCT | wc -w)
say "GPUs: ${#GIDS[@]} [idx ${GIDS[*]}] -> arch(s): $(echo $DISTINCT | sed 's/\([0-9]*\)/sm_\1/g')"

# ---- newest CUDA toolkit that can EMIT every arch in the set ----
pick_nvcc(){
  local c nv a ok
  for c in $(ls -d /usr/local/cuda-* /usr/local/cuda 2>/dev/null | sort -Vr); do
    nv="$c/bin/nvcc"; [ -x "$nv" ] || continue
    ok=1; for a in $DISTINCT; do "$nv" --list-gpu-arch 2>/dev/null | grep -qx "compute_$a" || { ok=0; break; }; done
    [ "$ok" = 1 ] && { echo "$nv"; return 0; }
  done
  if command -v nvcc >/dev/null 2>&1; then
    ok=1; for a in $DISTINCT; do nvcc --list-gpu-arch 2>/dev/null | grep -qx "compute_$a" || { ok=0; break; }; done
    [ "$ok" = 1 ] && { command -v nvcc; return 0; }
  fi
  return 1
}
if [ -n "$CUDA_OVR" ]; then
  NVCC="$CUDA_OVR/bin/nvcc"; { [ -f "$NVCC" ] && [ -x "$NVCC" ]; } || NVCC="$CUDA_OVR"
  { [ -f "$NVCC" ] && [ -x "$NVCC" ]; } || die "no nvcc under --cuda $CUDA_OVR"
else
  NVCC="$(pick_nvcc || true)"
fi
if [ -z "${NVCC:-}" ] || [ ! -x "$NVCC" ]; then
  m="no installed CUDA toolkit can compile for arch set { $DISTINCT }."
  case " $DISTINCT " in *" 60 "*|*" 61 "*|*" 70 "*) m="$m Pascal/Volta need CUDA <= 12.9 (CUDA 13 removed them) -- install one or pass --cuda. See PORTABILITY.md.";; esac
  die "$m"
fi
say "CUDA: $("$NVCC" --version | grep -oE 'release [0-9.]+' | head -1)  ($NVCC)"

# ---- build flags: single arch, or a fatbinary over the distinct set (mixed-GPU box) ----
if [ "$NDISTINCT" -le 1 ]; then BUILDFLAG=(ARCH="sm_$DISTINCT"); else BUILDFLAG=(SMS="$DISTINCT"); fi

# ---- build once (rebuild only on first run, --rebuild, or a changed GPU/toolkit set) ----
KEY="archset[$DISTINCT]|$NVCC"
if [ "$FORCE" = 1 ] || [ ! -f "$STAMP" ] || [ "$(cat "$STAMP" 2>/dev/null)" != "$KEY" ]; then
  say "first run / new GPU set -> compiling ($(echo $DISTINCT | sed 's/\([0-9]*\)/sm_\1/g')) ..."
  rm -f "$STAMP"
  make -C "$PIPE" clean >/dev/null 2>&1 || true
fi
say "make $TARGET"
make -C "$PIPE" "${BUILDFLAG[@]}" NVCC="$NVCC" "$TARGET" || die "build failed for '$TARGET'"
printf '%s' "$KEY" > "$STAMP"
[ "$BUILD_ONLY" = 1 ] && { say "built $TARGET ($(echo $DISTINCT | sed 's/\([0-9]*\)/sm_\1/g'))"; exit 0; }

# ---- make-only targets: nothing to run ----
case "$TARGET" in check|checkhf|all|clean) say "done ($TARGET)"; exit 0;; esac
[ -x "$PIPE/$TARGET" ] || { say "'$TARGET' built (no runnable binary)"; exit 0; }

BIN="$PIPE/$TARGET"
G="${#GIDS[@]}"
# Multi-GPU only helps an event-generation run with an output file (clean, mergeable). hadronize_mr*
# writes the per-event hadron dump as argv[2]; shard it across GPUs and merge. Everything else (stats /
# validators / single GPU) runs on one device.
if [ "$G" -gt 1 ] && [ "$SINGLE" != 1 ] && [ $# -ge 2 ] && case "$TARGET" in hadronize_mr*) true;; *) false;; esac; then
  N="$1"; OUT="$2"
  M=$(( (N + G - 1) / G ))            # ceil: per-shard cap (also the disjoint-stream stride)
  say "sharding $N events across $G GPUs (~$M each) -> $OUT"
  pids=(); used=(); rem="$N"; k=0
  for i in $(seq 0 $((G-1))); do
    cnt=$(( rem < M ? rem : M )); [ "$cnt" -le 0 ] && continue; rem=$((rem-cnt))
    CUDA_VISIBLE_DEVICES="${GIDS[$i]}" CUPYTHIA_SHARD="$i" CUPYTHIA_SHARD_N="$M" \
      "$BIN" "$cnt" "$OUT.part$i" > "$OUT.log$i" 2>&1 &
    pids+=("$!"); used+=("$i"); k=$((k+1))
  done
  [ "${#used[@]}" -ge 1 ] || die "nothing to generate (N=$N)"
  fail=0; for p in "${pids[@]}"; do wait "$p" || fail=1; done
  [ "$fail" = 0 ] || { for i in "${used[@]}"; do echo "--- GPU ${GIDS[$i]} (shard $i) ---"; cat "$OUT.log$i"; done >&2; die "one or more shards failed"; }
  # merge the per-event hadron dumps: sum the header counts, keep one header, concat the bodies
  total=0; for i in "${used[@]}"; do total=$((total + $(head -1 "$OUT.part$i" | awk '{print $1}'))); done
  mz=$(head -1 "$OUT.part${used[0]}" | awk '{print $2}')
  { echo "$total $mz"; for i in "${used[@]}"; do tail -n +2 "$OUT.part$i"; done; } > "$OUT"
  for i in "${used[@]}"; do rm -f "$OUT.part$i" "$OUT.log$i"; done
  say "done: $total events from $k GPU shard(s) -> $OUT"
  exit 0
fi

# ---- single-GPU run (pin to the first selected device) ----
say "run on GPU ${GIDS[0]}: $TARGET $*"
exec env CUDA_VISIBLE_DEVICES="${GIDS[0]}" "$BIN" "$@"
