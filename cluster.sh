#!/usr/bin/env bash
# cuPythia LAN cluster coordinator — pool GPUs across MULTIPLE machines for ONE generation run.
#
# Each host builds for its OWN GPU arch locally (a 5050 + a Jetson, a 4090 + an A100, any mix), computes
# a DISJOINT slice of the global counter-RNG event stream (event-index offset), and the per-host dumps
# are merged. Because every event's seed is a pure function of its GLOBAL index, the merged output is the
# same N-event sample no matter how it was split -- on same-arch hosts it is byte-identical to one GPU
# running all N; across mixed architectures it is statistically equivalent (see CLUSTER.md).
#
#   ./cluster.sh hosts.txt hadronize_mr_hf 1000000 events.dat
#
# hosts.txt -- one worker per line:  <ssh-target>   <cuPythia-dir-on-that-host>   [gpu-count]
#   localhost        /home/eylonk/cuPythia          # "localhost" runs locally (no SSH)
#   user@4090box     /home/user/cuPythia
#   nvidia@jetson    /home/nvidia/cuPythia      1   # optional: force a GPU count (skip the remote query)
#
# Requirements: passwordless SSH to each remote host; cuPythia checked out with a CUDA toolkit on each;
# the SAME <target> name. List the coordinator itself as bare "localhost" to avoid an SSH-to-self.
set -uo pipefail
say(){ printf '\033[1;35m[cluster]\033[0m %s\n' "$*" >&2; }
die(){ printf '\033[1;31m[cluster] %s\033[0m\n' "$*" >&2; exit 1; }

HOSTFILE="${1:-}"; TARGET="${2:-}"; N="${3:-}"; OUT="${4:-}"
{ [ -n "$HOSTFILE" ] && [ -n "$TARGET" ] && [ -n "$N" ] && [ -n "$OUT" ]; } \
  || die "usage: cluster.sh <hostfile> <target> <N> <out>   e.g. cluster.sh hosts.txt hadronize_mr_hf 1000000 events.dat"
[ -f "$HOSTFILE" ] || die "no hostfile: $HOSTFILE"
[ "$N" -ge 1 ] 2>/dev/null || die "N must be a positive integer (got '$N')"
case "$TARGET" in hadronize_mr*) ;; *) die "cluster runs only the hadronize_mr* generators (they write a mergeable per-event hadron dump)";; esac

is_local(){ case "$1" in
  localhost|127.0.0.1|127.0.1.1|::1|""|"$(hostname)"|"$(hostname -s 2>/dev/null)") return 0;;
  *@localhost|*@127.0.0.1|*@::1) return 0;;
  *) return 1;; esac; }
on(){ local h="$1"; shift; if is_local "$h"; then bash -lc "$*"; else ssh -o BatchMode=yes "$h" "$*"; fi; }
pull(){ local h="$1" r="$2" l="$3"; if is_local "$h"; then cp -f "$r" "$l"; else scp -q "$h":"$r" "$l"; fi; }

# ---- parse the hostfile ----
HOSTS=(); DIRS=(); FGPU=()
while read -r h d g _; do
  [ -z "${h:-}" ] && continue; case "$h" in \#*) continue;; esac
  [ -n "${d:-}" ] || die "hostfile line for '$h' is missing the cuPythia directory"
  HOSTS+=("$h"); DIRS+=("$d"); FGPU+=("${g:-auto}")
done < "$HOSTFILE"
[ "${#HOSTS[@]}" -ge 1 ] || die "no hosts in $HOSTFILE"

# ---- scratch dir + cleanup (removes the local scratch AND any remote shard parts, on every exit) ----
WORK="$(mktemp -d)"; META=()
cleanup(){
  if [ "${#META[@]}" -gt 0 ]; then
    for m in "${META[@]}"; do IFS='|' read -r ch _cd crp _ <<< "$m"; on "$ch" "rm -f '$crp'" >/dev/null 2>&1 || true; done
  fi
  rm -rf "$WORK" 2>/dev/null || true
}
trap cleanup EXIT

# ---- Phase 1 (parallel): per-host GPU count + build for that host's own arch ----
say "discovering ${#HOSTS[@]} host(s) and building each for its own GPU (parallel) ..."
for i in "${!HOSTS[@]}"; do
  h="${HOSTS[$i]}"; d="${DIRS[$i]}"; fg="${FGPU[$i]}"
  (
    if [ "$fg" = auto ]; then c=$(on "$h" "nvidia-smi --query-gpu=index --format=csv,noheader 2>/dev/null | wc -l" 2>/dev/null | tr -d ' '); else c="$fg"; fi
    printf '%s' "${c:-0}" > "$WORK/ng_$i"
    if on "$h" "cd '$d' && ./run.sh --build-only '$TARGET'" >"$WORK/log_$i" 2>&1; then printf ok > "$WORK/bd_$i"; else printf fail > "$WORK/bd_$i"; fi
  ) &
done
wait
NG=()
for i in "${!HOSTS[@]}"; do
  h="${HOSTS[$i]}"
  if [ "$(cat "$WORK/bd_$i" 2>/dev/null)" != ok ]; then
    sed 's/^/    /' "$WORK/log_$i" >&2 2>/dev/null || true
    die "host '$h': build failed. Run by hand:  ssh $h \"cd ${DIRS[$i]} && ./run.sh --build-only $TARGET\""
  fi
  c=$(cat "$WORK/ng_$i" 2>/dev/null)
  { [ -n "$c" ] && [ "$c" -ge 1 ] 2>/dev/null; } || die "host '$h': no GPU found (or unreachable). Check passwordless SSH + nvidia-smi."
  NG+=("$c"); say "  $h : $c GPU(s), built"
done

# ---- Phase 2: global plan (one contiguous slice of the global event stream per GPU) ----
G=0; for c in "${NG[@]}"; do G=$((G+c)); done
M=$(( (N + G - 1) / G ))
say "total $G GPU(s) across ${#HOSTS[@]} host(s); $N events, M=$M per shard"

# ---- Phase 3: launch every global shard (host, local GPU g, global shard index s) ----
pids=(); s=0; rem="$N"
for i in "${!HOSTS[@]}"; do
  h="${HOSTS[$i]}"; d="${DIRS[$i]}"
  for g in $(seq 0 $(( NG[i] - 1 ))); do
    cnt=$(( rem < M ? rem : M )); [ "$cnt" -le 0 ] && continue; rem=$((rem - cnt))
    rp="$d/.cluster_shard_$s.dat"
    on "$h" "cd '$d' && CUDA_VISIBLE_DEVICES=$g CUPYTHIA_SHARD=$s CUPYTHIA_SHARD_N=$M cuPythia/pipeline/$TARGET $cnt '$rp'" >/dev/null 2>&1 &
    pids+=("$!"); META+=("$h|$d|$rp|$s"); s=$((s + 1))
  done
done
[ "${#META[@]}" -ge 1 ] || die "nothing to generate (N=$N)"
say "running $s shard(s) ..."
fail=0; for p in "${pids[@]}"; do wait "$p" || fail=1; done
[ "$fail" = 0 ] || die "a shard failed (check per-host SSH / build / GPU / disk)"

# ---- Phase 4: collect each shard's dump back, validating it ----
locals=()
for m in "${META[@]}"; do
  IFS='|' read -r h d rp sidx <<< "$m"
  lp="$WORK/part_$sidx.dat"
  pull "$h" "$rp" "$lp" || die "could not collect $h:$rp"
  [ -s "$lp" ] || die "collected part is empty/missing: $h:$rp (shard $sidx short-wrote?)"
  hc=$(head -1 "$lp" | awk '{print $1}'); { [ -n "$hc" ] && [ "$hc" -ge 0 ] 2>/dev/null; } || die "part for shard $sidx has a bad header"
  locals+=("$lp")
done

# ---- Phase 5: merge in global-shard order (sum the VALID-event header counts, one header, concat bodies) ----
# NB the header counts VALID events (after refragment drops), so total < N is expected; sanity-bound it.
total=0; for lp in "${locals[@]}"; do total=$((total + $(head -1 "$lp" | awk '{print $1}'))); done
{ [ "$total" -gt 0 ] && [ "$total" -le "$N" ]; } || die "merged event count $total is implausible for N=$N (corrupt shard?)"
mz=$(head -1 "${locals[0]}" | awk '{print $2}')
{ echo "$total $mz"; for lp in "${locals[@]}"; do tail -n +2 "$lp"; done; } > "$OUT"
say "done: $total events (of $N requested; rest are refragment drops) from $G GPU(s) across ${#HOSTS[@]} host(s) -> $OUT"
