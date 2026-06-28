# cuPythia LAN cluster — pool GPUs across machines for one run

`cluster.sh` lets you generate one event sample using the GPUs of **several machines on your LAN at
once** — any mix of architectures. A 5050 laptop + a Jetson, a 4090 box + an A100 node, a pile of
heterogeneous workstations: each host builds cuPythia for **its own** GPU, computes a disjoint slice of
the global event stream, and the coordinator merges the slices into one output file.

```bash
./cluster.sh hosts.txt hadronize_mr_hf 1000000 events.dat
```

## The hostfile

One worker machine per line (`hosts.example.txt` is a template):

```
# <ssh-target>      <cuPythia-dir-on-that-host>     [gpu-count]
localhost           /home/eylonk/cuPythia
eylonk@4090-box     /home/eylonk/cuPythia
nvidia@jetson       /home/nvidia/cuPythia       1
eylonk@a100-node    /opt/cuPythia
```

- `localhost` runs locally (no SSH).
- Remote hosts need **passwordless SSH** (key auth; `ssh -o BatchMode=yes` is used), cuPythia checked
  out at the given path, and a CUDA toolkit. The same `<target>` name must be buildable on each.
- The optional 3rd column forces a GPU count; otherwise it is queried with `nvidia-smi`.

## How it works

1. **Discover + build** — for each host it queries the GPU count and runs `./run.sh --build-only
   <target>`, so every host compiles for *its own* arch (a Jetson builds `sm_87`, a 4090 `sm_89`, an
   A100 `sm_80`, a 5050 `sm_120` — automatically; mixed CUDA toolkits are fine).
2. **Plan** — with `G` total GPUs it sets `M = ceil(N/G)` and gives each GPU a contiguous block of the
   global event index space.
3. **Run** — each GPU runs its binary with `CUPYTHIA_SHARD=<global index>`, `CUPYTHIA_SHARD_N=M`,
   `CUDA_VISIBLE_DEVICES=<local gpu>`; the kernel adds `shard·M` to its local event index, so shard *s*
   computes exactly global events `[s·M, s·M+cnt)`.
4. **Collect + merge** — each host's per-event hadron dump is pulled back (`scp`, or a local copy for
   `localhost`) and merged in global-shard order (sum the header counts, one header, concat the bodies).

Because every event's seed is a pure function of its **global** index, the split is irrelevant to the
result — the same `N` events come out no matter how many machines or GPUs you used.

## Reproducibility — honest

- **Same-architecture hosts** (e.g. several 4090s, or 5050 + 5050): the merged output is
  **byte-for-byte identical** to one GPU running all `N` events — *verified* (a 2-shard cluster run and
  the single-GPU run produce the same sha256).
- **Mixed-architecture hosts** (e.g. 5050 + Jetson, 4090 + A100): the counter-RNG *seed* of each event
  is integer math, bit-identical on every arch — so the sharding itself is exact. But FP64
  transcendentals (`sin/cos/exp/log/pow`) differ by a few ULPs between architectures' math libraries, so
  momenta agree only to **~1e-12** (the GPU↔CPU transcendental limit cuPythia already documents). The
  **vast majority of events are identical**; however, because the shower and fragmentation use
  veto / accept-reject steps, a ULP-level difference can occasionally flip a branch, so a **small
  fraction of events differ in particle content** (not just momenta). The merged result is therefore a
  valid, statistically-equivalent physics sample — reproducible *given the same host assignment* — but
  **not bit-identical** across a heterogeneous split. (Same-arch hosts have no such ULP divergence →
  byte-identical, above.) This is a property of IEEE FP64 across vendors' math libraries, not the sharding.

## Scope / limitations

- Only the `hadronize_mr*` generators (they write the mergeable per-event hadron dump).
- No auto-deploy: each host must already have the repo + a CUDA toolkit (clone once; `run.sh` builds).
- Validated here on `localhost` acting as multiple hosts (proves the orchestration, sharding and merge,
  and the byte-identical same-arch property); real multi-machine SSH transport is the same code path.
- Transport is plain SSH/`scp` over the LAN — fine for the modest dump sizes; no shared filesystem
  required (use one if you have it by pointing every host's dir at it).
