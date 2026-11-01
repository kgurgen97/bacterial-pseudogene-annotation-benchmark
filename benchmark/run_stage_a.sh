#!/usr/bin/env bash
# Stage A: annotate every probe assembly with every runnable caller, STRICTLY
# one at a time.
#
# Bakta peaked at 10.3 GB resident on a 4.6 Mb genome. With 23 GB visible to
# WSL, two concurrent Bakta runs would leave nothing for the page cache and the
# second would be competing with the first for the 84 GB database. So this
# loops rather than parallelising, which is also what the owner asked for:
# one memory-intensive annotation process at a time.
#
# Resumable. run_callers.sh skips a run whose expected output already exists,
# so re-running this after an interruption costs only the remaining work.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
ASM=${ASM:-$HOME/research/datasets/n2-benchmark/assemblies}
OUT=${OUT:-$HOME/research/datasets/n2-benchmark/calls}
THREADS=${THREADS:-8}

# The disk that matters is the WINDOWS one. `df /` inside WSL reports the
# ext4.vhdx's virtual capacity -- 771 GB free of a 1 TB maximum -- while the
# file itself grows into whatever is left on C:, which is 14 GB. An earlier
# version of this guard read `df /` and would therefore never have fired until
# the host disk was already full.
FLOOR_GB=${FLOOR_GB:-8}

host_avail_gb() {
  if df --output=avail -BG /mnt/c > /dev/null 2>&1; then
    df --output=avail -BG /mnt/c | tail -1 | tr -dc '0-9'
  else
    echo ""      # unknown: treated as a stop condition below
  fi
}

echo "Stage A caller runs, serial, started $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "  assemblies: $ASM"
echo "  output:     $OUT"
echo "  wsl /    : $(df -h / | tail -1 | awk '{print $4}') free (virtual)"
echo "  host C:  : $(host_avail_gb)G free (the real budget, floor ${FLOOR_GB}G)"

fail=0
for caller in dfast bakta; do
  for fna in "$ASM"/*.fna; do
    id=$(basename "$fna" .fna)
    echo
    echo "--- $caller / $id"
    if ! bash "$HERE/run_callers.sh" "$caller" "$fna" "$OUT" "$THREADS"; then
      echo "    recorded as FAILED in results/tool_runs.tsv"
      fail=$((fail + 1))
    fi
    # a caller that fills the host disk would take Windows down with it, not
    # merely spoil a run, so stop rather than press on
    avail=$(host_avail_gb)
    if [ -z "$avail" ]; then
      echo "STOPPING: cannot read the host disk's free space, so the budget"
      echo "  cannot be checked. Refusing to continue."
      exit 3
    fi
    if [ "$avail" -lt "$FLOOR_GB" ]; then
      echo "STOPPING: only ${avail}G free on the host disk C:, floor is"
      echo "  ${FLOOR_GB}G. Refusing to continue."
      exit 3
    fi
  done
done

echo
echo "serial Stage A runs finished $(date -u +%Y-%m-%dT%H:%M:%SZ), $fail failures"
echo "  host C: now $(host_avail_gb)G free"
