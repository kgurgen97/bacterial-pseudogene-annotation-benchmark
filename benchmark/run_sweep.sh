#!/usr/bin/env bash
# Gate the 40-seed sweep on complete fidelity, then run it.
#
# The sweep is only interpretable if the pool covers every pair the original
# used and the replay reproduces the original's control sets exactly. Both are
# checked first and the sweep does not run otherwise.
set -u
FLOOR_GB=${FLOOR_GB:-8}
cd "$(dirname "$0")/.."
avail() { df --output=avail -BG /mnt/c 2>/dev/null | tail -1 | tr -dc '0-9'; }
A=$(avail)
[ -n "$A" ] || { echo "cannot read host disk free space; refusing"; exit 3; }
echo "host C: free ${A}G, floor ${FLOOR_GB}G"
[ "$A" -ge "$FLOOR_GB" ] || { echo "below floor; refusing"; exit 3; }

source ~/miniforge3/etc/profile.d/conda.sh
conda activate hpwh
export PATH="$HOME/.local/bin-analysis:$PATH"

echo
echo "=== fidelity gate: all 99 pairs must be covered and reproduce exactly"
if ! python3 benchmark/check_pool_fidelity.py 2>&1 | tail -16; then
  echo
  echo "FIDELITY GATE FAILED. Not running the sweep."
  exit 4
fi

echo
echo "=== 40-seed sweep"
python3 analysis_seed_sensitivity.py --sweep 2>&1 | tee benchmark/logs/seed_sweep.log
