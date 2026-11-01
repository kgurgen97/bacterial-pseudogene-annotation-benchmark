#!/usr/bin/env bash
# Launch the deterministic omega pool build, checking the HOST disk first.
#
# The host disk is the constraint, not WSL's virtual one: `df /` inside WSL
# reports the vhdx's 1 TB maximum while the file grows into whatever is left on
# C:. The floor is checked BEFORE the job starts, not after it finishes.
set -u
FLOOR_GB=${FLOOR_GB:-8}
cd "$(dirname "$0")/.."

avail() { df --output=avail -BG /mnt/c 2>/dev/null | tail -1 | tr -dc '0-9'; }
A=$(avail)
if [ -z "$A" ]; then
  echo "cannot read the host disk's free space; refusing to start"; exit 3
fi
echo "host C: free ${A}G, floor ${FLOOR_GB}G"
if [ "$A" -lt "$FLOOR_GB" ]; then
  echo "below the floor; refusing to start"; exit 3
fi

source ~/miniforge3/etc/profile.d/conda.sh
conda activate hpwh
export PATH="$HOME/.local/bin-analysis:$PATH"
command -v diamond >/dev/null || { echo "diamond not on PATH"; exit 2; }

mkdir -p benchmark/logs
echo "started $(date -u +%Y-%m-%dT%H:%M:%SZ)"
python3 analysis_seed_sensitivity.py --pool 2>&1 | tee benchmark/logs/omega_pool.log | tail -30
echo "finished $(date -u +%Y-%m-%dT%H:%M:%SZ), host C: now $(avail)G free"
