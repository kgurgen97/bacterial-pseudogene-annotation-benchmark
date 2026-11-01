#!/usr/bin/env bash
# Recompute every interval-based analysis after the classifier fix.
#
# Dependency order matters: analysis_clustering.py reads dnds_pairs.tsv and
# analysis_dfast.py reads loci.tsv, so each must follow its producer.
#
# The pre-fix results are copied to results/pre_fix/ FIRST, untouched, so the
# before/after comparison is auditable and the old numbers are not lost. They
# are the numbers the write-ups published.
#
# Count-based analyses are deliberately NOT rerun: dfast_counts.tsv and
# bakta_counts.tsv come from each caller's own summary file, the runners got
# that right, and nothing about them changed. Rerunning them would only add
# noise to the diff.
set -u
cd ~/research/projects/n2-pseudogene-caller-bias
source ~/miniforge3/etc/profile.d/conda.sh
conda activate hpwh

# analysis_dnds.py shells out to a bare `diamond`, which lives in the callers
# env, not the analysis env. Prepending that whole env to PATH would shadow
# python and the rest of the toolchain, so instead a directory holding only
# symlinks to the binaries the analyses actually invoke goes first. Only
# diamond is invoked; blastp and makeblastdb are linked for completeness.
export PATH="$HOME/.local/bin-analysis:$PATH"
command -v diamond >/dev/null || { echo "diamond still not on PATH"; exit 2; }
echo "diamond: $(diamond --version 2>&1 | head -1) at $(command -v diamond)"

mkdir -p benchmark/logs results/pre_fix

echo "=== preserving the pre-fix results"
for f in dnds.json frame.json locus.json clustering.json dfast.json \
         dnds_pairs.tsv frame_pairs.tsv loci.tsv locus_diagnostics.tsv; do
  if [ -s "results/$f" ] && [ ! -e "results/pre_fix/$f" ]; then
    cp "results/$f" "results/pre_fix/$f"
    echo "  saved results/pre_fix/$f"
  else
    echo "  already saved or absent: $f"
  fi
done

run() {
  local script=$1 log=benchmark/logs/rerun_${1%.py}.log
  echo
  echo "=== $script  started $(date -u +%H:%M:%SZ)"
  if python3 "$script" > "$log" 2>&1; then
    echo "    ok, $(date -u +%H:%M:%SZ). tail:"
    tail -6 "$log" | sed 's/^/      /'
  else
    echo "    FAILED, see $log"
    tail -20 "$log" | sed 's/^/      /'
    return 1
  fi
}

run analysis_dnds.py      || exit 1
run analysis_clustering.py || exit 1
run analysis_locus.py     || exit 1
run analysis_dfast.py     || exit 1
run analysis_frame.py     || exit 1

echo
echo "=== all interval analyses recomputed $(date -u +%Y-%m-%dT%H:%M:%SZ)"
