#!/usr/bin/env bash
# Wait for Bakta to finish the 523-genome corpus, then rerun the two analyses that
# the plans require on all four callers.
#
# Why a watcher rather than running the analyses now: both scripts admit a caller
# only when it covers the genome set the test uses, so running early would simply
# reproduce the three-caller answer. This fires once, when the fourth caller is
# actually complete.
#
# It will NOT run the analyses on a partial corpus. If Bakta's runner exits while
# short of 523, that is reported as a failure and nothing is recomputed.
set -u
P=~/research/projects/n2-pseudogene-caller-bias
R=$P/results
LOG=$P/logs/four_caller_rerun.log
TARGET=523

mkdir -p "$P/logs"
exec >> "$LOG" 2>&1
echo "=================================================================="
echo "$(date '+%F %H:%M')  watcher armed, waiting for bakta to reach $TARGET"

rows() { echo $(( $(wc -l < "$R/bakta_counts.tsv") - 1 )); }

last=-1
while :; do
  n=$(rows)
  if [ "$n" -ge "$TARGET" ]; then
    echo "$(date '+%F %H:%M')  bakta reached $n of $TARGET"
    break
  fi
  # If no bakta process remains and we are still short, the run died.
  if ! pgrep -f "run_bakta_all.sh|bin/bakta " > /dev/null 2>&1; then
    sleep 120   # let a final write land before concluding
    n=$(rows)
    if [ "$n" -lt "$TARGET" ]; then
      echo "$(date '+%F %H:%M')  STOPPED SHORT: bakta is at $n of $TARGET and no"
      echo "  bakta process is running. NOT rerunning the analyses on a partial"
      echo "  corpus. Investigate \$HOME/research/datasets/bakta_expansion.log,"
      echo "  then rerun this watcher or the two analyses by hand."
      exit 1
    fi
    break
  fi
  if [ "$n" -ne "$last" ]; then
    echo "$(date '+%F %H:%M')  bakta $n / $TARGET"
    last=$n
  fi
  sleep 300
done

cd "$P"
source ~/miniforge3/etc/profile.d/conda.sh
conda activate hpwh

echo
echo "$(date '+%F %H:%M')  ===== analysis_c3_strata.py"
python3 analysis_c3_strata.py
echo
echo "$(date '+%F %H:%M')  ===== analysis_literature_retest.py"
python3 analysis_literature_retest.py
echo
echo "$(date '+%F %H:%M')  ===== numbers in the write-ups vs the json"
python3 check_quoted_numbers.py || echo "  (write-ups still quote three-caller numbers; expected until they are updated)"

echo
echo "$(date '+%F %H:%M')  ===== did all four callers get admitted?"
python3 - <<'PY'
import json, pathlib
R = pathlib.Path("results")
c3 = json.load(open(R / "c3_strata.json"))
lit = json.load(open(R / "literature_retest.json"))
print("  C3 callers        :", c3["callers"])
print("  literature callers:", lit["callers"])
for label, got in (("C3", c3["callers"]), ("literature", lit["callers"])):
    if "Bakta" not in got:
        print(f"  WARNING: {label} still excluded Bakta")
PY
echo
echo "$(date '+%F %H:%M')  DONE. Four-caller numbers are in results/; write-ups not yet updated."
