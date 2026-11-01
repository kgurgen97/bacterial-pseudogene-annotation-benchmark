#!/usr/bin/env bash
# Annotate one benchmark assembly with one caller, recording everything.
#
# Plan section 5. Two flags carry the whole design: Bakta's
# --keep-contig-headers and DFAST's --use_original_name make each caller's
# output contig names identical to the input, so a lesion's known coordinate
# can be matched inside that caller's own output with no cross-caller name
# reconciliation. The parent study had to resolve 26% of 886,176 interval
# lookups by file order; this is how that is designed out instead of inherited.
#
# A run whose output contig names do not match the input is a FAILED run, not a
# run to be reconciled. check_contigs.py enforces that.
#
# Resumable: an existing, non-empty expected output means skip.
#
#   run_callers.sh <caller> <assembly.fna> <outdir> [threads]
set -eu
CALLER=$1
ASM=$(readlink -f "$2")
OUT=$3
THREADS=${4:-6}
ID=$(basename "$ASM" .fna)
DEST="$OUT/$CALLER/$ID"
BAKTA_DB=${BAKTA_DB:-$HOME/research/datasets/bakta-db/db}
RUNS=${RUNS:-$(cd "$(dirname "$0")" && pwd)/results/tool_runs.tsv}

source ~/miniforge3/etc/profile.d/conda.sh
conda activate callers

expected() {
  case "$CALLER" in
    bakta) echo "$DEST/$ID.gff3" ;;
    dfast) echo "$DEST/genome.gff" ;;
    *) echo "$DEST/done" ;;
  esac
}

if [ -s "$(expected)" ]; then
  echo "  skip $CALLER/$ID (already done)"
  exit 0
fi
mkdir -p "$DEST"

case "$CALLER" in
  bakta)
    VER=$(bakta --version 2>&1 | head -1)
    CMD="bakta --db $BAKTA_DB --output $DEST --prefix $ID --threads $THREADS
         --keep-contig-headers --translation-table 11 --force --skip-plot $ASM"
    ;;
  dfast)
    VER=$(dfast --version 2>&1 | head -1)
    CMD="dfast --genome $ASM --out $DEST --cpu $THREADS --use_original_name t
         --gcode 11 --force"
    ;;
  *)
    echo "unknown caller: $CALLER" >&2; exit 2 ;;
esac

# The timing file must NOT live in $DEST: dfast --force clears its output
# directory, which silently deleted it and left peak RSS unrecorded.
META="$OUT/.meta/$CALLER"
mkdir -p "$META"
LOGF="$META/$ID.run.log"
TIMEF="$META/$ID.time.txt"

START=$(date +%s)
/usr/bin/time -v -o "$TIMEF" $CMD > "$LOGF" 2>&1 || {
  echo "  FAILED $CALLER/$ID -- see $LOGF" >&2
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$ID" "$CALLER" "$VER" "FAILED" \
    "$(( $(date +%s) - START ))" "" "$CMD" >> "$RUNS"
  exit 1
}
SECS=$(( $(date +%s) - START ))
RSS=$(awk '/Maximum resident set size/ {print $NF}' "$TIMEF" 2>/dev/null)

if [ ! -s "$RUNS" ]; then
  printf 'assembly\tcaller\tversion\tstatus\twall_s\tmax_rss_kb\tcommand\n' > "$RUNS"
fi
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$ID" "$CALLER" "$VER" "ok" \
  "$SECS" "${RSS:-}" "$(echo $CMD)" >> "$RUNS"
echo "  done $CALLER/$ID in ${SECS}s, peak RSS ${RSS:-?} kB"
