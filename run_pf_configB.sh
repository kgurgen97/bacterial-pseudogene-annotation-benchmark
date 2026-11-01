#!/usr/bin/env bash
# Pseudofinder configuration B, as PF_CONFIG_PLAN.md specifies.
#
# Identical to config A in every respect except --length_pseudo 0.50 (default
# 0.65). Same 523 genomes, same Swiss-Prot DIAMOND database, same --diamond, same
# 4 threads per genome.
#
# Config A ran serially at ~3 min a genome, 27 hours for the corpus. This runs 5
# at a time on 20 cores, so about 6 hours. Resumable: a genome with a finished
# pf_pseudos.gff is skipped, so an interruption costs at most the five in flight.
set -u

P=~/research/projects/n2-pseudogene-caller-bias
D=~/research/datasets
DB="$D/n2-genomes/sprot.dmnd"
OUT="$D/n2-pseudofinder-B"
LOG="$P/logs/pf_configB.log"
JOBS=5
LENGTH_PSEUDO=0.50          # the one change; plan section 2

mkdir -p "$OUT" "$P/logs"
exec >> "$LOG" 2>&1
echo "=================================================================="
echo "$(date '+%F %H:%M')  config B start, --length_pseudo $LENGTH_PSEUDO"

test -s "$DB" || { echo "FATAL: database missing at $DB"; exit 1; }

# Both genome locations, deduplicated BY ACCESSION rather than by path. The
# expansion genomes appear under n2-genomes/extracted as well as under
# n2-genomes-expansion, so a plain path sort counts 160 of them twice and gives
# 683. -L because these trees contain symlinks; without it find silently returns
# a subset, which cost this project a day once already.
list=$(mktemp)
python3 - "$list" <<'PY'
import pathlib, sys
D = pathlib.Path.home() / "research/datasets"
seen = {}
for root in ("n2-genomes/extracted", "n2-genomes-expansion"):
    base = D / root
    if not base.exists():
        continue
    for g in base.rglob("genomic.gbff"):
        seen.setdefault(g.parent.name, str(g))   # first location wins
pathlib.Path(sys.argv[1]).write_text("\n".join(seen[a] for a in sorted(seen)) + "\n")
print(f"  unique accessions: {len(seen)}")
PY
n=$(grep -c . "$list")
echo "  genomes to process: $n"
if [ "$n" -ne 523 ]; then
  echo "FATAL: expected 523 genomes, found $n. Not starting."
  rm -f "$list"; exit 1
fi

source ~/miniforge3/etc/profile.d/conda.sh
conda activate pseudofinder

run_one() {
  g="$1"
  acc=$(basename "$(dirname "$g")")
  out="$OUT/$acc"
  [ -s "$out/pf_pseudos.gff" ] && return 0
  mkdir -p "$out"
  python "$P/tools/pseudofinder/pseudofinder.py" annotate \
    --genome "$g" --outprefix "$out/pf" --database "$DB" \
    --diamond --threads 4 --length_pseudo "$LENGTH_PSEUDO" \
    > "$out/log.txt" 2>&1 \
    || echo "FAILED $acc" >> "$OUT/failures.txt"
}
export -f run_one
export P D DB OUT LENGTH_PSEUDO

echo "$(date '+%F %H:%M')  running $JOBS at a time"
xargs -a "$list" -I{} -P "$JOBS" bash -c 'run_one "$@"' _ {}
rm -f "$list"

done_n=$(find -L "$OUT" -maxdepth 2 -name pf_pseudos.gff | wc -l)
echo "$(date '+%F %H:%M')  finished: $done_n of 523"
[ -s "$OUT/failures.txt" ] && { echo "  failures:"; sort -u "$OUT/failures.txt" | sed 's/^/    /'; }

# Counts, by the same rule config A used and this project verified:
# non-comment lines in pf_pseudos.gff, FASTA headers in pf_intact.faa.
echo "$(date '+%F %H:%M')  counting"
python3 - <<'PY'
import pathlib, csv
OUT = pathlib.Path.home() / "research/datasets/n2-pseudofinder-B"
P = pathlib.Path.home() / "research/projects/n2-pseudogene-caller-bias"
rows = []
for d in sorted(OUT.glob("GCF_*")):
    gff, faa = d / "pf_pseudos.gff", d / "pf_intact.faa"
    if not (gff.exists() and faa.exists()):
        continue
    ps = sum(1 for l in gff.read_text().splitlines()
             if l.strip() and not l.startswith("#"))
    it = sum(1 for l in faa.read_text().splitlines() if l.startswith(">"))
    rows.append({"accession": d.name, "pfB_pseudo": ps, "pfB_intact": it})
out = P / "results/pf_configB_counts.tsv"
with out.open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["accession", "pfB_pseudo", "pfB_intact"],
                       delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print(f"  wrote {out} with {len(rows)} genomes")
PY
echo "$(date '+%F %H:%M')  DONE"
