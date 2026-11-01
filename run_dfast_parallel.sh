#!/usr/bin/env bash
# DFAST across 4 workers, 4 CPUs each.
#
# Why: DFAST averages only about 4 busy cores of 20 - much of its pipeline
# (MetaGeneAnnotator, Aragorn) is single-threaded - so one worker leaves the machine
# mostly idle while Bakta waits its turn. Four workers use ~16 cores and 2 GB total,
# cutting DFAST from roughly 12.6 h to 3.5 h and handing Bakta the full machine sooner.
# Bakta is NOT run alongside: its DIAMOND wants all 20 cores and its 46 GB database
# already exceeds RAM, which is what made the first parallel attempt 2.4x slower.
#
# Workers take genomes by index modulo 4, so none collide. Each appends to its own
# counts file; they are merged at the end. Resume is unchanged: a genome whose
# pseudogene_summary.tsv exists is skipped.
#
# Counting rule is the one fixed in run_dfast_all.sh: disrupted classes only
# (frameshift + internal_stop_codon + both). `partial` is recorded, never added.
set -u
source ~/miniforge3/etc/profile.d/conda.sh
conda activate callers

NW=${NW:-4}
CPU=${CPU:-4}
E=~/research/datasets/n2-genomes/extracted
OUT=~/research/datasets/n2-dfast
LOG=~/research/datasets/n2-dfast-logs
R=~/research/projects/n2-pseudogene-caller-bias/results
mkdir -p "$OUT" "$LOG"

count() { awk -F'\t' -v c="$2" 'NR>1 && $4==c {n++} END {print n+0}' "$1"; }

worker() {
  local id=$1
  local mine=0 new=0 fail=0 i=0
  local part="$R/dfast_counts.w$id.tsv"
  : > "$part"
  while IFS= read -r G; do
    if [ $(( i % NW )) -ne "$id" ]; then i=$((i+1)); continue; fi
    i=$((i+1))
    local ACC O S
    ACC=$(basename "$(dirname "$G")")
    O="$OUT/$ACC"; S="$O/pseudogene_summary.tsv"
    [ -s "$S" ] && continue
    mine=$((mine+1))
    # Go through dfast_runtime_fix.py, not plain dfast: DFAST 1.5.3 aborts with
    # AttributeError: 'HmmHit' object has no attribute 'id' on any genome where an
    # rRNA/tRNA overlaps an HMM-named CDS. Calling plain dfast here lost 3 of the 160
    # additions, two of them extreme-reduction endosymbionts - the stratum that matters
    # most to C3 - so the loss would have been biased, as it was for the original 21.
    if python "$(dirname "$0")/dfast_runtime_fix.py" -g "$G" -o "$O" --cpu "$CPU" --force \
         > "$LOG/$ACC.out" 2> "$LOG/$ACC.err"; then
      local fs is bo pa cd_
      fs=$(count "$S" frameshift)
      is=$(count "$S" internal_stop_codon)
      bo=$(count "$S" 'frameshift,internal_stop_codon')
      pa=$(count "$S" partial)
      cd_=$(awk -F'\t' '/^Number of CDSs/{print $2}' "$O/statistics.txt")
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$ACC" "$((fs+is+bo))" "$fs" "$is" "$bo" "$pa" "${cd_:-NA}" >> "$part"
      rm -f "$O/protein.faa" "$O/cds.fna" "$O/rna.fna" "$O/genome.fna" \
            "$O/genome.embl" "$O/genome.gbk"
      rm -rf "$O/ddbj"
      new=$((new+1))
    else
      fail=$((fail+1)); echo "FAILED $ACC" >> "$OUT/failures.txt"
    fi
    if [ $(( new % 20 )) -eq 0 ] && [ "$new" -gt 0 ]; then
      echo "$(date +%H:%M)  worker $id: new=$new fail=$fail"
    fi
  # -L is required: the expansion genomes are symlinked directories under $E and plain
  # find will not descend into a symlink, silently skipping them.
  done < <(find -L "$E" -name '*.fna' | sort)
  echo "worker $id finished: attempted=$mine new=$new fail=$fail"
}

echo "$(date '+%F %H:%M')  starting $NW DFAST workers at $CPU cpus each"
for w in $(seq 0 $((NW-1))); do worker "$w" & done
wait
echo "$(date '+%F %H:%M')  all workers done; merging"

MAIN="$R/dfast_counts.tsv"
[ -s "$MAIN" ] || printf 'accession\tdfast_disrupted\tdfast_frameshift\tdfast_internal_stop\tdfast_both\tdfast_partial\tdfast_cds\n' > "$MAIN"
cat "$R"/dfast_counts.w*.tsv >> "$MAIN"
rm -f "$R"/dfast_counts.w*.tsv
# de-duplicate on accession, keeping the header
python3 - "$MAIN" <<'EOF'
import csv, sys
p = sys.argv[1]
rows = list(csv.reader(open(p), delimiter="\t"))
hdr, seen, out = rows[0], set(), []
for r in rows[1:]:
    if not r or r[0] in seen or r[0] == "accession":
        continue
    seen.add(r[0]); out.append(r)
out.sort(key=lambda r: r[0])
with open(p, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(hdr); w.writerows(out)
print(f"  dfast_counts.tsv: {len(out)} unique accessions")
EOF
echo "$(date '+%F %H:%M')  DFAST COMPLETE"
