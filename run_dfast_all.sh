#!/usr/bin/env bash
# DFAST over all 363 genomes as the fourth caller. ~2m15s and 1.8 GB per genome.
#
# COUNTING RULE, fixed before this was run and not to be changed because of a result:
#
#   DFAST's pseudogene_summary.tsv has a `classification` column with four values:
#   partial, frameshift, internal_stop_codon, and "frameshift,internal_stop_codon".
#   On the probe genome those were 307 / 25 / 6 / 4 = 342 rows in total.
#
#   `partial` is NOT a pseudogene call. It marks a CDS with low coverage against its
#   reference or truncated at a contig edge - largely an assembly property. Counting it
#   would put DFAST at 342 on a genome where PGAP says 43 and Pseudofinder 194, which
#   would be an artefact of definition, not a finding.
#
#   So the primary DFAST count is the DISRUPTED classes only:
#       frameshift + internal_stop_codon + both      (35 on the probe genome)
#   and the `partial` count is recorded alongside it, never added to it.
#
# Note also: do not grep the GFF for "pseudo". It matches Pseudomonas, pseudomallei,
# pseudouridine and pseudopilin. The summary TSV is the only authoritative source.
#
# Resumable: a genome whose pseudogene_summary.tsv already exists is skipped.
set -u
source ~/miniforge3/etc/profile.d/conda.sh
conda activate callers

E=~/research/datasets/n2-genomes/extracted
OUT=~/research/datasets/n2-dfast
LOG=~/research/datasets/n2-dfast-logs
TSV=~/research/projects/n2-pseudogene-caller-bias/results/dfast_counts.tsv
mkdir -p "$OUT" "$LOG"

[ -s "$TSV" ] || printf 'accession\tdfast_disrupted\tdfast_frameshift\tdfast_internal_stop\tdfast_both\tdfast_partial\tdfast_cds\n' > "$TSV"

count() {  # $1 = summary tsv, $2 = exact classification value
  awk -F'\t' -v c="$2" 'NR>1 && $4==c {n++} END {print n+0}' "$1"
}

# -L is required: the expansion genomes are symlinked directories under $E and
# plain find does not descend into a symlink, silently skipping them.
mapfile -t FNAS < <(find -L "$E" -name '*.fna' | sort)
echo "genomes found: ${#FNAS[@]}"
cached=0; new=0; fail=0

for G in "${FNAS[@]}"; do
  ACC=$(basename "$(dirname "$G")")
  O="$OUT/$ACC"
  S="$O/pseudogene_summary.tsv"
  if [ -s "$S" ]; then cached=$((cached+1)); continue; fi
  if dfast -g "$G" -o "$O" --cpu 4 --force > "$LOG/$ACC.out" 2> "$LOG/$ACC.err"; then
    fs=$(count "$S" frameshift)
    is=$(count "$S" internal_stop_codon)
    bo=$(count "$S" 'frameshift,internal_stop_codon')
    pa=$(count "$S" partial)
    cd_=$(awk -F'\t' '/^Number of CDSs/{print $2}' "$O/statistics.txt")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$ACC" "$((fs+is+bo))" "$fs" "$is" "$bo" "$pa" "${cd_:-NA}" >> "$TSV"
    # keep the summary and statistics; drop the bulky sequence output
    rm -f "$O/protein.faa" "$O/cds.fna" "$O/rna.fna" "$O/genome.fna" \
          "$O/genome.embl" "$O/genome.gbk"
    rm -rf "$O/ddbj"
    new=$((new+1))
  else
    fail=$((fail+1)); echo "FAILED $ACC" >> "$OUT/failures.txt"
  fi
  if [ $(( (new + cached) % 25 )) -eq 0 ]; then
    echo "$(date +%H:%M)  cached=$cached new=$new fail=$fail  of ${#FNAS[@]}"
  fi
done

echo "FINISHED  cached=$cached new=$new failed=$fail"
echo "rows in $TSV : $(( $(wc -l < "$TSV") - 1 ))"
