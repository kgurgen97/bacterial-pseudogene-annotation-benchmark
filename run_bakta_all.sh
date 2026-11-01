#!/usr/bin/env bash
# Bakta over all 363 genomes, one at a time with 20 threads.
#
# The pseudogene count is taken from Bakta's own summary line in <prefix>.txt.
# Do NOT grep the TSV for "pseudo": that matches products such as
# "23S rRNA (pseudouridine(1915)-N(3))-methyltransferase", which gave 21 instead of 7
# on the probe genome. The authoritative places are the .txt summary, the JSON
# features carrying a "pseudogene" key, and the GFF3 "pseudogene=" attribute.
#
# Resumable: a genome whose out/<acc>/<acc>.txt already exists is skipped.
set -u
source ~/miniforge3/etc/profile.d/conda.sh
conda activate callers

DB=~/research/datasets/bakta-db/db
E=~/research/datasets/n2-genomes/extracted
OUT=~/research/datasets/n2-bakta
TSV=~/research/projects/n2-pseudogene-caller-bias/results/bakta_counts.tsv
mkdir -p "$OUT"

[ -s "$TSV" ] || printf 'accession\tbakta_pseudo\tbakta_cds\tbakta_sorf\tbakta_length\tbakta_coding_density\n' > "$TSV"

# -L is required: the 160 expansion genomes are symlinked directories under $E, and
# plain find does not descend into a symlink, so it would silently skip all of them.
mapfile -t FNAS < <(find -L "$E" -name '*.fna' | sort)
echo "genomes found: ${#FNAS[@]}"
done_n=0; new_n=0; fail_n=0

for G in "${FNAS[@]}"; do
  ACC=$(basename "$(dirname "$G")")
  O="$OUT/$ACC"
  if [ -s "$O/$ACC.txt" ]; then done_n=$((done_n+1)); continue; fi
  mkdir -p "$O"
  if bakta --db "$DB" --output "$O" --prefix "$ACC" --threads 20 \
           --skip-plot --force "$G" > "$O/run.log" 2>&1; then
    ps=$(awk -F': *' '/^pseudogenes:/{print $2}'      "$O/$ACC.txt")
    cd_=$(awk -F': *' '/^CDSs:/{print $2}'            "$O/$ACC.txt")
    so=$(awk -F': *' '/^sORFs:/{print $2}'            "$O/$ACC.txt")
    ln=$(awk -F': *' '/^Length:/{print $2}'           "$O/$ACC.txt")
    de=$(awk -F': *' '/^coding density:/{print $2}'   "$O/$ACC.txt")
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$ACC" "${ps:-NA}" "${cd_:-NA}" "${so:-NA}" "${ln:-NA}" "${de:-NA}" >> "$TSV"
    # keep the small, informative outputs; drop the bulky sequence files
    rm -f "$O/$ACC.faa" "$O/$ACC.ffn" "$O/$ACC.fna" "$O/$ACC.embl" \
          "$O/$ACC.gbff" "$O/$ACC.hypotheticals.faa"
    new_n=$((new_n+1))
  else
    fail_n=$((fail_n+1))
    echo "FAILED $ACC" >> "$OUT/failures.txt"
  fi
  if [ $(( (new_n + done_n) % 10 )) -eq 0 ]; then
    echo "$(date +%H:%M)  done=$done_n new=$new_n fail=$fail_n  of ${#FNAS[@]}"
  fi
done

echo "FINISHED  cached=$done_n new=$new_n failed=$fail_n"
echo "counts in $TSV : $(( $(wc -l < "$TSV") - 1 )) rows"
