#!/usr/bin/env bash
# Bring all four callers up to the full 523-genome corpus, cheapest first.
#
# PGAP needs nothing: its counts for the 160 additions come from NCBI's own
# annotation_info gene_counts (manifests/expansion_pgap.tsv), which is the same
# RefSeq/PGAP annotation the original 363 used.
#
# Order is deliberate. DFAST (~2.3 min/genome across 4 workers) and the GenBank
# download run together because one is CPU-bound and the other network-bound. Then
# Pseudofinder, which needs those GenBank files. Bakta goes last and alone: it claims
# all 20 cores for DIAMOND and its 46 GB database already exceeds RAM, so sharing makes
# everything 2.4x slower.
#
# Every stage is resumable, so an interruption costs at most the genome in flight.
set -u
P=~/research/projects/n2-pseudogene-caller-bias
R=$P/results
EXPD=~/research/datasets/n2-genomes-expansion
D=~/research/datasets/n2-genomes
L=~/research/datasets

stamp() { date '+%F %H:%M'; }
rows() { echo $(( $(wc -l < "$1") - 1 )); }

echo "$(stamp)  START.  dfast=$(rows $R/dfast_counts.tsv) bakta=$(rows $R/bakta_counts.tsv) of 523"

echo "$(stamp)  [1/4] GenBank download for the 160 additions, in the background"
source ~/miniforge3/etc/profile.d/conda.sh
conda activate hpwh
python "$P/download_expansion_gbff.py" > "$L/gbff_download.log" 2>&1 &
GBFF_PID=$!

echo "$(stamp)  [2/4] DFAST on the additions"
bash "$P/run_dfast_parallel.sh" > "$L/dfast_expansion.log" 2>&1
echo "$(stamp)        dfast now $(rows $R/dfast_counts.tsv) of 523"

echo "$(stamp)  waiting for the GenBank download"
wait $GBFF_PID || true
tail -3 "$L/gbff_download.log"

echo "$(stamp)  [3/4] Pseudofinder on the additions"
conda activate pseudofinder
new=0; skip=0; fail=0
for g in "$EXPD"/GCF_*/genomic.gbff; do
  [ -s "$g" ] || continue
  acc=$(basename "$(dirname "$g")")
  if [ -s "$L/n2-pseudofinder/$acc/pf_pseudos.gff" ]; then skip=$((skip+1)); continue; fi
  if bash "$P/run_one.sh" "$g" "$D" "$L/n2-pseudofinder" "$P"; then
    new=$((new+1))
  else
    fail=$((fail+1)); echo "PF FAILED $acc" >> "$L/n2-pseudofinder/failures.txt"
  fi
  if [ $(( (new+skip) % 20 )) -eq 0 ]; then
    echo "$(stamp)        pseudofinder new=$new cached=$skip fail=$fail"
  fi
done
echo "$(stamp)        pseudofinder done: new=$new cached=$skip fail=$fail"
echo "          total pf dirs: $(ls "$L/n2-pseudofinder" | grep -c '^GCF_')"

echo "$(stamp)  [4/4] Bakta on everything outstanding, with the machine to itself"
bash "$P/run_bakta_all.sh" > "$L/bakta_expansion.log" 2>&1
echo "$(stamp)        bakta now $(rows $R/bakta_counts.tsv) of 523"

echo "$(stamp)  ALL FOUR CALLERS DONE"
