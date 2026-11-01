#!/usr/bin/env bash
# DFAST first, parallelised across workers; then Bakta with the machine to itself.
#
# Measured on this box:
#   both callers at once      Bakta 14m25s, DFAST 5m56s per genome   (~84 h total)
#   one at a time             Bakta 6m10s,  DFAST 2m14s              (~49 h total)
#   DFAST across 4 workers    uses ~16 of 20 cores, 2 GB             (~40 h total)
# Bakta never shares: its DIAMOND claims all 20 cores and its 46 GB database already
# exceeds the 23 GB of RAM, so a second job steals cores and page cache at once.
#
# Safe to interrupt anywhere: both runners skip a genome only when its finished summary
# file exists, and that file is written last.
set -u
P=~/research/projects/n2-pseudogene-caller-bias

echo "$(date '+%F %H:%M')  DFAST, parallel workers"
bash "$P/run_dfast_parallel.sh"

echo "$(date '+%F %H:%M')  DFAST done; Bakta now has the machine to itself"
bash "$P/run_bakta_all.sh"

echo "$(date '+%F %H:%M')  BOTH CALLERS DONE"
printf 'bakta rows: %s\ndfast rows: %s\n' \
  "$(( $(wc -l < "$P/results/bakta_counts.tsv") - 1 ))" \
  "$(( $(wc -l < "$P/results/dfast_counts.tsv") - 1 ))"
