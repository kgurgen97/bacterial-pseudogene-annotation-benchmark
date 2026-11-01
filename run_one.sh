#!/usr/bin/env bash
g="$1"; D="$2"; O="$3"; P="$4"
acc=$(basename "$(dirname "$g")")
out="$O/$acc"
[ -s "$out/pf_pseudos.gff" ] && exit 0
mkdir -p "$out"
python "$P/tools/pseudofinder/pseudofinder.py" annotate \
  --genome "$g" --outprefix "$out/pf" --database "$D/sprot.dmnd" \
  --diamond --threads 4 > "$out/log.txt" 2>&1
