#!/usr/bin/env bash
# Reproduces the retrieval behind SIMULATION_NOVELTY_AUDIT.md.
#
# Every bibliographic detail and every "appears N times" count in that document
# comes from this script. It writes a short log to logs/lit_audit.log and leaves
# the retrieved full texts in a scratch directory outside the repository, since
# they are third-party copyrighted material and are not redistributed here.
#
# The publisher sites for the two BMC papers redirect to a Springer identity
# provider and refuse automated access. Europe PMC is used instead, where all
# three primary sources are open access.
set -eu
E=https://www.ebi.ac.uk/europepmc/webservices/rest
OUT=${OUT:-$HOME/lit_audit}
LOG=$(cd "$(dirname "$0")" && pwd)/logs/lit_audit.log
mkdir -p "$OUT" "$(dirname "$LOG")"
exec > >(tee "$LOG") 2>&1

echo "retrieval for SIMULATION_NOVELTY_AUDIT.md"
echo "date $(date -u +%Y-%m-%dT%H:%M:%SZ)  scratch $OUT"
echo

# ---------------------------------------------------------------- records
echo "=== bibliographic records, from Europe PMC core"
for doi in 10.1186/s13059-026-04262-0 10.1093/molbev/msac153 \
           10.1186/s12864-024-10137-0 10.1099/mgen.0.000685 \
           10.1093/bioinformatics/btx713; do
  curl -sfG "$E/search" --data-urlencode "query=DOI:\"$doi\"" \
       --data-urlencode "resultType=core" --data-urlencode "format=json" \
       -o "$OUT/rec.json"
  python3 -I - "$OUT/rec.json" <<'PY'
import json, sys, pathlib
d = json.loads(pathlib.Path(sys.argv[1]).read_text())
for x in d["resultList"]["result"]:
    ji = x.get("journalInfo", {})
    print(f"  {x.get('doi')}")
    print(f"    {x.get('title')}")
    print(f"    {x.get('authorString')}")
    print(f"    {ji.get('journal',{}).get('medlineAbbreviation')} "
          f"{ji.get('volume')}({ji.get('issue')}) {x.get('pubYear')} "
          f"art {x.get('pageInfo')}  PMID {x.get('pmid')}  {x.get('pmcid')}  "
          f"OA {x.get('isOpenAccess')}")
PY
done

# ---------------------------------------------------------------- full texts
echo
echo "=== full texts, and the term counts quoted in the audit"
for pmc in PMC13548610 PMC11017572; do
  curl -sfL "$E/$pmc/fullTextXML" -o "$OUT/$pmc.xml"
  echo "  $pmc  $(wc -c < "$OUT/$pmc.xml") bytes"
  python3 -I - "$OUT/$pmc.xml" <<'PY'
import re, sys, pathlib
t = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ",
      pathlib.Path(sys.argv[1]).read_text(errors="replace")))
for kw in ("pseudogene", "Pseudofinder", "DFAST", "Bakta",
           "ground truth", "sensitivit", "precision", "frameshift"):
    print(f"    {kw!r:16s} {len(re.findall(kw, t, re.I)):4d}")
PY
done
echo "  PMC9336565 (Pseudofinder) was read for Results, Methods, figure captions"
echo "    and Conclusions; Supplementary Fig. S1-S3 and the project wiki were not."

# ---------------------------------------------------------------- prior art
echo
echo "=== prior-art searches"
search() {
  curl -sfG "$E/search" --data-urlencode "query=$2" --data-urlencode "format=json" \
       --data-urlencode "pageSize=25" --data-urlencode "sort=CITED desc" -o "$OUT/q.json"
  python3 -I - "$OUT/q.json" "$1" <<'PY'
import json, sys, pathlib
d = json.loads(pathlib.Path(sys.argv[1]).read_text())
print(f"  [{d['hitCount']:>5} hits] {sys.argv[2]}")
for x in d["resultList"]["result"][:6]:
    j = (x.get("journalInfo", {}).get("journal", {}) or {}).get("medlineAbbreviation", "?")
    print(f"      {x.get('pubYear')} {j:16s} cited {x.get('citedByCount',0):4}  "
          f"{(x.get('title') or '')[:92]}")
PY
}
search "pseudogene benchmarking in bacteria" \
  'TITLE_ABS:("pseudogene" AND (benchmark OR benchmarking OR "performance evaluation")) AND (bacteri* OR prokaryot*)'
search "in silico pseudogenization" \
  'TITLE_ABS:("pseudogene*" AND (simulat* OR "in silico") AND (frameshift OR "stop codon" OR disrupt* OR inactivat*))'
search "prokaryotic pseudogene detection tools, 2018-2026" \
  'TITLE_ABS:("pseudogene" AND (detect* OR identif* OR annotat*) AND (prokaryot* OR bacterial)) AND PUB_YEAR:[2018 TO 2026]'
search "annotation pipeline disagreement" \
  'TITLE_ABS:(("genome annotation" OR "annotation pipeline") AND (discrepan* OR disagree* OR inconsisten* OR concordance))'

echo
echo "done. full texts left in $OUT and deliberately not committed."
