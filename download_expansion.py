#!/usr/bin/env python3
"""Download the 160 selected genomes and collect their PGAP counts from the same API.

PGAP's pseudogene and protein-coding counts come from NCBI's own annotation_info
gene_counts, which is the same RefSeq/PGAP annotation the original corpus used - so the
PGAP arm needs no local rerun, only a lookup.
"""
import io, json, sys, time, urllib.parse, urllib.request, zipfile
from pathlib import Path
import pandas as pd

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
HDR = {"User-Agent": "research/1.0"}
P = Path.home() / "research/projects/n2-pseudogene-caller-bias"
OUT = Path.home() / "research/datasets/n2-genomes-expansion"
OUT.mkdir(parents=True, exist_ok=True)

sel = pd.read_csv(P / "manifests/expansion_selected.tsv", sep="\t")
accs = sel.accession.tolist()
print(f"{len(accs)} accessions to fetch")


def get(url, tries=4, timeout=600):
    for a in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=HDR),
                                        timeout=timeout) as f:
                return f.read()
        except Exception as e:
            if a == tries - 1:
                print(f"    failed: {type(e).__name__}", file=sys.stderr)
                return None
            time.sleep(3 * (a + 1))


# ---------------------------------------------------------------- FASTA
print("\n=== genomic FASTA, one accession per request")
# The batch /genome/download endpoint answers 405 to GET - it needs POST. The
# per-accession form works on GET and returns a zip, so use that.
have = {d.name for d in OUT.iterdir() if d.is_dir() and (d / "genomic.fna").exists()}
print(f"  already on disk: {len(have)}")
todo = [a for a in accs if a not in have]
for n_i, acc in enumerate(todo, 1):
    url = (f"{API}/genome/accession/{acc}/download?"
           + urllib.parse.urlencode({"include_annotation_type": "GENOME_FASTA"}))
    blob = get(url, timeout=300)
    if blob is None:
        continue
    try:
        z = zipfile.ZipFile(io.BytesIO(blob))
    except zipfile.BadZipFile:
        print(f"    {acc}: not a zip ({len(blob)} bytes)", file=sys.stderr)
        continue
    wrote = False
    for name in z.namelist():
        if name.endswith((".fna", ".fa", ".fasta")):
            d = OUT / acc
            d.mkdir(exist_ok=True)
            (d / "genomic.fna").write_bytes(z.read(name))
            wrote = True
            break
    if not wrote:
        print(f"    {acc}: no FASTA inside; members {z.namelist()[:4]}", file=sys.stderr)
    if n_i % 20 == 0 or n_i == len(todo):
        done = sum(1 for a in accs if (OUT / a / "genomic.fna").exists())
        print(f"  {n_i}/{len(todo)} requested, {done}/{len(accs)} on disk", flush=True)
    time.sleep(0.34)

ok = [a for a in accs if (OUT / a / "genomic.fna").exists()]
missing = [a for a in accs if a not in ok]
print(f"\nFASTA present: {len(ok)}/{len(accs)}")
if missing:
    print(f"  missing: {missing[:12]}{' ...' if len(missing) > 12 else ''}")

# ---------------------------------------------------------------- PGAP counts
print("\n=== PGAP gene counts from annotation_info")
rows = []
for i in range(0, len(accs), 20):
    chunk = accs[i:i + 20]
    url = (f"{API}/genome/accession/{','.join(chunk)}/dataset_report?"
           + urllib.parse.urlencode({"page_size": "100"}))
    blob = get(url, timeout=180)
    if blob is None:
        continue
    for r in json.loads(blob).get("reports", []):
        ai = r.get("annotation_info", {}) or {}
        gc = (ai.get("stats", {}) or {}).get("gene_counts", {}) or {}
        st = r.get("assembly_stats", {}) or {}
        # NCBI OMITS the pseudogene key when the count is zero. Verified on
        # GCF_000287275.1 (Carsonella), GCF_000340795.1 and GCF_053760525.1: each has
        # total == protein_coding + non_coding exactly, so there are genuinely no
        # pseudogenes. Treating those as missing would discard the strongest evidence
        # for C3, since a zero in an extremely reduced endosymbiont is the literature's
        # claim about the reduced endpoint. Infer 0 only when the arithmetic closes.
        pseudo = gc.get("pseudogene")
        inferred = False
        if pseudo is None:
            tot, pc, nc = gc.get("total"), gc.get("protein_coding"), gc.get("non_coding")
            try:
                if None not in (tot, pc, nc) and int(tot) == int(pc) + int(nc):
                    pseudo, inferred = 0, True
            except (TypeError, ValueError):
                pass
        rows.append(dict(accession=r.get("accession"),
                         organism=(r.get("organism", {}) or {}).get("organism_name"),
                         pgap_pseudo=pseudo,
                         pgap_pseudo_inferred_zero=inferred,
                         pgap_cds=gc.get("protein_coding"),
                         pgap_total=gc.get("total"),
                         annotation_name=ai.get("name"),
                         annotation_date=ai.get("release_date"),
                         size_bp=st.get("total_sequence_length"),
                         gc_percent=st.get("gc_percent")))
    time.sleep(0.4)
    print(f"  {min(i+20, len(accs))}/{len(accs)}", flush=True)

ann = pd.DataFrame(rows).drop_duplicates("accession")
for c in ("pgap_pseudo", "pgap_cds", "pgap_total", "size_bp", "gc_percent"):
    ann[c] = pd.to_numeric(ann[c], errors="coerce")
print(f"\nannotation rows: {len(ann)}")
print(f"  with a pseudogene count: {ann.pgap_pseudo.notna().sum()}")
print(f"  with a protein_coding count: {ann.pgap_cds.notna().sum()}")
miss = ann[ann.pgap_pseudo.isna() | ann.pgap_cds.isna()]
if len(miss):
    print(f"  incomplete annotation for {len(miss)}:")
    for _, r in miss.head(12).iterrows():
        print(f"    {r.accession}  {str(r.organism)[:46]}")

ann = ann.merge(sel[["accession", "stratum", "species", "genus", "query_taxon"]],
                on="accession", how="left")
ann["pgap_frac"] = ann.pgap_pseudo / (ann.pgap_pseudo + ann.pgap_cds)
ann.to_csv(P / "manifests/expansion_pgap.tsv", sep="\t", index=False)

print("\n=== PGAP fraction by stratum (first look at the gradient)")
for st, g in ann.groupby("stratum"):
    v = g.pgap_frac.dropna()
    if len(v):
        print(f"  {st:28s} n={len(v):3d}  median {v.median():.4f}  "
              f"IQR {v.quantile(.25):.4f}-{v.quantile(.75):.4f}")

tot = sum(f.stat().st_size for f in OUT.rglob("genomic.fna"))
print(f"\ndisk: {tot/1e6:.0f} MB in {OUT}")
print(f"wrote manifests/expansion_pgap.tsv")
