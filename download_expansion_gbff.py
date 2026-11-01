#!/usr/bin/env python3
"""Fetch GenBank flat files for the 160 expansion genomes.

Pseudofinder's `annotate` mode takes a GenBank file, not a FASTA ("Please provide your
genome file in the genbank format"), and every original corpus directory holds a
genomic.gbff beside its FASTA. The expansion download only took GENOME_FASTA, so without
this the Pseudofinder arm would simply have no data for the new strata.

Bakta and DFAST take the FASTA and are unaffected.
"""
import io, sys, time, urllib.parse, urllib.request, zipfile
from pathlib import Path
import pandas as pd

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
HDR = {"User-Agent": "research/1.0"}
P = Path.home() / "research/projects/n2-pseudogene-caller-bias"
OUT = Path.home() / "research/datasets/n2-genomes-expansion"

accs = pd.read_csv(P / "manifests/expansion_selected.tsv", sep="\t").accession.tolist()
todo = [a for a in accs if not (OUT / a / "genomic.gbff").exists()]
print(f"{len(accs)} accessions, {len(todo)} still need a gbff")

fail = []
for i, acc in enumerate(todo, 1):
    url = (f"{API}/genome/accession/{acc}/download?"
           + urllib.parse.urlencode({"include_annotation_type": "GENOME_GBFF"}))
    blob = None
    for a in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=HDR),
                                        timeout=600) as f:
                blob = f.read()
            break
        except Exception as e:
            if a == 3:
                print(f"  {acc}: {type(e).__name__}", file=sys.stderr)
            else:
                time.sleep(3 * (a + 1))
    if blob is None:
        fail.append(acc); continue
    try:
        z = zipfile.ZipFile(io.BytesIO(blob))
    except zipfile.BadZipFile:
        print(f"  {acc}: not a zip ({len(blob)} bytes)", file=sys.stderr)
        fail.append(acc); continue
    wrote = False
    for name in z.namelist():
        if name.endswith((".gbff", ".gbk", ".gb")):
            (OUT / acc).mkdir(parents=True, exist_ok=True)
            (OUT / acc / "genomic.gbff").write_bytes(z.read(name))
            wrote = True
            break
    if not wrote:
        print(f"  {acc}: no gbff inside; members {z.namelist()[:4]}", file=sys.stderr)
        fail.append(acc)
    if i % 20 == 0 or i == len(todo):
        have = sum(1 for a in accs if (OUT / a / "genomic.gbff").exists())
        print(f"  {i}/{len(todo)} requested, {have}/{len(accs)} on disk", flush=True)
    time.sleep(0.34)

have = [a for a in accs if (OUT / a / "genomic.gbff").exists()]
print(f"\ngbff present: {len(have)}/{len(accs)}")
if fail:
    print(f"failed ({len(fail)}): {fail[:15]}")
size = sum((OUT / a / "genomic.gbff").stat().st_size for a in have)
print(f"disk: {size/1e6:.0f} MB")

# a GenBank file Pseudofinder can actually read must carry CDS features
bad = []
for a in have[:len(have)]:
    p = OUT / a / "genomic.gbff"
    head = p.open("rb").read(400_000).decode("utf-8", "replace")
    if "     CDS " not in head:
        bad.append(a)
print(f"files with no CDS feature in the first 400 kB: {len(bad)}")
if bad:
    print(f"  {bad[:12]}")
