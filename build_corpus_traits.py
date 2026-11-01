#!/usr/bin/env python3
"""Can an objective host-association label be attached to the 363 genomes?

Join the corpus to Madin et al. 2020 on tax_id, then species name, then genus, and
report coverage at each step. Also look at what the isolation_source vocabulary offers
and at doubling_h, which the literature links to effective population size.
"""
import collections
from pathlib import Path
import numpy as np, pandas as pd

P = Path.home() / "research/projects/n2-pseudogene-caller-bias"
T = Path.home() / "research/datasets/bacteria-traits"

sp = pd.read_csv(T / "condensed_species_NCBI.csv", dtype=str, low_memory=False)
print(f"Madin species table: {sp.shape}")
print("\n=== isolation_source vocabulary, top 30")
iso = sp.isolation_source.dropna()
for k, v in collections.Counter(iso).most_common(30):
    print(f"  {v:6d}  {k}")
print(f"  distinct values: {iso.nunique()}  coverage {len(iso)}/{len(sp)} = {len(iso)/len(sp):.1%}")

print("\n=== doubling_h coverage")
dh = pd.to_numeric(sp.doubling_h, errors="coerce").dropna()
print(f"  {len(dh)}/{len(sp)} = {len(dh)/len(sp):.1%}  median {dh.median():.2f} h  "
      f"range {dh.min():.2f}-{dh.max():.1f}")

# ---- build the corpus table
frames = []
for f in ("h10_manifest.tsv", "h12_manifest.tsv"):
    d = pd.read_csv(P / "manifests" / f, sep="\t", dtype=str)
    keep = [c for c in ("accession", "tax_id", "organism", "species", "genus",
                        "family", "order", "class", "phylum", "gc_percent",
                        "genome_size", "isolation_source") if c in d.columns]
    frames.append(d[keep].assign(src=f.split("_")[0]))
h11 = pd.read_csv(P / "manifests" / "h11_primary_genomes.tsv", sep="\t", dtype=str)
frames.append(h11[[c for c in ("accession", "genus", "family", "order",
                               "gc_percent", "genome_size_mb", "isolation_source")
                   if c in h11.columns]].assign(src="h11"))
cor = pd.concat(frames, ignore_index=True)
cor = cor.drop_duplicates("accession")
print(f"\n=== corpus: {len(cor)} accessions")
print(f"  with tax_id   : {cor.tax_id.notna().sum() if 'tax_id' in cor else 0}")
print(f"  with species  : {cor.species.notna().sum() if 'species' in cor else 0}")
print(f"  with genus    : {cor.genus.notna().sum()}")
print(f"  phylum split  : {dict(collections.Counter(cor.phylum.dropna())) if 'phylum' in cor else 'n/a'}")

# ---- join cascade
sp_by_tax = sp.set_index("species_tax_id")
sp_by_name = sp.drop_duplicates("species").set_index("species")
gen = (sp.assign(_n=1).groupby("genus")
         .agg(isolation_source=("isolation_source",
                                lambda s: s.dropna().mode().iat[0] if s.notna().any() else np.nan),
              doubling_h=("doubling_h", lambda s: pd.to_numeric(s, errors="coerce").median()),
              n_species=("_n", "sum")))

def lookup(row):
    tid = row.get("tax_id")
    if pd.notna(tid) and str(tid) in sp_by_tax.index:
        r = sp_by_tax.loc[str(tid)]
        if isinstance(r, pd.DataFrame): r = r.iloc[0]
        return r.isolation_source, r.doubling_h, "tax_id"
    s = row.get("species")
    if pd.notna(s) and s in sp_by_name.index:
        r = sp_by_name.loc[s]
        return r.isolation_source, r.doubling_h, "species"
    g = row.get("genus")
    if pd.notna(g) and g in gen.index:
        r = gen.loc[g]
        return r.isolation_source, r.doubling_h, "genus"
    return np.nan, np.nan, "none"

res = cor.apply(lookup, axis=1, result_type="expand")
cor["iso_madin"], cor["doubling_h"], cor["join_level"] = res[0], res[1], res[2]
print(f"\n=== join level")
for k, v in collections.Counter(cor.join_level).most_common():
    print(f"  {k:8s} {v:4d}")
print(f"  isolation_source obtained for {cor.iso_madin.notna().sum()}/{len(cor)}")
print(f"  doubling_h obtained for       {pd.to_numeric(cor.doubling_h, errors='coerce').notna().sum()}/{len(cor)}")

print("\n=== isolation_source values in OUR corpus")
for k, v in collections.Counter(cor.iso_madin.dropna()).most_common(25):
    print(f"  {v:4d}  {k}")

cor.to_csv(P / "results/corpus_traits.tsv", sep="\t", index=False)
print(f"\nwrote results/corpus_traits.tsv")
