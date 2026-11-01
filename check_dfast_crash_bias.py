#!/usr/bin/env python3
"""Were the 21 DFAST crashes a biased loss?

The crash fired only where an rRNA/tRNA overlapped an HMM-named CDS, which is not a
random property of a genome. If the 21 differ systematically from the other 342, then
dropping them - the easy option - would have biased the caller comparison, and patching
was necessary rather than tidy.
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

P = Path.home() / "research/projects/n2-pseudogene-caller-bias"
RECOVERED = """GCF_000015025.1 GCF_000024025.1 GCF_000024865.1 GCF_000166135.1
GCF_000819565.1 GCF_001402875.1 GCF_001941625.1 GCF_002096055.1 GCF_002240355.1
GCF_003568845.1 GCF_009299385.1 GCF_019397265.1 GCF_029814805.1 GCF_038502305.1
GCF_042926695.1 GCF_043575855.1 GCF_050942415.1 GCF_053646815.1 GCF_054052205.1
GCF_900637955.1 GCF_948472415.1""".split()

df = pd.read_csv(P / "results/dfast_counts.tsv", sep="\t")
pg = pd.read_csv(P / "results/pseudogene_counts.tsv", sep="\t")
cor = pd.read_csv(P / "results/corpus_traits.tsv", sep="\t", dtype=str)
size = pd.to_numeric(cor.genome_size, errors="coerce")
mb = pd.to_numeric(cor.get("genome_size_mb"), errors="coerce") * 1e6
cor["size_bp"] = size.fillna(mb)
m = df.merge(pg[["accession", "pgap_cds", "pgap_pseudo", "pf_pseudo"]], on="accession") \
      .merge(cor[["accession", "size_bp", "phylum"]], on="accession", how="left")
m["grp"] = np.where(m.accession.isin(RECOVERED), "crashed", "ran first time")
print(f"merged {len(m)} genomes; crashed {int((m.grp=='crashed').sum())}")

for col, label in (("dfast_cds", "DFAST CDS count"),
                   ("size_bp", "genome size (bp)"),
                   ("dfast_disrupted", "DFAST disrupted"),
                   ("dfast_partial", "DFAST partial"),
                   ("pgap_pseudo", "PGAP pseudogenes")):
    a = pd.to_numeric(m.loc[m.grp == "crashed", col], errors="coerce").dropna()
    b = pd.to_numeric(m.loc[m.grp != "crashed", col], errors="coerce").dropna()
    if len(a) < 3 or len(b) < 3:
        continue
    u = stats.mannwhitneyu(a, b, alternative="two-sided")
    print(f"\n  {label}")
    print(f"    crashed        n={len(a):3d}  median {a.median():12,.0f}")
    print(f"    ran first time n={len(b):3d}  median {b.median():12,.0f}")
    print(f"    Mann-Whitney p = {u.pvalue:.4g}"
          f"   ratio of medians {a.median()/b.median():.2f}")

print("\n  by phylum:")
print(pd.crosstab(m.phylum.fillna("unknown"), m.grp))

print("\n=== so, would dropping them have mattered?")
a = pd.to_numeric(m.loc[m.grp == 'crashed', 'dfast_cds'], errors='coerce').dropna()
b = pd.to_numeric(m.loc[m.grp != 'crashed', 'dfast_cds'], errors='coerce').dropna()
u = stats.mannwhitneyu(a, b, alternative='two-sided')
if u.pvalue < 0.05:
    print(f"  YES - the crashed genomes are systematically different in size "
          f"(p = {u.pvalue:.3g}), so dropping 6% of the corpus would have been a "
          f"biased loss, not a random one.")
else:
    print(f"  No detectable difference (p = {u.pvalue:.3g}); dropping them would "
          f"probably have been harmless, but recovering them costs nothing now.")
