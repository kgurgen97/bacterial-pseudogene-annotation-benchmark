#!/usr/bin/env python3
"""Select the expansion genomes exactly as CORPUS_EXPANSION_PLAN.md specifies.

Plan sha256: 1b1ee35c19bd66388ba7b952de41abab840e08ebbd4656bf08e47aef83e8c5fa

Deterministic: no randomness anywhere, so the selected set is reproducible from the plan.
Writes manifests/expansion_candidates.tsv (everything seen) and
manifests/expansion_selected.tsv (what will be downloaded), plus the reasons for
every exclusion.
"""
import collections, json, re, sys, time, urllib.parse, urllib.request
from pathlib import Path
import pandas as pd

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
HDR = {"User-Agent": "research/1.0", "Accept": "application/json"}
P = Path.home() / "research/projects/n2-pseudogene-caller-bias"
M = P / "manifests"

STRATA = {
    # S1 and S3 lists extended per the 2026-10-22 amendment in the plan (supply
    # shortfall; no pseudogene value had been computed for any new genome).
    "S1_extreme_obligate": (25, 5, ["Buchnera", "Carsonella", "Candidatus Sulcia",
                                    "Candidatus Tremblaya", "Candidatus Zinderia",
                                    "Candidatus Portiera", "Candidatus Blochmannia",
                                    "Wigglesworthia",
                                    "Candidatus Moranella", "Candidatus Profftella",
                                    "Candidatus Evansia", "Candidatus Vidania",
                                    "Candidatus Purcelliella", "Candidatus Uzinura",
                                    "Candidatus Walczuchella", "Candidatus Westeberhardia",
                                    "Candidatus Nardonella",
                                    "Candidatus Kinetoplastibacterium"]),
    "S2_obligate_intracellular": (35, 5, ["Rickettsia", "Orientia", "Anaplasma", "Ehrlichia",
                                          "Wolbachia", "Chlamydia", "Coxiella", "Mycoplasma",
                                          "Mycoplasmoides", "Ureaplasma"]),
    "S3_recent_restriction": (30, 5, ["Serratia symbiotica", "Sodalis", "Mycobacterium leprae",
                                      "Yersinia pestis", "Shigella", "Bordetella pertussis",
                                      "Treponema pallidum", "Francisella tularensis",
                                      "Burkholderia mallei",
                                      "Salmonella enterica subsp. enterica serovar Typhi",
                                      "Mycobacterium ulcerans", "Bordetella parapertussis",
                                      "Bordetella holmesii",
                                      "Salmonella enterica subsp. enterica serovar Paratyphi A",
                                      "Arsenophonus", "Hamiltonella", "Regiella",
                                      "Candidatus Riesia", "Baumannia"]),
    "S4_free_living_breadth": (70, 5, ["Bacillota", "Bacteroidota", "Cyanobacteriota",
                                       "Spirochaetota", "Campylobacterota", "Thermotogota",
                                       "Deinococcota", "Fusobacteriota", "Verrucomicrobiota",
                                       "Planctomycetota", "Acidobacteriota", "Chloroflexota",
                                       "Aquificota", "Nitrospirota", "Synergistota",
                                       "Thermodesulfobacteriota"]),
}
RANK = {"reference genome": 0, "representative genome": 1, "": 2, None: 2}


def fetch(taxon, want=4000):
    """All complete RefSeq reports for a taxon, paged."""
    out, token = [], None
    while len(out) < want:
        q = {"filters.assembly_level": "complete_genome",
             "filters.assembly_source": "refseq",
             "filters.exclude_atypical": "true",
             "page_size": "1000"}
        if token:
            q["page_token"] = token
        url = f"{API}/genome/taxon/{urllib.parse.quote(taxon)}/dataset_report?" + \
              urllib.parse.urlencode(q)
        for attempt in range(4):
            try:
                with urllib.request.urlopen(urllib.request.Request(url, headers=HDR),
                                            timeout=180) as f:
                    d = json.load(f)
                break
            except Exception as e:
                if attempt == 3:
                    print(f"    {taxon}: giving up ({type(e).__name__})", file=sys.stderr)
                    return out
                time.sleep(2 * (attempt + 1))
        out += d.get("reports", [])
        token = d.get("next_page_token")
        if not token:
            break
        time.sleep(0.34)
    return out


def species_of(name):
    """'Buchnera aphidicola str. APS' -> 'Buchnera aphidicola'; keeps Candidatus binomials."""
    t = name.replace("[", "").replace("]", "").split()
    if not t:
        return name
    if t[0] == "Candidatus":
        return " ".join(t[:3])
    return " ".join(t[:2])


def flat(r):
    org = r.get("organism", {}) or {}
    ai = r.get("assembly_info", {}) or {}
    st = r.get("assembly_stats", {}) or {}
    ck = r.get("checkm_info", {}) or {}
    nm = org.get("organism_name", "") or ""
    return dict(accession=r.get("accession", ""), organism=nm,
                species=species_of(nm),
                genus=species_of(nm).split()[1] if nm.startswith("Candidatus ")
                      else (nm.split()[0] if nm else ""),
                tax_id=org.get("tax_id"),
                refseq_category=ai.get("refseq_category") or "",
                level=ai.get("assembly_level") or "",
                size_bp=st.get("total_sequence_length"),
                n_contigs=st.get("number_of_contigs"),
                gc=st.get("gc_percent"),
                completeness=ck.get("completeness"),
                contamination=ck.get("contamination"))


# The ORIGINAL corpus manifests only, named explicitly. A glob over manifests/*.tsv
# would also pick up expansion_candidates.tsv written by a previous run of this very
# script, which then counts every candidate as "already in the corpus" - that happened,
# and it silently reduced a 140-genome selection to 27.
ORIGINAL = ["accession_asm.tsv", "h10_manifest.tsv", "h11_primary_genomes.tsv",
            "h12_manifest.tsv"]
existing = set()
for name in ORIGINAL:
    f = M / name
    if not f.exists():
        continue
    d = pd.read_csv(f, sep="\t", dtype=str)
    c = next((x for x in d.columns if "accession" in x.lower()), None)
    if c:
        existing |= set(d[c].dropna())
t = M / "all_accessions.txt"
if t.exists():
    existing |= {l.strip() for l in t.read_text().split() if l.strip()}
print(f"already in the corpus: {len(existing)} accessions "
      f"(from {', '.join(ORIGINAL)} + all_accessions.txt)\n")
if len(existing) != 363:
    print(f"  WARNING: expected 363, got {len(existing)}", file=sys.stderr)

cand_rows, sel_rows, drops = [], [], collections.Counter()

for stratum, (target, cap_genus, taxa) in STRATA.items():
    print(f"=== {stratum}  target {target}, cap {cap_genus} per genus")
    pool = []
    for taxon in taxa:
        rs = fetch(taxon)
        for r in rs:
            f_ = flat(r)
            f_["stratum"] = stratum
            f_["query_taxon"] = taxon
            pool.append(f_)
        print(f"    {len(rs):5d} reports  {taxon}")
        time.sleep(0.34)
    df = pd.DataFrame(pool)
    if df.empty:
        continue
    df = df.drop_duplicates("accession")
    # the API returns these as strings
    for c in ("size_bp", "n_contigs", "gc", "completeness", "contamination"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    cand_rows.append(df)

    before = len(df)
    df = df[~df.accession.isin(existing)]
    drops[f"{stratum}: already in corpus"] += before - len(df)
    before = len(df)
    df = df[df.size_bp.notna() & (df.size_bp > 0)]
    drops[f"{stratum}: no size"] += before - len(df)
    before = len(df)
    df = df[(df.contamination.isna()) | (df.contamination <= 5)]
    drops[f"{stratum}: contamination > 5"] += before - len(df)

    # one genome per species, preferring reference > representative > smallest accession
    df["rank"] = df.refseq_category.map(lambda v: RANK.get(v, 2))
    df = df.sort_values(["species", "rank", "accession"])
    df = df.drop_duplicates("species", keep="first")

    # cap per genus, then take the target in a deterministic order
    df = df.sort_values(["rank", "genus", "species", "accession"])
    keep, per_genus = [], collections.Counter()
    for _, r in df.iterrows():
        if per_genus[r.genus] >= cap_genus:
            continue
        per_genus[r.genus] += 1
        keep.append(r)
        if len(keep) >= target:
            break
    k = pd.DataFrame(keep)
    print(f"    -> {len(k)} selected from {len(df)} species "
          f"({k.genus.nunique() if len(k) else 0} genera)")
    if len(k):
        print(f"       size range {k.size_bp.min()/1e6:.2f}-{k.size_bp.max()/1e6:.2f} Mb, "
              f"median {k.size_bp.median()/1e6:.2f} Mb")
    sel_rows.append(k)
    print()

cand = pd.concat(cand_rows, ignore_index=True).drop_duplicates("accession")
sel = pd.concat(sel_rows, ignore_index=True)
cand.to_csv(M / "expansion_candidates.tsv", sep="\t", index=False)
sel.to_csv(M / "expansion_selected.tsv", sep="\t", index=False)

print("=" * 78)
print(f"candidates seen : {len(cand)}")
print(f"SELECTED        : {len(sel)}")
print(f"\nby stratum:")
for s, g in sel.groupby("stratum"):
    print(f"  {s:28s} {len(g):4d}  size {g.size_bp.min()/1e6:5.2f}-"
          f"{g.size_bp.max()/1e6:5.2f} Mb (median {g.size_bp.median()/1e6:.2f})")
print(f"\nexclusions:")
for k, v in drops.most_common():
    if v:
        print(f"  {k}: {v}")
print(f"\ndistinct species {sel.species.nunique()}, genera {sel.genus.nunique()}")
print(f"overlap with the existing corpus: {len(set(sel.accession) & existing)} (must be 0)")
print(f"\nwrote {M/'expansion_selected.tsv'}")
