#!/usr/bin/env python3
"""Does selection agree with the caller? As DNDS_PLAN.md specifies.

A pseudogene is a sequence no longer held by purifying selection, which leaves a
signature no annotation pipeline produces. Our Pseudofinder runs did not use
--reference, so they computed no dN/dS and the measure is independent of all four
callers; that was checked before the plan was written.

Run with --probe N to do the first N target genomes only.
"""
import argparse
import collections
import gzip
import hashlib
import json
import pathlib
import re
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
from Bio import Align, SeqIO
from Bio.Seq import Seq
from scipy import stats

P = pathlib.Path.home() / "research/projects/n2-pseudogene-caller-bias"
D = pathlib.Path.home() / "research/datasets"
R = P / "results"

PLAN = P / "DNDS_PLAN.md"
print(f"plan sha256: {hashlib.sha256(PLAN.read_bytes()).hexdigest()[:16]}…")

MIN_CODONS = 30          # plan section 2
MIN_IDENT = 50.0         # percent
MIN_COV = 0.50
CALLERS = ["PGAP", "Pseudofinder", "Bakta", "DFAST"]
RNG = np.random.default_rng(20261028)

# ----------------------------------------------------------------- NG86
SYN_CACHE: dict = {}


def syn_sites(codon):
    """Expected synonymous sites in a codon, Nei-Gojobori."""
    if codon in SYN_CACHE:
        return SYN_CACHE[codon]
    if len(codon) != 3 or any(c not in "ACGT" for c in codon):
        SYN_CACHE[codon] = None
        return None
    aa = str(Seq(codon).translate())
    if aa == "*":
        SYN_CACHE[codon] = None
        return None
    s = 0.0
    for i in range(3):
        for b in "ACGT":
            if b == codon[i]:
                continue
            alt = codon[:i] + b + codon[i + 1:]
            if str(Seq(alt).translate()) == aa:
                s += 1.0 / 3.0
    SYN_CACHE[codon] = s
    return s


def pairwise_diffs(c1, c2):
    """Synonymous and non-synonymous differences between two codons.

    Averaged over all shortest mutational paths, which is what NG86 specifies
    when the codons differ at more than one position.
    """
    pos = [i for i in range(3) if c1[i] != c2[i]]
    if not pos:
        return 0.0, 0.0
    import itertools
    paths, sd, nd = 0, 0.0, 0.0
    for order in itertools.permutations(pos):
        cur = c1
        s = n = 0
        ok = True
        for i in order:
            nxt = cur[:i] + c2[i] + cur[i + 1:]
            a1, a2 = str(Seq(cur).translate()), str(Seq(nxt).translate())
            if a1 == "*" or a2 == "*":
                ok = False
                break
            if a1 == a2:
                s += 1
            else:
                n += 1
            cur = nxt
        if ok:
            paths += 1
            sd += s
            nd += n
    if paths == 0:
        return None, None
    return sd / paths, nd / paths


def ng86(codons1, codons2, min_codons=MIN_CODONS):
    """-> (dN, dS, n_codons) with the Jukes-Cantor correction, or Nones.

    min_codons is a parameter so the self-test can use a short sequence; the
    analysis always uses the plan's floor of 30.
    """
    S = N = Sd = Nd = 0.0
    used = 0
    for c1, c2 in zip(codons1, codons2):
        s1, s2 = syn_sites(c1), syn_sites(c2)
        if s1 is None or s2 is None:
            continue
        sd, nd = pairwise_diffs(c1, c2)
        if sd is None:
            continue
        s = (s1 + s2) / 2.0
        S += s
        N += 3.0 - s
        Sd += sd
        Nd += nd
        used += 1
    if used < min_codons or S <= 0 or N <= 0:
        return None, None, used

    def jc(p):
        if p is None or p >= 0.749999:
            return None
        return -0.75 * np.log(1.0 - 4.0 * p / 3.0)

    ds = jc(Sd / S)
    dn = jc(Nd / N)
    return dn, ds, used


def selftest():
    """NG86 must give zero on identical sequences, dS only for a silent change,
    and dN only for a replacement change."""
    cut = lambda s: [s[i:i+3] for i in range(0, len(s), 3)]
    a = "ATGGCTTCTAAAGGTGAAGAATTATTCACTGGTGTT"      # 12 codons
    dn, ds, n = ng86(cut(a), cut(a), min_codons=5)
    assert n == 12, n
    assert abs(dn) < 1e-9 and abs(ds) < 1e-9, (dn, ds)

    # GCT -> GCC is Ala to Ala: synonymous only
    syn = a[:5] + "C" + a[6:]
    assert str(Seq(a[3:6]).translate()) == str(Seq(syn[3:6]).translate())
    dn_s, ds_s, _ = ng86(cut(a), cut(syn), min_codons=5)
    assert abs(dn_s) < 1e-9 and ds_s > 0, (dn_s, ds_s)

    # GCT -> CCT is Ala to Pro: non-synonymous only
    non = a[:3] + "C" + a[4:]
    assert str(Seq(a[3:6]).translate()) != str(Seq(non[3:6]).translate())
    dn_n, ds_n, _ = ng86(cut(a), cut(non), min_codons=5)
    assert dn_n > 0 and abs(ds_n) < 1e-9, (dn_n, ds_n)

    # the floor is respected
    assert ng86(cut(a), cut(a))[0] is None, "the 30-codon floor is not applied"

    print(f"  NG86 self-test: identical dN {dn:.4f} dS {ds:.4f}; "
          f"silent dN {dn_s:.4f} dS {ds_s:.4f}; "
          f"replacement dN {dn_n:.4f} dS {ds_n:.4f}; floor applied")


# ----------------------------------------------------------------- inputs
def genus_of(org):
    if not isinstance(org, str) or not org.strip():
        return None
    w = re.sub(r"[\[\]]", "", org).split()
    if not w:
        return None
    if w[0] == "Candidatus" and len(w) > 1:
        return w[1]        # Candidatus is not a genus; the next word is
    return w[0]


def species_of(org):
    if not isinstance(org, str):
        return None
    w = re.sub(r"[\[\]]", "", org).split()
    if w and w[0] == "Candidatus":
        w = w[1:]
    return " ".join(w[:2]) if len(w) >= 2 else None


def read_fasta(p):
    if not p.exists():
        return {}
    op = gzip.open if str(p).endswith(".gz") else open
    with op(p, "rt") as fh:
        return {r.id: str(r.seq).upper() for r in SeqIO.parse(fh, "fasta")}


def genome_seqs(acc):
    """Contig sequences for a genome, from the Pseudofinder working copy."""
    return read_fasta(D / f"n2-pseudofinder/{acc}/pf_contigs.fasta")


def contig_key(contigs):
    """Map every caller's contig naming onto the pf_contigs names.

    The four callers rename contigs: Bakta uses contig_1.., DFAST sequence1..,
    and GenBank carries a version suffix (NC_004337.2) that Pseudofinder drops.
    Matching on LENGTH rather than on order, because order would have to be taken
    on trust and length can be checked. Lengths that are not unique fall back to
    order, and the fallback is counted so it cannot pass unnoticed.
    """
    by_len = collections.defaultdict(list)
    for name, s in contigs.items():
        by_len[len(s)].append(name)
    order = list(contigs)

    def resolve(name, length=None, index=None):
        if name in contigs:
            return name, False
        base = name.rsplit(".", 1)[0]            # NC_004337.2 -> NC_004337
        if base in contigs:
            return base, False
        if length is not None and len(by_len.get(length, [])) == 1:
            return by_len[length][0], False
        if index is not None and 0 <= index < len(order):
            return order[index], True            # ambiguous: resolved by order
        return None, False

    return resolve


def gff_region_lengths(path):
    """-> {contig name: length} from a GFF's ##sequence-region lines."""
    out = {}
    if not path.exists():
        return out
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("##sequence-region"):
            f = line.split()
            if len(f) >= 4:
                try:
                    out[f[1]] = int(f[3]) - int(f[2]) + 1
                except ValueError:
                    pass
        elif not line.startswith("#"):
            break
    return out


def sub(seq, start, end, strand):
    s = seq[start:end]
    return str(Seq(s).reverse_complement()) if strand == "-" else s


def gff_intervals(path, want_pseudo):
    """-> [(contig, start, end, strand)] from a GFF.

    want_pseudo None means take every feature, for a file that already holds only
    one kind.
    """
    out = []
    if not path.exists():
        return out
    for line in path.read_text(errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] not in ("CDS", "gene", "pseudogene"):
            continue
        # Coordinates can be fuzzy where a feature runs off a contig end, as in
        # ">21344". Strip the marker and keep the bound.
        try:
            start = int(f[3].lstrip("<>")) - 1
            end = int(f[4].lstrip("<>"))
        except ValueError:
            continue
        if want_pseudo is None:
            out.append((f[0], start, end, f[6]))
            continue
        if gff_asserts_pseudo(f) == want_pseudo:
            out.append((f[0], start, end, f[6]))
    return out


# A pseudogene ASSERTION in GFF3: an attribute whose KEY is pseudo or
# pseudogene, or the feature type itself.
_PSEUDO_ATTR = re.compile(r"(?:^|;)\s*pseudo(?:gene)?\s*=", re.I)


def gff_asserts_pseudo(f):
    """Does this GFF row assert the feature is a pseudogene?

    This was `"pseudo" in f[8].lower()`, a substring search over the whole
    attributes column -- which holds product, Name, gene, Dbxref and free-text
    note. It matched the NAMES of perfectly intact genes. Over the 523-genome
    corpus it produced 125,307 DFAST "pseudogenes", none carrying an assertion:
    the hits were `pseudomonas` (91,676 occurrences), `pseudotuberculosis`
    (15,843), `pseudouridine` (11,342), `pseudomallei` (10,463). Bakta's set
    was 39.1% contaminated the same way, nearly all tRNA pseudouridine
    synthases.

    run_dfast_all.sh and run_bakta_all.sh both warn against precisely this in
    their headers and take their COUNTS from each caller's own summary file.
    The warning never reached this interval path, so the count analyses were
    right and the interval analyses were wrong. See
    benchmark/SUBSTRING_CLASSIFIER_DEFECT.md.
    """
    return bool(_PSEUDO_ATTR.search(f[8])) or f[2] == "pseudogene"


# DFAST asserts its pseudogenes in pseudogene_summary.tsv, not in the GFF: the
# GFF carries no pseudogene feature type and no pseudo= attribute, and its
# disruption wording sits in free-text note values. The summary's
# classification column has four values, and `partial` is NOT a pseudogene
# call -- it marks low coverage against a reference, a length judgement.
# run_dfast_all.sh counts the disrupted classes only and records partial
# alongside, never added to them. This mirrors that decision exactly.
DFAST_DISRUPTED = frozenset({"frameshift", "internal_stop_codon",
                             "frameshift,internal_stop_codon"})
_DFAST_LOC = re.compile(r"^([^:]+):<?(\d+)\.\.>?(\d+)\(([+-])\)")


def dfast_pseudo(acc):
    """-> (locus ids, intervals) DFAST asserts as disrupted pseudogenes."""
    p = D / "n2-dfast" / acc / "pseudogene_summary.tsv"
    ids, ivs = set(), []
    if not p.exists():
        return ids, ivs
    for line in p.read_text(errors="replace").splitlines()[1:]:
        f = line.split("\t")
        if len(f) < 4 or f[3] not in DFAST_DISRUPTED:
            continue
        m = _DFAST_LOC.match(f[2])
        if not m:
            continue
        ids.add(f[0])
        ivs.append((m.group(1), int(m.group(2)) - 1, int(m.group(3)), m.group(4)))
    return ids, ivs


def dfast_intervals(acc, want_pseudo):
    """DFAST's pseudogenes from its summary, its intact genes from its GFF.

    Intact is every CDS the GFF reports whose locus id is not in the disrupted
    set, matched by ID rather than by coordinate so the two files cannot drift
    apart.
    """
    ids, ivs = dfast_pseudo(acc)
    if want_pseudo:
        return ivs
    out = []
    g = D / "n2-dfast" / acc / "genome.gff"
    if not g.exists():
        return out
    for line in g.read_text(errors="replace").splitlines():
        if line.startswith("##FASTA"):
            break
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] not in ("CDS", "gene"):
            continue
        fid = next((kv.split("=", 1)[1] for kv in f[8].split(";")
                    if kv.startswith("ID=")), None)
        if fid in ids:
            continue
        try:
            out.append((f[0], int(f[3].lstrip("<>")) - 1,
                        int(f[4].lstrip("<>")), f[6]))
        except ValueError:
            continue
    return out


def gbff_intervals(acc, want_pseudo):
    """PGAP's calls, from the RefSeq GenBank file.

    A RefSeq pseudogene is a `gene` feature carrying /pseudo or /pseudogene, with
    NO CDS feature beside it, because there is no protein to translate. Looking
    only at CDS features therefore finds none of them: a probe counted 571 grep
    hits for /pseudo in one genome and zero pseudogene CDS features. Pseudogenes
    come from gene features here, and intact genes from CDS features.
    """
    g = next((D / "n2-genomes/extracted").rglob(f"*{acc}*/genomic.gbff"), None)
    if g is None and (D / f"n2-genomes-expansion/{acc}").exists():
        g = next((D / f"n2-genomes-expansion/{acc}").rglob("genomic.gbff"), None)
    out = []
    if g is None:
        return out
    want_type = "gene" if want_pseudo else "CDS"
    for rec in SeqIO.parse(str(g), "genbank"):
        for ft in rec.features:
            if ft.type != want_type:
                continue
            is_pseudo = ("pseudo" in ft.qualifiers) or ("pseudogene" in ft.qualifiers)
            if is_pseudo != want_pseudo:
                continue
            out.append((rec.id, int(ft.location.start), int(ft.location.end),
                        "-" if ft.location.strand == -1 else "+"))
    return out


def calls_for(acc, caller, want_pseudo):
    """Genomic intervals a caller labelled pseudogene, or intact.

    Pseudofinder splits them into two files, so there the flag picks the file and
    every entry in it counts. The others mix both in one annotation, so there the
    flag filters within the file.
    """
    if caller == "PGAP":
        return gbff_intervals(acc, want_pseudo)
    if caller == "Pseudofinder":
        name = "pf_pseudos.gff" if want_pseudo else "pf_intact.gff"
        return gff_intervals(D / "n2-pseudofinder" / acc / name, want_pseudo=None)
    if caller == "Bakta":
        return gff_intervals(D / "n2-bakta" / acc / f"{acc}.gff3", want_pseudo)
    if caller == "DFAST":
        return dfast_intervals(acc, want_pseudo)
    return []


# ----------------------------------------------------------------- orthologues
ALIGNER = Align.PairwiseAligner(scoring="blastp", mode="global",
                                open_gap_score=-11, extend_gap_score=-1)


def best_hits(query_fa, db_faa, threads=8):
    """DIAMOND blastx: query nucleotide against partner proteins."""
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        dmnd = td / "db.dmnd"
        subprocess.run(["diamond", "makedb", "--in", str(db_faa), "-d", str(dmnd),
                        "--quiet", "--threads", str(threads)], check=True)
        out = td / "hits.tsv"
        subprocess.run(["diamond", "blastx", "-q", str(query_fa), "-d", str(dmnd),
                        "-o", str(out), "--quiet", "--threads", str(threads),
                        "--max-target-seqs", "1", "--outfmt", "6", "qseqid",
                        "sseqid", "pident", "length", "qstart", "qend", "slen",
                        "--evalue", "1e-5"], check=True)
        best = {}
        for line in out.read_text().splitlines():
            q, s, pid, ln, qs, qe, sl = line.split("\t")
            pid, ln, sl = float(pid), int(ln), int(sl)
            if pid < MIN_IDENT or ln / sl < MIN_COV:
                continue
            if q not in best:
                best[q] = (s, int(qs), int(qe))
        return best


def codon_pairs(cand_nt, partner_nt, qstart, qend):
    """Align candidate to partner in the PARTNER's frame; return codon pairs."""
    span = cand_nt[min(qstart, qend) - 1:max(qstart, qend)]
    if len(span) < MIN_CODONS * 3:
        return None
    if qstart > qend:
        span = str(Seq(span).reverse_complement())
    best = None
    for frame in (0, 1, 2):
        trimmed = span[frame:]
        trimmed = trimmed[:len(trimmed) - len(trimmed) % 3]
        if len(trimmed) < MIN_CODONS * 3:
            continue
        prot = str(Seq(trimmed).translate())
        ref = str(Seq(partner_nt[:len(partner_nt) - len(partner_nt) % 3]).translate())
        try:
            aln = ALIGNER.align(prot, ref)[0]
        except Exception:                                     # noqa: BLE001
            continue
        if best is None or aln.score > best[0]:
            best = (aln.score, aln, trimmed, partner_nt)
    if best is None:
        return None
    _, aln, qnt, snt = best
    qa, sa = aln[0], aln[1]
    c1, c2, qi, si = [], [], 0, 0
    for a, b in zip(qa, sa):
        if a != "-" and b != "-":
            x, y = qnt[qi*3:qi*3+3], snt[si*3:si*3+3]
            if len(x) == 3 and len(y) == 3:
                c1.append(x); c2.append(y)
        if a != "-":
            qi += 1
        if b != "-":
            si += 1
    return (c1, c2) if len(c1) >= MIN_CODONS else None


# ----------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", type=int, default=0)
    args = ap.parse_args()

    selftest()

    d = pd.read_csv(R / "corpus523.tsv", sep="\t")
    d["genus"] = d.organism.map(genus_of)
    d["species"] = d.organism.map(species_of)
    have = d[d.genus.notna() & d.species.notna()]
    pairs = []
    for g, grp in have.groupby("genus"):
        if grp.species.nunique() < 2:
            continue
        for row in grp.itertuples():
            cands = grp[grp.species != row.species].accession.sort_values()
            if len(cands):
                pairs.append((row.accession, cands.iloc[0], row.stratum, g))
    pairs.sort()
    print(f"\ntarget/partner pairs: {len(pairs)}")
    print(f"  genera: {len({g for *_, g in pairs})}")
    if args.probe:
        pairs = pairs[:args.probe]
        print(f"  PROBE: first {len(pairs)}")

    rows = []
    n_unmapped = collections.Counter()
    n_by_order = collections.Counter()
    for k, (acc, partner, stratum, gen) in enumerate(pairs, 1):
        contigs = genome_seqs(acc)
        pcds = read_fasta(D / f"n2-pseudofinder/{partner}/pf_intact.ffn")
        pfaa = D / f"n2-pseudofinder/{partner}/pf_intact.faa"
        if not contigs or not pcds or not pfaa.exists():
            continue
        resolve = contig_key(contigs)
        reglen = {
            "Bakta": gff_region_lengths(D / "n2-bakta" / acc / f"{acc}.gff3"),
            "DFAST": gff_region_lengths(D / "n2-dfast" / acc / "genome.gff"),
        }
        with tempfile.TemporaryDirectory() as td:
            qf = pathlib.Path(td) / "q.fna"
            meta = {}
            with qf.open("w") as fh:
                for caller in CALLERS:
                    ps = calls_for(acc, caller, True)
                    it = calls_for(acc, caller, False)
                    if not ps:
                        continue
                    if len(it) > len(ps):
                        idx = RNG.choice(len(it), len(ps), replace=False)
                        it = [it[i] for i in idx]
                    seen = {}
                    for tag, lst in (("pseudo", ps), ("intact", it)):
                        for i, (c, s, e, st) in enumerate(lst):
                            if c not in seen:
                                seen[c] = len(seen)
                            real, by_order = resolve(
                                c, reglen.get(caller, {}).get(c), seen[c])
                            if by_order:
                                n_by_order[caller] += 1
                            if real is None:
                                n_unmapped[caller] += 1
                                continue
                            if e - s < MIN_CODONS * 3 or e > len(contigs[real]):
                                continue
                            name = f"{caller}|{tag}|{i}"
                            seq = sub(contigs[real], s, e, st)
                            fh.write(f">{name}\n{seq}\n")
                            meta[name] = seq
            if not meta:
                continue
            hits = best_hits(qf, pfaa)
        for name, (sid, qs, qe) in hits.items():
            pnt = pcds.get(sid)
            if pnt is None:
                continue
            cp = codon_pairs(meta[name], pnt, qs, qe)
            if cp is None:
                continue
            dn, ds, used = ng86(*cp)
            if dn is None or ds is None or ds <= 0:
                continue
            caller, tag, _ = name.split("|")
            rows.append(dict(accession=acc, partner=partner, genus=gen,
                             stratum=stratum, caller=caller, call=tag,
                             dN=dn, dS=ds, omega=dn / ds, codons=used))
        if k % 10 == 0 or k == len(pairs):
            print(f"  {k}/{len(pairs)} genomes, {len(rows)} usable sequences")

    if n_unmapped or n_by_order:
        print("\ncontig-name reconciliation")
        for c in CALLERS:
            if n_unmapped[c] or n_by_order[c]:
                print(f"  {c}: {n_unmapped[c]} unmapped, "
                      f"{n_by_order[c]} resolved by order not length")

    t = pd.DataFrame(rows)
    if t.empty:
        print("no usable sequences; nothing to report")
        sys.exit(1)
    t.to_csv(R / "dnds_pairs.tsv", sep="\t", index=False)
    print(f"\nwrote results/dnds_pairs.tsv, {len(t)} rows")

    print(f"\n{'caller':14s} {'n pseudo':>9s} {'n intact':>9s} "
          f"{'med w pseudo':>13s} {'med w intact':>13s} {'gap':>7s} {'p':>10s}")
    summary, gaps = {}, {}
    for c in CALLERS:
        a = t[(t.caller == c) & (t.call == "pseudo")].omega
        b = t[(t.caller == c) & (t.call == "intact")].omega
        if len(a) < 5 or len(b) < 5:
            print(f"{c:14s} {len(a):9d} {len(b):9d}   too few")
            continue
        p = stats.mannwhitneyu(a, b, alternative="greater").pvalue
        gap = float(np.median(a) - np.median(b))
        gaps[c] = gap
        summary[c] = dict(n_pseudo=int(len(a)), n_intact=int(len(b)),
                          median_pseudo=float(np.median(a)),
                          median_intact=float(np.median(b)), gap=gap, p=float(p))
        print(f"{c:14s} {len(a):9d} {len(b):9d} {np.median(a):13.3f} "
              f"{np.median(b):13.3f} {gap:+7.3f} {p:10.3g}")

    verdicts = {}

    def verdict(k, text, ok, bad, val, force=None):
        """force: a label to record instead of judging, for a test that the
        design could not run as the plan wrote it. Without it a deviation gets
        stamped CONFIRMED on the strength of whatever was computed instead."""
        v = force or ("CONFIRMED" if ok else ("FAILED" if bad else "INCONCLUSIVE"))
        verdicts[k] = dict(prediction=text, value=val, verdict=v)
        print(f"  {k}  {v:13s} {text}\n        {val}")

    print("\n=== pre-committed predictions")
    higher = [c for c, s in summary.items() if s["gap"] > 0]
    verdict("G1", "dN/dS carries signal: pseudogene calls exceed intact controls",
            len(higher) == len(summary) and len(summary) >= 3,
            len(summary) - len(higher) >= 2,
            f"{len(higher)} of {len(summary)} callers have a positive gap: "
            + ", ".join(f"{c} {s['gap']:+.3f}" for c, s in summary.items()))

    if len(gaps) >= 2:
        spread = max(gaps.values()) - min(gaps.values())
        verdict("G2", "the callers differ in how well their calls match selection",
                spread >= 0.15, spread <= 0.05,
                f"best-worst gap spread {spread:.3f} "
                f"({max(gaps, key=gaps.get)} {max(gaps.values()):+.3f} vs "
                f"{min(gaps, key=gaps.get)} {min(gaps.values()):+.3f})")

    # G3 CANNOT BE RUN AS WRITTEN. The plan asks for the sequences one caller
    # calls pseudogene and another calls intact, which needs the four callers'
    # coordinates reconciled into one gene set; they do not share one. What is
    # computed here instead is every call against every control pooled, which
    # is close to a restatement of G1 and says nothing about the disagreement
    # cases. It is recorded as NOT_TESTED, not as a confirmation: an earlier
    # version of this file stamped it CONFIRMED on the strength of the
    # substitute, and the write-up had to contradict its own json.
    base = t[t.call == "intact"].omega
    allp = t[t.call == "pseudo"].omega
    if len(allp) > 5 and len(base) > 5:
        pd_ = stats.mannwhitneyu(allp, base, alternative="greater").pvalue
        verdict("G3", "where callers disagree, selection takes a side",
                False, False,
                f"NOT TESTED as written: no caller-neutral gene set exists in "
                f"these annotations. A pooled substitute gives all pseudogene "
                f"calls median {np.median(allp):.3f} against the all-intact "
                f"baseline {np.median(base):.3f}, p {pd_:.3g}, which does not "
                f"address the disagreement cases the plan specified.",
                force="NOT_TESTED")

    st_gap = {}
    for s in t.stratum.dropna().unique():
        a = t[(t.stratum == s) & (t.call == "pseudo")].omega
        b = t[(t.stratum == s) & (t.call == "intact")].omega
        if len(a) >= 5 and len(b) >= 5:
            st_gap[s] = float(np.median(a) - np.median(b))
    if "S3_recent_restriction" in st_gap and "S1_extreme_obligate" in st_gap:
        verdict("G4", "the signal is weakest where restriction is most recent",
                st_gap["S3_recent_restriction"] < st_gap["S1_extreme_obligate"],
                st_gap["S3_recent_restriction"] == max(st_gap.values()),
                "separation by stratum: "
                + ", ".join(f"{k.split('_')[0]} {v:+.3f}" for k, v in sorted(st_gap.items())))

    (R / "dnds.json").write_text(json.dumps(dict(
        plan_sha256=hashlib.sha256(PLAN.read_bytes()).hexdigest(),
        n_pairs=len(pairs), n_rows=int(len(t)), per_caller=summary,
        stratum_gap=st_gap, verdicts=verdicts), indent=1))
    nf = sum(1 for v in verdicts.values() if v["verdict"] == "FAILED")
    print(f"\nwrote results/dnds.json - {nf} of {len(verdicts)} FAILED")


if __name__ == "__main__":
    main()
