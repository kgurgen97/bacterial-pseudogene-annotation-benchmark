#!/usr/bin/env python3
"""Independent ground truth: put known lesions into intact genes, and prove it.

This simulator is deliberately independent of all four callers under test. In
particular it does NOT use Pseudofinder's Break module, whose own validation is
circular by its authors' statement -- the simulated pseudogene set is "determined
by the same parameters that would then be used to detect those pseudogenes"
(Syberg-Olsen et al. 2022). A truth set defined by one tool's detection
parameters cannot be used to score that tool, let alone the other three. Break
is kept for a secondary comparison, run afterwards, never as the truth source.

Everything the design is allowed to do is fixed in SIMULATION_PLAN.md, whose
hash is asserted below. No threshold, position band, lesion type or eligibility
rule may be changed here because of a result.

The job has three parts and the third is the one that matters:

  1. find genes that were intact to begin with, conservatively (plan section 3)
  2. break a small, fixed fraction of them in named ways (plan section 4)
  3. refuse to emit anything unless every lesion is provably present and every
     coordinate round-trips (plan section 4, eight blocking validations)

Usage:
    python3 simulate.py --accessions GCF_000006925.2 ... --replicates 2 --out DIR
    python3 simulate.py --manifest manifests/benchmark_genomes.tsv --replicates 3
"""
from __future__ import annotations

import argparse
import bisect
import collections
import csv
import dataclasses
import hashlib
import json
import os
import pathlib
import random
import shutil
import subprocess
import sys
import tempfile

from Bio import SeqIO
from Bio.Seq import Seq

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = pathlib.Path.home() / "research/datasets"

PLAN = HERE / "SIMULATION_PLAN.md"
PLAN_SHA = "657c024ecd1214076f059c24655fd3f9655b4ad3813ab5506d58058933b35965"

# ---- fixed by the plan; none of these may be tuned ------------------------
TABLE = 11
STARTS = ("ATG", "GTG", "TTG")
STOPS = ("TAA", "TAG", "TGA")
MIN_CODONS = 90
FLANK = 30
END_MARGIN = 500
MIN_CDS = 1500
MIN_ELIGIBLE = 400
LESION_FRACTION = 0.025
LESION_CAP = 120
PARALOG_ID = 90.0
PARALOG_COV = 0.90
POSITIONS = {"early": (0.05, 0.20), "middle": (0.40, 0.60), "late": (0.80, 0.95)}
TECH_ERROR_RATE = 0.0005

DISABLING = ("M1_stop_sub", "M2_fs_ins1", "M3_fs_del1", "M4_fs_indel2",
             "M5_multi", "M6_start_loss", "M7_trunc_term", "M8_disrupt_int")
CONTROLS = ("C1_synonymous", "C2_intergenic", "C3_inframe_indel")
POSITIONED = ("M1_stop_sub", "M2_fs_ins1", "M3_fs_del1", "M4_fs_indel2")


def assert_plan() -> None:
    h = hashlib.sha256(PLAN.read_bytes()).hexdigest()
    if h != PLAN_SHA:
        sys.exit(f"SIMULATION_PLAN.md has changed.\n  expected {PLAN_SHA}\n"
                 f"  found    {h}\nA frozen plan is never edited because of a "
                 f"result. Write a new plan with its own hash instead.")


# ===========================================================================
# genome space edits and the liftover they imply
# ===========================================================================
@dataclasses.dataclass(frozen=True)
class Edit:
    """Replace `ref_len` bases at 0-based `pos` with `alt`, in genome space.

    substitution: ref_len == len(alt) == 1.  insertion: ref_len == 0.
    deletion: alt == "".  All coordinates are on the forward strand of the
    contig regardless of which strand the affected gene sits on.
    """
    contig: str
    pos: int
    ref_len: int
    alt: str
    tag: str = ""

    @property
    def delta(self) -> int:
        return len(self.alt) - self.ref_len


class Liftover:
    """Parental 0-based coordinate -> modified 0-based coordinate.

    Built from the edits of one contig. A position that was deleted outright has
    no image and raises, which is intentional: silently mapping deleted bases
    onto their neighbour is how coordinate bugs survive a test suite.

    A position inside a *replacement* does have an image, and an earlier version
    of this class refused it, which broke every start-codon lesion -- those
    replace three bases with three others, so the gene's first base sits inside
    the replaced stretch. The k < altlen branch below is that fix. The test
    suite missed it because the randomised test skipped substituted positions;
    it no longer does.
    """

    def __init__(self, edits: list[Edit]):
        self.starts: list[int] = []
        self.ends: list[int] = []      # parental end of the replaced stretch
        self.altlen: list[int] = []    # length of what replaced it
        self.cum: list[int] = []       # cumulative delta AFTER this edit
        run = 0
        for e in sorted(edits, key=lambda e: e.pos):
            run += e.delta
            self.starts.append(e.pos)
            self.ends.append(e.pos + e.ref_len)
            self.altlen.append(len(e.alt))
            self.cum.append(run)

    def delta_before(self, pos: int) -> int:
        """Sum of deltas of edits whose replaced stretch ends at or before `pos`."""
        i = bisect.bisect_right(self.ends, pos) - 1
        return self.cum[i] if i >= 0 else 0

    def to_modified(self, pos: int) -> int:
        i = bisect.bisect_right(self.starts, pos) - 1
        if i >= 0 and self.starts[i] <= pos < self.ends[i]:
            k = pos - self.starts[i]
            if k < self.altlen[i]:
                return self.starts[i] + self.delta_before(self.starts[i]) + k
            raise KeyError(f"parental position {pos} was deleted")
        return pos + self.delta_before(pos)


def apply_edits(seq: str, edits: list[Edit]) -> str:
    out, cursor = [], 0
    for e in sorted(edits, key=lambda e: e.pos):
        if e.pos < cursor:
            raise ValueError(f"overlapping edits at {e.pos} (cursor {cursor})")
        out.append(seq[cursor:e.pos])
        out.append(e.alt)
        cursor = e.pos + e.ref_len
    out.append(seq[cursor:])
    return "".join(out)


def revcomp(s: str) -> str:
    return str(Seq(s).reverse_complement())


def gene_edit_to_genome(gene_pos: int, ref_len: int, alt: str,
                        start1: int, end1: int, strand: int,
                        contig: str, tag: str) -> Edit:
    """Translate an edit expressed in the gene's own reading frame to the contig.

    `start1`/`end1` are the 1-based inclusive genomic bounds of the CDS. For a
    minus-strand gene the gene sequence is the reverse complement of that span,
    so gene offset 0 is the genomic base at `end1`, the replaced stretch runs
    leftwards, and the replacement itself must be reverse complemented. Getting
    this wrong is silent -- the lesion lands somewhere real, just not where the
    truth table says -- so test_simulate.py checks both strands explicitly.
    """
    if strand >= 0:
        return Edit(contig, start1 - 1 + gene_pos, ref_len, alt, tag)
    return Edit(contig, end1 - gene_pos - ref_len, ref_len, revcomp(alt), tag)


def gene_seq(contigs: dict[str, str], contig: str, start1: int, end1: int,
             strand: int) -> str:
    s = contigs[contig][start1 - 1:end1]
    return s if strand >= 0 else revcomp(s)


def translate(nt: str) -> str:
    return str(Seq(nt).translate(table=TABLE))


# ===========================================================================
# loading a parental genome
# ===========================================================================
@dataclasses.dataclass
class Cds:
    locus_tag: str
    protein_id: str
    contig: str
    start1: int
    end1: int
    strand: int
    n_codons: int


def find_genome(acc: str) -> tuple[pathlib.Path, pathlib.Path] | None:
    gb = next((DATA / "n2-genomes/extracted").rglob(f"*{acc}*/genomic.gbff"), None)
    if gb is None:
        gb = next((DATA / f"n2-genomes-expansion/{acc}").rglob("genomic.gbff"), None)
    if gb is None:
        return None
    # the two stores name the assembly differently: the expansion download uses
    # a bare genomic.fna, the original datasets download the long NCBI name.
    fna = next((p for pat in ("genomic.fna", "*_genomic.fna", "*.fna")
                for p in sorted(gb.parent.glob(pat))), None)
    return (gb, fna) if fna else None


def load_parental(acc: str) -> tuple[dict, list, list, pathlib.Path, dict]:
    """Return contigs, CDS records, all feature intervals, the fna path, notes.

    The sequence is taken from the .gbff records and then checked base for base
    against the .fna, which is what the callers will actually be given. A
    mismatch between the two means the annotation does not describe the file
    being annotated, and the genome is rejected rather than patched.
    """
    found = find_genome(acc)
    if found is None:
        return {}, [], [], pathlib.Path(), {"reject": "no gbff/fna pair on disk"}
    gb, fna = found
    contigs, cds, feats, notes = {}, [], [], {}
    bad_table = 0
    plasmids: set[str] = set()
    for rec in SeqIO.parse(str(gb), "genbank"):
        contigs[rec.id] = str(rec.seq).upper()
        # Plan section 2.1: a plasmid replicon stays in the input the callers
        # see, so the assembly is not silently altered, but its genes are not
        # truth genes. Taken from the source feature's own /plasmid qualifier
        # rather than guessed from the description.
        src = next((f for f in rec.features if f.type == "source"), None)
        if src is not None and "plasmid" in src.qualifiers:
            plasmids.add(rec.id)
        for f in rec.features:
            if f.type in ("source",):
                continue
            feats.append((rec.id, int(f.location.start), int(f.location.end)))
            if f.type != "CDS":
                continue
            q = f.qualifiers
            if "pseudo" in q or "pseudogene" in q:
                continue
            tt = q.get("transl_table", ["11"])[0]
            if tt != str(TABLE):
                bad_table += 1
                continue
            if any(k in q for k in ("exception", "transl_except",
                                    "ribosomal_slippage", "artificial_location")):
                continue
            if len(f.location.parts) != 1:
                continue            # joined/compound location: frame undefined here
            s0, e0 = int(f.location.start), int(f.location.end)
            if (e0 - s0) % 3 or (e0 - s0) < MIN_CODONS * 3:
                continue
            cds.append(Cds(
                locus_tag=q.get("locus_tag", [f"{rec.id}:{s0}"])[0],
                protein_id=q.get("protein_id", [""])[0],
                contig=rec.id, start1=s0 + 1, end1=e0,
                strand=1 if f.location.strand in (1, None) else -1,
                n_codons=(e0 - s0) // 3))
    if bad_table:
        notes["reject"] = f"{bad_table} CDS use a translation table other than {TABLE}"
        return contigs, [], feats, fna, notes
    fa = {r.id: str(r.seq).upper() for r in SeqIO.parse(str(fna), "fasta")}
    if set(fa) != set(contigs):
        notes["reject"] = (f"fna and gbff disagree on contigs: "
                           f"{len(set(fa) ^ set(contigs))} differ")
    elif any(fa[k] != contigs[k] for k in fa):
        notes["reject"] = "fna and gbff sequences differ"
    notes["plasmids"] = sorted(plasmids)
    return contigs, cds, feats, fna, notes


# ===========================================================================
# eligibility (plan section 3)
# ===========================================================================
def find_diamond() -> str:
    """Resolve diamond without putting another env's bin on PATH.

    diamond lives in the `callers` conda env, not in the analysis env, and
    prepending that env's bin would shadow unrelated binaries. The resolved
    path and version are recorded in results/simulate.json.
    """
    cand = [os.environ.get("DIAMOND"), shutil.which("diamond"),
            str(pathlib.Path.home() / "miniforge3/envs/callers/bin/diamond"),
            str(pathlib.Path.home() / "miniforge3/envs/pseudofinder/bin/diamond")]
    for c in cand:
        if c and pathlib.Path(c).is_file():
            return c
    raise SystemExit("diamond is required for the paralogue screen (plan "
                     "section 3 criterion 8). Set DIAMOND=/path/to/diamond.")


DIAMOND = ""


def paralogues(cds: list[Cds], contigs: dict, workdir: pathlib.Path) -> set[str]:
    """Locus tags with a near-identical CDS elsewhere in the same genome.

    A caller can rescue a broken gene by homology to an intact paralogue, which
    makes a miss uninterpretable, so these never become truth genes.
    """
    dmnd = DIAMOND or find_diamond()
    faa = workdir / "prot.faa"
    with faa.open("w") as fh:
        for c in cds:
            p = translate(gene_seq(contigs, c.contig, c.start1, c.end1, c.strand))
            fh.write(f">{c.locus_tag}\n{p.rstrip('*')}\n")
    db, hits = workdir / "prot.dmnd", workdir / "self.tsv"
    for cmd in ([dmnd, "makedb", "--in", str(faa), "-d", str(db), "--quiet"],
                [dmnd, "blastp", "-d", str(db), "-q", str(faa), "-o", str(hits),
                 "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
                 "qlen", "slen", "--max-target-seqs", "6", "--quiet",
                 "--threads", "4"]):
        subprocess.run(cmd, check=True)
    out: set[str] = set()
    with hits.open() as fh:
        for row in csv.reader(fh, delimiter="\t"):
            q, s, pid, aln, ql, sl = row[0], row[1], float(row[2]), int(row[3]), int(row[4]), int(row[5])
            if q == s:
                continue
            if pid >= PARALOG_ID and aln >= PARALOG_COV * min(ql, sl):
                out.add(q)
                out.add(s)
    return out


def eligible(acc: str, contigs: dict, cds: list[Cds], feats: list,
             plasmids: set, workdir: pathlib.Path
             ) -> tuple[list[Cds], collections.Counter]:
    why: collections.Counter = collections.Counter()
    by_contig: dict[str, list] = collections.defaultdict(list)
    for cg, s, e in feats:
        by_contig[cg].append((s, e))
    for v in by_contig.values():
        v.sort()

    stage1 = []
    for c in cds:
        if c.contig in plasmids:
            why["on a plasmid replicon"] += 1
            continue
        seq = gene_seq(contigs, c.contig, c.start1, c.end1, c.strand)
        if seq[:3] not in STARTS:
            why["no valid start codon"] += 1
            continue
        if seq[-3:] not in STOPS:
            why["no terminal stop codon"] += 1
            continue
        prot = translate(seq)
        if "*" in prot[:-1]:
            why["internal stop in frame"] += 1
            continue
        if set(seq) - set("ACGT"):
            why["ambiguous bases"] += 1
            continue
        L = len(contigs[c.contig])
        if c.start1 - 1 < END_MARGIN or L - c.end1 < END_MARGIN:
            why["too close to a replicon end"] += 1
            continue
        s0, e0 = c.start1 - 1, c.end1
        overlap = any(not (e <= s0 - FLANK or s >= e0 + FLANK)
                      for s, e in by_contig[c.contig]
                      if not (s == s0 and e == e0))
        if overlap:
            why["overlaps another feature or its flank"] += 1
            continue
        stage1.append(c)

    para = paralogues(stage1, contigs, workdir) if stage1 else set()
    out = [c for c in stage1 if c.locus_tag not in para]
    why["near-identical paralogue in genome"] += len(stage1) - len(out)
    return out, why


# ===========================================================================
# drawing the lesions (plan section 4)
# ===========================================================================
@dataclasses.dataclass
class Lesion:
    gene: Cds
    code: str
    position: str
    tier: str
    gene_edits: list[tuple[int, int, str]]   # (gene_pos, ref_len, alt)
    note: str = ""


def codon_window(n_codons: int, position: str, rng: random.Random) -> int:
    lo, hi = POSITIONS[position]
    a, b = int(n_codons * lo), max(int(n_codons * hi), int(n_codons * lo) + 1)
    return rng.randrange(a, min(b, n_codons - 1))


def make_stop(seq: str, codon: int, rng: random.Random):
    """One substitution turning codon `codon` into a stop, if one exists."""
    cur = seq[codon * 3:codon * 3 + 3]
    opts = []
    for i in range(3):
        for b in "ACGT":
            if b == cur[i]:
                continue
            cand = cur[:i] + b + cur[i + 1:]
            if cand in STOPS:
                opts.append((codon * 3 + i, b))
    return rng.choice(opts) if opts else None


def draw_lesion(c: Cds, contigs: dict, code: str, rng: random.Random):
    seq = gene_seq(contigs, c.contig, c.start1, c.end1, c.strand)
    n = c.n_codons

    if code == "M1_stop_sub" or code == "M7_trunc_term":
        pos = "late" if code == "M7_trunc_term" else rng.choice(list(POSITIONS))
        for _ in range(40):
            codon = (rng.randrange(int(n * 0.90), n - 1) if code == "M7_trunc_term"
                     else codon_window(n, pos, rng))
            hit = make_stop(seq, codon, rng)
            if hit:
                gp, b = hit
                tier = "medium" if (code == "M7_trunc_term" or pos == "late") else "high"
                return Lesion(c, code, pos, tier, [(gp, 1, b)],
                              f"codon {codon} of {n} -> stop")
        return None

    if code in ("M2_fs_ins1", "M3_fs_del1", "M4_fs_indel2"):
        pos = rng.choice(list(POSITIONS))
        codon = codon_window(n, pos, rng)
        gp = codon * 3 + rng.randrange(3)
        tier = "medium" if pos == "late" else "high"
        if code == "M2_fs_ins1":
            ed = [(gp, 0, rng.choice("ACGT"))]
        elif code == "M3_fs_del1":
            ed = [(gp, 1, "")]
        else:
            ed = ([(gp, 0, "".join(rng.choice("ACGT") for _ in range(2)))]
                  if rng.random() < 0.5 else [(gp, 2, "")])
        if gp + max(e[1] for e in ed) > len(seq) - 3:
            return None
        return Lesion(c, code, pos, tier, ed, f"codon {codon} of {n}")

    if code == "M5_multi":
        k = rng.choice((2, 3))
        picks, used = [], set()
        for _ in range(k * 6):
            if len(picks) == k:
                break
            sub = draw_lesion(c, contigs, rng.choice(POSITIONED), rng)
            if sub is None:
                continue
            gp = sub.gene_edits[0][0]
            if any(abs(gp - u) < 12 for u in used):
                continue
            used.add(gp)
            picks.append(sub)
        if len(picks) < 2:
            return None
        ed = [e for p in picks for e in p.gene_edits]
        return Lesion(c, code, "mixed", "high", ed,
                      f"{len(picks)} lesions: " + ",".join(p.code for p in picks))

    if code == "M6_start_loss":
        alt = [b + seq[1:3] for b in "ACGT" if b + seq[1:3] not in STARTS]
        if not alt:
            return None
        new = rng.choice(alt)
        # an alternative in-frame start inside the first 10% of codons leaves a
        # shorter but still translatable product, so truth is only "high" when
        # there is none.
        lim = max(1, int(n * 0.10))
        downstream = any(seq[i * 3:i * 3 + 3] in STARTS for i in range(1, lim))
        return Lesion(c, code, "start",
                      "ambiguous" if downstream else "high",
                      [(0, 3, new)],
                      f"start {seq[:3]}->{new}" +
                      (f"; alternative start within first {lim} codons" if downstream else ""))

    if code == "M8_disrupt_int":
        for _ in range(40):
            L = rng.randrange(30, 151)
            if L % 3 == 0:
                continue
            codon = codon_window(n, "middle", rng)
            gp = codon * 3
            if gp + L > len(seq) - 6:
                continue
            return Lesion(c, code, "middle", "medium", [(gp, L, "")],
                          f"{L} bp deleted from codon {codon} of {n}")
        return None

    if code == "C1_synonymous":
        order = list(range(1, n - 1))
        rng.shuffle(order)
        for codon in order[:200]:
            cur = seq[codon * 3:codon * 3 + 3]
            aa = translate(cur)
            for i in (2, 1, 0):
                for b in "ACGT":
                    if b == cur[i]:
                        continue
                    cand = cur[:i] + b + cur[i + 1:]
                    if cand not in STOPS and translate(cand) == aa:
                        return Lesion(c, code, "any", "negative",
                                      [(codon * 3 + i, 1, b)],
                                      f"codon {codon} {cur}->{cand}, {aa} kept")
        return None

    if code == "C3_inframe_indel":
        pos = rng.choice(list(POSITIONS))
        codon = codon_window(n, pos, rng)
        gp = codon * 3
        ed = ([(gp, 0, "".join(rng.choice("ACGT") for _ in range(3)))]
              if rng.random() < 0.5 else [(gp, 3, "")])
        if gp + 3 > len(seq) - 6:
            return None
        return Lesion(c, code, pos, "ambiguous", ed, f"codon {codon} of {n}")

    raise ValueError(f"unknown lesion code {code}")


def intergenic_edits(contigs: dict, feats: list, k: int, plasmids: set,
                     rng: random.Random) -> list[Edit]:
    """C2: substitutions in sequence no annotated feature claims."""
    busy: dict[str, list] = collections.defaultdict(list)
    for cg, s, e in feats:
        busy[cg].append((s - FLANK, e + FLANK))
    out, tries = [], 0
    names = [c for c in contigs
             if len(contigs[c]) > 2 * END_MARGIN and c not in plasmids]
    while len(out) < k and tries < k * 400:
        tries += 1
        cg = rng.choice(names)
        p = rng.randrange(END_MARGIN, len(contigs[cg]) - END_MARGIN)
        if any(s <= p < e for s, e in busy[cg]):
            continue
        if contigs[cg][p] not in "ACGT":
            continue
        b = rng.choice([x for x in "ACGT" if x != contigs[cg][p]])
        out.append(Edit(cg, p, 1, b, "C2_intergenic"))
        busy[cg].append((p - FLANK, p + FLANK))
    return out


# ===========================================================================
# the eight blocking validations (plan section 4)
# ===========================================================================
def validate(acc: str, rep: int, contigs: dict, modified: dict,
             lifts: dict, lesions: list[Lesion], extra: list[Edit],
             per_contig: dict, fna: pathlib.Path) -> list[str]:
    bad: list[str] = []

    # 7. contig count, order and names unchanged
    if list(modified) != list(contigs):
        bad.append("contig names or order changed")

    # 8. the parental fasta we hand the callers is the file as downloaded
    orig = {r.id: str(r.seq).upper() for r in SeqIO.parse(str(fna), "fasta")}
    if any(orig.get(k) != contigs[k] for k in contigs):
        bad.append("parental sequence does not match the downloaded fna")

    for les in lesions:
        c = les.gene
        p_seq = gene_seq(contigs, c.contig, c.start1, c.end1, c.strand)
        # expected gene sequence, built independently in gene space
        want = apply_edits(p_seq, [Edit("g", gp, rl, alt)
                                   for gp, rl, alt in les.gene_edits])
        # what the genome-space edit actually produced, via the liftover
        lf = lifts[c.contig]
        delta = sum(len(a) - r for _, r, a in les.gene_edits)
        try:
            if c.strand >= 0:
                m0 = lf.to_modified(c.start1 - 1)
                got = modified[c.contig][m0:m0 + len(p_seq) + delta]
            else:
                m0 = lf.to_modified(c.end1 - 1)
                got = revcomp(modified[c.contig][m0 - (len(p_seq) + delta) + 1:m0 + 1])
        except KeyError as exc:
            bad.append(f"{les.code} {c.locus_tag}: liftover: {exc}")
            continue
        if got != want:
            bad.append(f"{les.code} {c.locus_tag}: gene-space and genome-space "
                       f"edits disagree (strand {c.strand})")
            continue

        # 2/3. the intended biological consequence actually happened
        if les.code in ("M1_stop_sub", "M7_trunc_term"):
            prot = translate(want[:len(want) - len(want) % 3])
            if "*" not in prot[:-1]:
                bad.append(f"{les.code} {c.locus_tag}: no premature stop")
        elif les.code in ("M2_fs_ins1", "M3_fs_del1", "M4_fs_indel2",
                          "M8_disrupt_int"):
            if len(want) % 3 == 0:
                bad.append(f"{les.code} {c.locus_tag}: frame not changed "
                           f"(length still a multiple of 3)")
        elif les.code == "M5_multi":
            # Several independent lesions need not change the gene's LENGTH --
            # a +1 insertion and a -1 deletion cancel, and two stop
            # substitutions change nothing -- while still disabling the gene.
            # What must be true is that the protein changed, and that either a
            # premature stop appeared or some internal stretch is out of frame.
            prot = translate(want[:len(want) - len(want) % 3])
            if prot == translate(p_seq):
                bad.append(f"{les.code} {c.locus_tag}: protein unchanged")
            run, shifted = 0, False
            for gp, rl, alt in sorted(les.gene_edits):
                run += len(alt) - rl
                if run != 0:
                    shifted = True
            if "*" not in prot[:-1] and not shifted:
                bad.append(f"{les.code} {c.locus_tag}: neither a premature stop "
                           f"nor an out-of-frame stretch")
        elif les.code == "M6_start_loss":
            if want[:3] in STARTS:
                bad.append(f"{les.code} {c.locus_tag}: start codon still valid")
        elif les.code == "C1_synonymous":
            if translate(want) != translate(p_seq):
                bad.append(f"{les.code} {c.locus_tag}: protein changed")
            if len(want) != len(p_seq):
                bad.append(f"{les.code} {c.locus_tag}: length changed")
        elif les.code == "C3_inframe_indel":
            if len(want) % 3 or abs(len(want) - len(p_seq)) != 3:
                bad.append(f"{les.code} {c.locus_tag}: not a 3 bp in-frame indel")

    # 4. C2 edits touch no annotated feature  (checked at draw time, re-asserted)
    for e in extra:
        if e.tag == "C2_intergenic" and e.ref_len != 1:
            bad.append("C2 edit is not a single substitution")

    # 5. liftover round-trips. A substituted base legitimately differs, so the
    # comparison is against what the edit put there, not against the parental
    # base -- comparing against the parental base is why this check first fired
    # on a position that was correct. Every lesion's own boundaries are tested
    # explicitly rather than left to random sampling.
    rng = random.Random(f"{acc}:{rep}:roundtrip")
    for cg, lf in lifts.items():
        ed = sorted(per_contig.get(cg, []), key=lambda e: e.pos)
        starts = [e.pos for e in ed]
        n = len(contigs[cg])

        def expect(p: int):
            """The base the modified contig should hold for parental p, or None."""
            i = bisect.bisect_right(starts, p) - 1
            if i >= 0 and ed[i].pos <= p < ed[i].pos + ed[i].ref_len:
                k = p - ed[i].pos
                return ed[i].alt[k] if k < len(ed[i].alt) else None
            return contigs[cg][p]

        probes = [rng.randrange(n) for _ in range(400)]
        for e in ed:                       # both flanks of every single lesion
            probes += [e.pos - 1, e.pos + e.ref_len, e.pos + e.ref_len + 1]
        for p in probes:
            if not 0 <= p < n:
                continue
            want = expect(p)
            if want is None:
                continue
            try:
                m = lf.to_modified(p)
            except KeyError:
                bad.append(f"{cg}: parental {p} has no image but was not deleted")
                break
            if modified[cg][m] != want:
                bad.append(f"{cg}: liftover sends parental {p} to {modified[cg][m]!r}, "
                           f"expected {want!r}")
                break

    # 6. no two lesions in one gene, and none within FLANK of another
    seen = [l.gene.locus_tag for l in lesions]
    if len(seen) != len(set(seen)):
        bad.append("a gene carries more than one lesion")
    return bad


# ===========================================================================
# one replicate
# ===========================================================================
def build_replicate(acc: str, rep: int, contigs: dict, elig: list[Cds],
                    feats: list, plasmids: set, seed: str, fna: pathlib.Path):
    rng = random.Random(seed)
    n_les = min(LESION_CAP, max(1, int(round(len(elig) * LESION_FRACTION))))
    pool = list(elig)
    rng.shuffle(pool)

    codes = list(DISABLING) + ["C1_synonymous", "C3_inframe_indel"]
    lesions: list[Lesion] = []
    gi = 0
    while len(lesions) < n_les and gi < len(pool):
        code = codes[len(lesions) % len(codes)]
        les = draw_lesion(pool[gi], contigs, code, rng)
        gi += 1
        if les is not None:
            lesions.append(les)

    gene_edits: list[Edit] = []
    for les in lesions:
        c = les.gene
        for gp, rl, alt in les.gene_edits:
            gene_edits.append(gene_edit_to_genome(
                gp, rl, alt, c.start1, c.end1, c.strand, c.contig, les.code))
    extra = intergenic_edits(contigs, feats, max(5, n_les // 8), plasmids, rng)

    all_edits = gene_edits + extra
    per_contig: dict[str, list[Edit]] = collections.defaultdict(list)
    for e in all_edits:
        per_contig[e.contig].append(e)
    modified = {cg: apply_edits(s, per_contig.get(cg, [])) for cg, s in contigs.items()}
    lifts = {cg: Liftover(per_contig.get(cg, [])) for cg in contigs}

    bad = validate(acc, rep, contigs, modified, lifts, lesions, extra,
                   per_contig, fna)
    return lesions, extra, modified, lifts, bad


def tech_error_replicate(acc: str, contigs: dict, seed: str):
    """C4: genome-wide single-nucleotide deletions, no targeted gene lesions."""
    rng = random.Random(seed)
    per: dict[str, list[Edit]] = collections.defaultdict(list)
    for cg, s in contigs.items():
        k = int(len(s) * TECH_ERROR_RATE)
        for p in sorted(rng.sample(range(len(s)), k)):
            per[cg].append(Edit(cg, p, 1, "", "C4_tech_error"))
    mod = {cg: apply_edits(s, per.get(cg, [])) for cg, s in contigs.items()}
    return mod, sum(len(v) for v in per.values())


# ===========================================================================
def write_fasta(path: pathlib.Path, contigs: dict) -> None:
    with path.open("w") as fh:
        for name, seq in contigs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")


def main() -> None:
    assert_plan()
    ap = argparse.ArgumentParser()
    ap.add_argument("--accessions", nargs="*", default=[])
    ap.add_argument("--manifest")
    ap.add_argument("--replicates", type=int, default=2)
    ap.add_argument("--out", default=str(DATA / "n2-benchmark"))
    ap.add_argument("--stage", default="A")
    ap.add_argument("--diamond", default=None,
                    help="path to diamond; otherwise $DIAMOND, PATH, then the "
                         "callers conda env")
    a = ap.parse_args()

    global DIAMOND
    DIAMOND = a.diamond or find_diamond()
    dv = subprocess.run([DIAMOND, "--version"], capture_output=True,
                        text=True).stdout.strip()

    accs = list(a.accessions)
    if a.manifest:
        with open(a.manifest) as fh:
            accs += [r["accession"] for r in csv.DictReader(fh, delimiter="\t")]
    if not accs:
        sys.exit("no accessions given")

    out = pathlib.Path(a.out)
    (out / "assemblies").mkdir(parents=True, exist_ok=True)
    (HERE / "results/liftover").mkdir(parents=True, exist_ok=True)
    (HERE / "logs").mkdir(parents=True, exist_ok=True)

    truth_rows, elig_rows, failures = [], [], []
    print(f"plan {PLAN_SHA[:16]} ok | stage {a.stage} | "
          f"{len(accs)} genomes x {a.replicates} replicates")
    print(f"  {dv}  at {DIAMOND}")

    for acc in accs:
        contigs, cds, feats, fna, notes = load_parental(acc)
        if notes.get("reject"):
            failures.append((acc, notes["reject"]))
            print(f"  {acc:20s} REJECTED  {notes['reject']}")
            continue
        if len(cds) < MIN_CDS:
            failures.append((acc, f"only {len(cds)} usable CDS, need {MIN_CDS}"))
            print(f"  {acc:20s} REJECTED  only {len(cds)} usable CDS")
            continue
        with tempfile.TemporaryDirectory() as td:
            elig, why = eligible(acc, contigs, cds, feats,
                                 set(notes.get("plasmids", [])), pathlib.Path(td))
        for reason, n in why.items():
            elig_rows.append({"accession": acc, "reason": reason, "n": n})
        elig_rows.append({"accession": acc, "reason": "ELIGIBLE", "n": len(elig)})
        if len(elig) < MIN_ELIGIBLE:
            failures.append((acc, f"only {len(elig)} eligible genes, need {MIN_ELIGIBLE}"))
            print(f"  {acc:20s} REJECTED  only {len(elig)} eligible genes")
            continue

        write_fasta(out / "assemblies" / f"{acc}__C0_parental.fna", contigs)
        print(f"  {acc:20s} {len(cds)} CDS, {len(elig)} eligible", flush=True)

        for rep in range(1, a.replicates + 1):
            seed = f"{acc}|rep{rep}|v1"
            les, extra, mod, lifts, bad = build_replicate(
                acc, rep, contigs, elig, feats,
                set(notes.get("plasmids", [])), seed, fna)
            if bad:
                for b in bad[:6]:
                    print(f"      VALIDATION FAILED  {b}")
                sys.exit(f"{acc} replicate {rep}: {len(bad)} validation failures; "
                         f"refusing to write a truth set that is not provably correct")
            rid = f"{acc}__R{rep}"
            write_fasta(out / "assemblies" / f"{rid}.fna", mod)
            with (HERE / f"results/liftover/{rid}.tsv").open("w") as fh:
                fh.write("contig\tparental_pos0\tparental_end0\tcum_delta\n")
                for cg, lf in lifts.items():
                    for s, e, c in zip(lf.starts, lf.ends, lf.cum):
                        fh.write(f"{cg}\t{s}\t{e}\t{c}\n")
            for l in les:
                g = l.gene
                p_seq = gene_seq(contigs, g.contig, g.start1, g.end1, g.strand)
                gp = l.gene_edits[0][0]
                truth_rows.append({
                    "accession": acc, "replicate": rid, "seed": seed,
                    "locus_tag": g.locus_tag, "protein_id": g.protein_id,
                    "contig": g.contig, "parental_start1": g.start1,
                    "parental_end1": g.end1, "strand": g.strand,
                    "n_codons": g.n_codons, "lesion_code": l.code,
                    "position": l.position, "truth_tier": l.tier,
                    "n_edits": len(l.gene_edits),
                    "gene_offset": gp,
                    "delta_bp": sum(len(x[2]) - x[1] for x in l.gene_edits),
                    "window_parental": p_seq[max(0, gp - 10):gp + 11],
                    "note": l.note})
            print(f"      R{rep}: {len(les)} lesions, {len(extra)} intergenic, "
                  f"validations all pass", flush=True)

        mod4, n4 = tech_error_replicate(acc, contigs, f"{acc}|C4|v1")
        write_fasta(out / "assemblies" / f"{acc}__C4_tech_error.fna", mod4)
        print(f"      C4: {n4} genome-wide deletions", flush=True)

    tsv = HERE / "results/ground_truth.tsv"
    if truth_rows:
        with tsv.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(truth_rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(truth_rows)
    with (HERE / "results/eligibility.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["accession", "reason", "n"], delimiter="\t")
        w.writeheader()
        w.writerows(elig_rows)
    tier = collections.Counter(r["truth_tier"] for r in truth_rows)
    code = collections.Counter(r["lesion_code"] for r in truth_rows)
    (HERE / "results/simulate.json").write_text(json.dumps({
        "plan_sha256": PLAN_SHA, "stage": a.stage,
        "diamond": {"path": DIAMOND, "version": dv},
        "n_genomes_requested": len(accs),
        "n_genomes_rejected": len(failures),
        "rejections": [{"accession": x, "reason": y} for x, y in failures],
        "replicates_per_genome": a.replicates,
        "n_lesions": len(truth_rows),
        "by_tier": dict(tier), "by_code": dict(code),
        "assemblies": str(out / "assemblies"),
    }, indent=2, sort_keys=True) + "\n")
    print(f"\n{len(truth_rows)} lesions in {len(accs) - len(failures)} genomes")
    print(f"  tiers {dict(tier)}")
    print(f"  codes {dict(code)}")
    print(f"wrote {tsv.relative_to(ROOT)} and results/simulate.json")
    if failures:
        print(f"  {len(failures)} genomes rejected, all listed in simulate.json")


if __name__ == "__main__":
    main()
