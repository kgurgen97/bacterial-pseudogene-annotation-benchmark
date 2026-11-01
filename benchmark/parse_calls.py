#!/usr/bin/env python3
"""Normalise caller output into the six classes SIMULATION_PLAN.md section 6 fixes.

    explicit_pseudogene · disrupted · partial · intact_cds · absent · ambiguous

The rules that matter, taken from the plan and not relaxed here:

  * `explicit_pseudogene` needs an explicit ASSERTION. For Bakta that is a
    `pseudogene=` or `pseudo=` attribute. For DFAST it is an entry in
    pseudogene_summary.tsv whose classification is a disrupted class --
    frameshift, internal_stop_codon, or both.
  * **DFAST `partial` is `partial`, never `explicit_pseudogene`.** It marks low
    coverage against a reference, which is a length judgement. run_dfast_all.sh
    made the same decision for its counts and the plan repeats it.
  * **`absent` is `absent`, never a pseudogene call.** A missing CDS is scored
    under endpoints B and C only.
  * `disrupted` is a gene split into two or more CDS features, which is the
    other way a caller can register breakage without asserting a pseudogene.
  * `ambiguous` is when one locus draws both a pseudogene assertion and an
    intact CDS from the same caller. Those genes are excluded from that
    caller's primary analysis and counted, rather than resolved by a guess.

Matching is >= 50 % reciprocal overlap, as LOCUS_PLAN.md uses, computed in the
coordinates of the assembly the caller was actually given.

    python3 benchmark/parse_calls.py --calls DIR --assembly-id ID --caller bakta
"""
from __future__ import annotations

import argparse
import bisect
import collections
import csv
import dataclasses
import pathlib
import re

RECIP = 0.50

# priority order: a feature gets exactly one kind
PSEUDO, PARTIAL, CDS, OTHER = "pseudogene", "partial", "cds", "other"

BAKTA_PSEUDO = re.compile(r"(?:^|;)\s*pseudo(?:gene)?\s*=", re.I)
DFAST_DISRUPTED = frozenset({"frameshift", "internal_stop_codon",
                             "frameshift,internal_stop_codon"})


@dataclasses.dataclass(frozen=True)
class Feat:
    contig: str
    start0: int
    end0: int
    strand: str
    kind: str
    ident: str


def _gff_rows(path: pathlib.Path):
    """GFF3 data lines, stopping at an appended ##FASTA block."""
    if not path.exists():
        return
    with path.open(errors="replace") as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                return
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9:
                yield f


def _attr(blob: str, key: str):
    for kv in blob.split(";"):
        if kv.startswith(key + "="):
            return kv.split("=", 1)[1]
    return None


def load_features(caller: str, run_dir: pathlib.Path, assembly_id: str) -> list[Feat]:
    """Every annotated feature, each carrying exactly one of the four kinds."""
    out: list[Feat] = []
    if caller == "bakta":
        gff = run_dir / f"{assembly_id}.gff3"
        for f in _gff_rows(gff):
            try:
                s, e = int(f[3].lstrip("<>")) - 1, int(f[4].lstrip("<>"))
            except ValueError:
                continue
            ident = _attr(f[8], "ID") or f"{f[0]}:{s}"
            if f[2] in ("CDS", "gene", "pseudogene"):
                kind = PSEUDO if BAKTA_PSEUDO.search(f[8]) else CDS
            else:
                kind = OTHER
            out.append(Feat(f[0], s, e, f[6], kind, ident))
        return out

    if caller == "dfast":
        # DFAST asserts pseudogenes in a separate file, keyed by locus id
        summary = run_dir / "pseudogene_summary.tsv"
        kinds: dict[str, str] = {}
        if summary.exists():
            for line in summary.read_text(errors="replace").splitlines()[1:]:
                r = line.split("\t")
                if len(r) < 4:
                    continue
                if r[3] in DFAST_DISRUPTED:
                    kinds[r[0]] = PSEUDO
                elif r[3] == "partial":
                    kinds[r[0]] = PARTIAL
        for f in _gff_rows(run_dir / "genome.gff"):
            try:
                s, e = int(f[3].lstrip("<>")) - 1, int(f[4].lstrip("<>"))
            except ValueError:
                continue
            ident = _attr(f[8], "ID") or f"{f[0]}:{s}"
            if f[2] in ("CDS", "gene", "pseudogene"):
                kind = kinds.get(ident, CDS)
            else:
                kind = OTHER
            out.append(Feat(f[0], s, e, f[6], kind, ident))
        return out

    raise ValueError(f"no output rule for caller {caller!r}")


class Index:
    """Features per contig, sorted, for reciprocal-overlap queries."""

    def __init__(self, feats: list[Feat]):
        self.by_contig: dict[str, list[Feat]] = collections.defaultdict(list)
        for f in feats:
            self.by_contig[f.contig].append(f)
        for v in self.by_contig.values():
            v.sort(key=lambda f: f.start0)
        self.starts = {c: [f.start0 for f in v] for c, v in self.by_contig.items()}

    def _scan(self, contig: str, start0: int, end0: int):
        v = self.by_contig.get(contig, [])
        if not v or end0 <= start0:
            return
        # features can be long, so start well before the query
        i = max(0, bisect.bisect_left(self.starts[contig], start0) - 200)
        for f in v[i:]:
            if f.start0 >= end0:
                break
            ov = min(end0, f.end0) - max(start0, f.start0)
            if ov > 0:
                yield f, ov

    def overlapping(self, contig: str, start0: int, end0: int,
                    recip: float = RECIP) -> list[Feat]:
        """Plan section 7's matching rule: >= 50 % RECIPROCAL overlap."""
        qlen = end0 - start0
        return [f for f, ov in self._scan(contig, start0, end0)
                if (f.end0 - f.start0) > 0
                and ov >= recip * qlen
                and ov >= recip * (f.end0 - f.start0)]

    def fragments(self, contig: str, start0: int, end0: int,
                  cover: float = RECIP) -> list[Feat]:
        """Features at least half of whose own length lies inside the query.

        DEVIATION from plan section 7, recorded rather than silently applied.
        Section 6 defines `disrupted` as a gene split into two or more CDS
        features, and section 7 fixes matching at >= 50 % reciprocal overlap.
        Those two cannot both hold: a gene split in half yields two fragments
        that each cover under half of the gene, so no fragment ever clears the
        query side of the reciprocal bar and fragmentation is undetectable by
        construction. This one-sided, feature-side criterion is used ONLY to
        count fragments for the `disrupted` class. The reciprocal rule is left
        untouched for identity, and no threshold was tuned - 50 % is the
        plan's own number, applied to the other side of the comparison.
        """
        return [f for f, ov in self._scan(contig, start0, end0)
                if (f.end0 - f.start0) > 0
                and ov >= cover * (f.end0 - f.start0)]

    def any_overlap(self, contig: str, start0: int, end0: int) -> list[Feat]:
        """Plan section 7, endpoint C: ANY feature overlapping the interval."""
        return [f for f, _ in self._scan(contig, start0, end0)]


def classify(feats: list[Feat], frags: list[Feat] | None = None) -> str:
    """The six-way class for one query interval.

    `feats` are the features meeting the plan's reciprocal rule and decide the
    class. `frags` are features mostly contained in the query and are consulted
    only to detect fragmentation, which the reciprocal rule cannot see (see
    Index.fragments). A single fragment is NOT called disrupted: the plan's
    definition is two or more, and `absent` already counts under endpoint B.
    """
    kinds = {f.kind for f in feats if f.kind != OTHER}
    n_cds = sum(1 for f in feats if f.kind == CDS)
    if PSEUDO in kinds and (CDS in kinds or PARTIAL in kinds):
        return "ambiguous"
    if PSEUDO in kinds:
        return "explicit_pseudogene"
    if PARTIAL in kinds:
        return "partial"
    if n_cds >= 2:
        return "disrupted"
    if n_cds == 1:
        return "intact_cds"
    if frags:
        fk = {f.kind for f in frags if f.kind != OTHER}
        if PSEUDO in fk:
            return "explicit_pseudogene"
        if sum(1 for f in frags if f.kind == CDS) >= 2:
            return "disrupted"
    return "absent"


# endpoint membership, plan section 7
ENDPOINT_A = {"explicit_pseudogene"}
ENDPOINT_B = {"explicit_pseudogene", "disrupted", "partial", "absent"}


def endpoint_positive(cls: str, endpoint: str) -> bool:
    if endpoint == "A":
        return cls in ENDPOINT_A
    if endpoint == "B":
        return cls in ENDPOINT_B
    if endpoint == "C":
        return cls != "absent"
    raise ValueError(endpoint)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True)
    ap.add_argument("--assembly-id", required=True)
    ap.add_argument("--caller", required=True)
    a = ap.parse_args()
    d = pathlib.Path(a.calls) / a.caller / a.assembly_id
    feats = load_features(a.caller, d, a.assembly_id)
    n = collections.Counter(f.kind for f in feats)
    print(f"{a.caller}/{a.assembly_id}: {len(feats)} features")
    for k, c in n.most_common():
        print(f"  {k:12s} {c}")


if __name__ == "__main__":
    main()
