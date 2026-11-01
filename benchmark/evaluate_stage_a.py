#!/usr/bin/env python3
"""Stage A evaluation: what did the callers see, against independent truth?

Plan section 7. Three endpoints, kept strictly apart:

  A  explicit pseudogene classification
  B  disruption or abnormal structure detected (pseudogene, disrupted,
     partial, or absent)
  C  the gene recovered as any annotated feature at all

Every endpoint is PAIRED against the parental assembly: a gene counts as a true
positive only if the caller calls it positive in the modified replicate and did
NOT call it positive in the parental. A false positive is a gene called
positive in the modified replicate that carried no lesion and was not already
positive in the parental. That removes each caller's baseline call rate, which
is the only way a conservative caller and a liberal one can be compared.

Scoring is per TRUTH GENE, so one broken gene split into several annotations
cannot count as several true positives. The evaluation universe is the eligible
gene set, not the set of predictions.

Only Bakta and DFAST are evaluated. PGAP cannot run here and
Pseudofinder-on-PGAP needs PGAP's annotation of the modified genome, so both
are BLOCKED and reported as blocked rather than approximated.

    python3 benchmark/evaluate_stage_a.py
"""
from __future__ import annotations

import bisect
import collections
import csv
import json
import pathlib
import random

HERE = pathlib.Path(__file__).resolve().parent
ASM = pathlib.Path.home() / "research/datasets/n2-benchmark/assemblies"
CALLS = pathlib.Path.home() / "research/datasets/n2-benchmark/calls"

import parse_calls as PC  # noqa: E402
import simulate as S  # noqa: E402

CALLERS = ["bakta", "dfast"]
BLOCKED = {"PGAP": "no container runtime, and ~30 GB of data against 17 GB free",
           "Pseudofinder-on-PGAP": "requires PGAP's annotation of the modified "
                                   "genome, so blocked by the same thing"}


# ---------------------------------------------------------------- liftover
class Lift:
    """Parental -> modified coordinates, rebuilt exactly from the saved table.

    The table stores (contig, parental_pos0, parental_end0, cum_delta). The
    replacement length is recoverable: for edit i, delta_i = cum_i - cum_(i-1)
    and ref_len_i = end_i - start_i, so altlen_i = ref_len_i + delta_i. That
    makes the reconstruction exact rather than approximate, including for
    equal-length replacements such as a start-codon lesion.
    """

    def __init__(self, rows):
        rows = sorted(rows)
        self.starts, self.ends, self.altlen, self.cum = [], [], [], []
        prev = 0
        for s, e, c in rows:
            self.starts.append(s)
            self.ends.append(e)
            self.altlen.append((e - s) + (c - prev))
            self.cum.append(c)
            prev = c

    def delta_before(self, pos):
        i = bisect.bisect_right(self.ends, pos) - 1
        return self.cum[i] if i >= 0 else 0

    def to_modified(self, pos):
        i = bisect.bisect_right(self.starts, pos) - 1
        if i >= 0 and self.starts[i] <= pos < self.ends[i]:
            k = pos - self.starts[i]
            if k < self.altlen[i]:
                return self.starts[i] + self.delta_before(self.starts[i]) + k
            return None
        return pos + self.delta_before(pos)


def load_lift(rid):
    p = HERE / f"results/liftover/{rid}.tsv"
    if not p.exists():
        return None
    by = collections.defaultdict(list)
    with p.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            by[r["contig"]].append((int(r["parental_pos0"]),
                                    int(r["parental_end0"]),
                                    int(r["cum_delta"])))
    return {c: Lift(v) for c, v in by.items()}


def c4_lift(acc):
    """C4's liftover, recomputed deterministically.

    tech_error_replicate seeds random.Random(f"{acc}|C4|v1") and deletes
    int(len * 0.0005) single positions per contig, in contig order. The
    simulator never wrote a liftover for this arm, so without reproducing the
    draw the technical-error arm could not be scored per gene at all. The
    reconstruction is verified against the written FASTA before use.
    """
    contigs, *_ = S.load_parental(acc)
    if not contigs:
        return None, None
    rng = random.Random(f"{acc}|C4|v1")
    per = {}
    for cg, s in contigs.items():
        k = int(len(s) * S.TECH_ERROR_RATE)
        pos = sorted(rng.sample(range(len(s)), k))
        run, rows = 0, []
        for p in pos:
            run -= 1
            rows.append((p, p + 1, run))
        per[cg] = Lift(rows)
    return per, contigs


# ---------------------------------------------------------------- scoring
def classes_for(idx, contig, s0, e0):
    return PC.classify(idx.overlapping(contig, s0, e0),
                       idx.fragments(contig, s0, e0)), \
        bool(idx.any_overlap(contig, s0, e0))


def main() -> None:
    truth = list(csv.DictReader((HERE / "results/ground_truth.tsv").open(),
                                delimiter="\t"))
    elig = list(csv.DictReader((HERE / "results/eligible_genes.tsv").open(),
                               delimiter="\t"))
    by_gene = collections.defaultdict(dict)   # (acc, tag) -> row
    for g in elig:
        by_gene[g["accession"]][g["locus_tag"]] = g
    lesions = collections.defaultdict(dict)   # rid -> tag -> truth row
    for r in truth:
        lesions[r["replicate"]][r["locus_tag"]] = r
    reps = sorted(lesions)
    print(f"{len(truth)} lesions, {len(elig)} eligible genes, "
          f"{len(reps)} lesioned replicates")
    print(f"callers evaluated: {CALLERS}")
    for k, v in BLOCKED.items():
        print(f"  BLOCKED {k}: {v}")

    rows = []
    for caller in CALLERS:
        for rid in reps:
            acc = rid.split("__")[0]
            par_id = f"{acc}__C0_parental"
            pdir = CALLS / caller / par_id
            mdir = CALLS / caller / rid
            if not pdir.exists() or not mdir.exists():
                print(f"  missing output: {caller} {rid}")
                continue
            pidx = PC.Index(PC.load_features(caller, pdir, par_id))
            midx = PC.Index(PC.load_features(caller, mdir, rid))
            lift = load_lift(rid)
            if lift is None:
                print(f"  no liftover for {rid}; skipped")
                continue
            for tag, g in by_gene[acc].items():
                s1, e1 = int(g["start1"]), int(g["end1"])
                cg = g["contig"]
                # parental class, in parental coordinates
                p_cls, p_any = classes_for(pidx, cg, s1 - 1, e1)
                lf = lift.get(cg)
                if lf is None:
                    continue
                m0, m1 = lf.to_modified(s1 - 1), lf.to_modified(e1 - 1)
                if m0 is None or m1 is None or m1 <= m0:
                    continue
                m_cls, m_any = classes_for(midx, cg, m0, m1 + 1)
                les = lesions[rid].get(tag)
                rows.append({
                    "caller": caller, "replicate": rid, "accession": acc,
                    "locus_tag": tag,
                    "lesion_code": les["lesion_code"] if les else "",
                    "truth_tier": les["truth_tier"] if les else "negative_untouched",
                    "position": les["position"] if les else "",
                    "parental_class": p_cls, "modified_class": m_cls,
                    "parental_any": int(p_any), "modified_any": int(m_any)})

    out = HERE / "results/stage_a_predictions.tsv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out.relative_to(HERE.parent)} with {len(rows)} gene rows")

    # ------------------------------------------------------------ confusion
    def positive(cls, any_feat, ep):
        if ep == "C":
            return bool(any_feat)
        return PC.endpoint_positive(cls, ep)

    report = {"callers_evaluated": CALLERS, "blocked": BLOCKED,
              "n_gene_rows": len(rows), "endpoints": {}}

    # ------------------------------------------------- endpoint C, unpaired
    # C asks whether the gene was RECOVERED as a feature. That is not a
    # detection event and must not be paired against the parental: an intact
    # gene is recovered in the parental too, so "positive in modified and not
    # in parental" is empty by construction and gave sensitivity 0.000 for
    # both callers in a first run. C is a rate, reported directly.
    print("\n" + "=" * 74)
    print("endpoint C: recovery rate, UNPAIRED (see the note in the source)")
    report["endpoint_C_recovery"] = {}
    for caller in CALLERS:
        sub = [r for r in rows if r["caller"] == caller]
        les = [r for r in sub if r["truth_tier"] in ("high", "medium")]
        unt = [r for r in sub if r["truth_tier"] == "negative_untouched"]
        rl = sum(r["modified_any"] for r in les) / len(les) if les else None
        ru = sum(r["modified_any"] for r in unt) / len(unt) if unt else None
        rp = sum(r["parental_any"] for r in les) / len(les) if les else None
        report["endpoint_C_recovery"][caller] = {
            "n_lesioned": len(les), "recovered_lesioned": rl,
            "n_untouched": len(unt), "recovered_untouched": ru,
            "recovered_lesioned_in_parental": rp}
        print(f"  {caller:6s} lesioned genes recovered {rl:.3f} "
              f"(n {len(les)}), in the parental {rp:.3f};  "
              f"untouched recovered {ru:.3f} (n {len(unt)})")

    print("\n" + "=" * 74)
    print("PAIRED confusion matrices for A and B. A gene counts positive only")
    print("if the caller calls it in the modified replicate and did NOT in the")
    print("parental.")
    for ep in ("A", "B"):
        report["endpoints"][ep] = {}
        print(f"\nendpoint {ep}")
        for caller in CALLERS:
            sub = [r for r in rows if r["caller"] == caller]
            # ambiguous genes are excluded from that caller's primary analysis
            amb = sum(1 for r in sub
                      if "ambiguous" in (r["parental_class"], r["modified_class"]))
            sub = [r for r in sub
                   if "ambiguous" not in (r["parental_class"], r["modified_class"])]
            tp = fp = fn = tn = 0
            for r in sub:
                mp = positive(r["modified_class"], r["modified_any"], ep)
                pp = positive(r["parental_class"], r["parental_any"], ep)
                gained = mp and not pp
                lesioned = r["truth_tier"] in ("high", "medium")
                if lesioned:
                    tp += gained
                    fn += not gained
                elif r["truth_tier"] == "negative_untouched":
                    fp += gained
                    tn += not gained
            sens = tp / (tp + fn) if tp + fn else None
            spec = tn / (tn + fp) if tn + fp else None
            prec = tp / (tp + fp) if tp + fp else None
            f1 = (2 * prec * sens / (prec + sens)
                  if prec and sens and prec + sens else None)
            bal = ((sens + spec) / 2) if sens is not None and spec is not None else None
            fpr = fp / (fp + tn) if fp + tn else None
            report["endpoints"][ep][caller] = {
                "TP": tp, "FP": fp, "FN": fn, "TN": tn,
                "sensitivity": sens, "specificity": spec, "precision": prec,
                "f1": f1, "balanced_accuracy": bal, "fpr_untouched": fpr,
                "n_ambiguous_excluded": amb}
            def f(x):
                return "   n/a" if x is None else f"{x:6.3f}"
            print(f"  {caller:6s} TP {tp:4d} FP {fp:4d} FN {fn:4d} TN {tn:5d}"
                  f"  sens {f(sens)}  spec {f(spec)}  prec {f(prec)}"
                  f"  FPR {f(fpr)}  ambiguous excluded {amb}")

    # ------------------------------------------------- by lesion class
    print("\n" + "=" * 74)
    print("endpoint A and B sensitivity by lesion class (high+medium tiers)")
    bycode = {}
    codes = sorted({r["lesion_code"] for r in rows
                    if r["truth_tier"] in ("high", "medium")})
    print(f"  {'lesion':18s} {'n':>5s} " +
          "  ".join(f"{c:>12s}" for c in CALLERS))
    for code in codes:
        line, cell = f"  {code:18s}", {}
        n_shown = 0
        for caller in CALLERS:
            sub = [r for r in rows if r["caller"] == caller
                   and r["lesion_code"] == code
                   and r["truth_tier"] in ("high", "medium")
                   and "ambiguous" not in (r["parental_class"], r["modified_class"])]
            n_shown = len(sub)
            a = sum(1 for r in sub
                    if positive(r["modified_class"], r["modified_any"], "A")
                    and not positive(r["parental_class"], r["parental_any"], "A"))
            b = sum(1 for r in sub
                    if positive(r["modified_class"], r["modified_any"], "B")
                    and not positive(r["parental_class"], r["parental_any"], "B"))
            cell[caller] = {"n": len(sub),
                            "A": a / len(sub) if sub else None,
                            "B": b / len(sub) if sub else None}
        bycode[code] = cell
        line += f"{n_shown:5d} " + "  ".join(
            f"A {cell[c]['A']:.2f} B {cell[c]['B']:.2f}" if cell[c]["A"] is not None
            else "        n/a " for c in CALLERS)
        print(line)
    report["by_lesion_class"] = bycode

    # ------------------------------------------------- class distribution
    print("\n" + "=" * 74)
    print("what the callers actually emit for lesioned genes (modified class)")
    for caller in CALLERS:
        sub = [r for r in rows if r["caller"] == caller
               and r["truth_tier"] in ("high", "medium")]
        n = collections.Counter(r["modified_class"] for r in sub)
        print(f"  {caller:6s} " + "  ".join(f"{k} {v}" for k, v in n.most_common()))
        report.setdefault("modified_class_distribution", {})[caller] = dict(n)

    # ------------------------------------------------- B1, the plan's gate
    print("\n" + "=" * 74)
    print("B1 GATE, exactly as SIMULATION_PLAN.md section 9 fixes it:")
    print("  every executed caller detects >= 50 % of HIGH-tier M1_stop_sub")
    print("  under endpoint B. If NO caller clears it, the plan says the")
    print("  pipeline or parser is broken rather than the callers, the run is")
    print("  void and nothing downstream may be reported.")
    gate = {}
    for caller in CALLERS:
        sub = [r for r in rows if r["caller"] == caller
               and r["lesion_code"] == "M1_stop_sub"
               and r["truth_tier"] == "high"
               and "ambiguous" not in (r["parental_class"], r["modified_class"])]
        hit = sum(1 for r in sub
                  if positive(r["modified_class"], r["modified_any"], "B")
                  and not positive(r["parental_class"], r["parental_any"], "B"))
        rate = hit / len(sub) if sub else None
        gate[caller] = {"n": len(sub), "detected": hit, "rate": rate}
        print(f"  {caller:6s} {hit}/{len(sub)} = "
              f"{'n/a' if rate is None else f'{rate:.3f}'}"
              f"  {'clears' if rate and rate >= 0.5 else 'does NOT clear'}")
    cleared = [c for c, v in gate.items() if v["rate"] and v["rate"] >= 0.5]
    report["B1_gate"] = {"per_caller": gate, "callers_clearing": cleared,
                         "gate_fires": len(cleared) == 0}
    print(f"  -> callers clearing the bar: {cleared if cleared else 'NONE'}")
    print(f"  -> gate fires: {len(cleared) == 0}")

    # ------------------------------------------------- by position
    print("\n" + "=" * 74)
    print("endpoint A/B detection by lesion POSITION, high+medium tiers")
    bypos = {}
    for pos in ("early", "middle", "late", "start", "mixed"):
        cells = {}
        for caller in CALLERS:
            sub = [r for r in rows if r["caller"] == caller
                   and r["position"] == pos
                   and r["truth_tier"] in ("high", "medium")
                   and "ambiguous" not in (r["parental_class"], r["modified_class"])]
            if not sub:
                continue
            a = sum(1 for r in sub
                    if positive(r["modified_class"], r["modified_any"], "A")
                    and not positive(r["parental_class"], r["parental_any"], "A"))
            b = sum(1 for r in sub
                    if positive(r["modified_class"], r["modified_any"], "B")
                    and not positive(r["parental_class"], r["parental_any"], "B"))
            cells[caller] = {"n": len(sub), "A": a / len(sub), "B": b / len(sub)}
        if cells:
            bypos[pos] = cells
            print(f"  {pos:8s} " + "  ".join(
                f"{c}: n {v['n']:3d} A {v['A']:.2f} B {v['B']:.2f}"
                for c, v in cells.items()))
    report["by_position"] = bypos

    (HERE / "results/stage_a_evaluation.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("\nwrote benchmark/results/stage_a_evaluation.json")


if __name__ == "__main__":
    main()
