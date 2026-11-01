#!/usr/bin/env python3
"""Bound the disk and runtime footprint of reciprocal-best-hit orthology.

WHY THIS EXISTS. FINAL_INDEPENDENT_VALIDATION.md ranks unverified orthology as
the largest standing assumption in the dN/dS analysis, and says the obvious
remedy - reciprocal best hits with synteny confirmation - was out of scope
"because it requires an all-vs-all DIAMOND pass whose disk footprint has not
been bounded against the 10 GB guard". That sentence contains an assumption
("all-vs-all") that was never checked. This script checks it, by measuring one
genome's worth of RBH work end to end and reporting what it costs.

WHAT RBH ACTUALLY NEEDS, per (query genome, partner genome) pair:
  forward   locus nucleotides  -> partner proteome   (already run; blastx)
  reverse   the matched partner protein -> QUERY genome's own proteome (blastp)
  accept    iff the reverse best hit is the protein annotated at the locus
That is two directed searches per pair, not all-vs-all across the corpus. The
corpus has 92 eligible (query, partner) pairs, so the work is linear in pairs.

WHAT THIS SCRIPT DOES NOT DO. It does not run RBH and it does not touch any
validated output. It builds the query-side proteome and index for a few
representative genomes, runs one reverse search at realistic size, measures
every byte and second, and extrapolates. Everything it writes goes under
results/rbh_probe/ and the large transient artefacts are removed at the end.

CONCURRENCY NOTE. A full-corpus reference-sensitivity run may be in progress and
owns ~/research/datasets/n2-refseq-reference/. This script treats that directory
as READ ONLY and writes its own files elsewhere, so the two cannot race. It also
builds the query-side proteome itself rather than calling
analysis_reference_sensitivity.build_refseq_reference, for the same reason.

DISK GUARD. Refuses to start below 10 GB free on /mnt/c, and aborts mid-run if
free space falls below 8 GB. Host free space is read from /mnt/c, never from
`df /`, which reports WSL's virtual capacity and is not evidence of physical
space.
"""
import json
import pathlib
import shutil
import subprocess
import sys
import time

from Bio import SeqIO

P = pathlib.Path.home() / "research/projects/n2-pseudogene-caller-bias"
D = pathlib.Path.home() / "research/datasets"
R = P / "results"
OUT = R / "rbh_probe"
REFDIR = D / "n2-refseq-reference"             # READ ONLY here

FLOOR_START_GB = 10
FLOOR_ABORT_GB = 8
THREADS = 2          # deliberately low: a full-corpus run may hold the CPU

# One small, one mid, one large genome, from three strata, all in the
# 92-target manifest.
TARGETS = [
    ("GCF_024338845.1", "GCF_000723365.1", "Mycoplasma, smallest proteome"),
    ("GCF_004008975.1", "GCF_020735925.1", "Bordetella, mid"),
    ("GCF_000012005.1", "GCF_000006925.2", "Shigella, most loci (802 discordant)"),
]


# --------------------------------------------------------------- disk guard
def host_free_gb():
    """Host C: free space, from /mnt/c. Never `df /`."""
    out = subprocess.run(["df", "--output=avail", "-BG", "/mnt/c"],
                         capture_output=True, text=True, check=True).stdout
    return int("".join(c for c in out.splitlines()[-1] if c.isdigit()))


def guard(stage):
    free = host_free_gb()
    print(f"  [guard] {stage}: host C: free {free} GB")
    if free < FLOOR_ABORT_GB:
        sys.exit(f"EMERGENCY STOP: {free} GB < {FLOOR_ABORT_GB} GB floor")
    return free


def du_bytes(p: pathlib.Path):
    if not p.exists():
        return 0
    if p.is_file():
        return p.stat().st_size
    return sum(f.stat().st_size for f in p.rglob("*") if f.is_file())


# --------------------------------------------------------------- proteomes
def gbff_path(acc):
    g = next((D / "n2-genomes/extracted").rglob(f"*{acc}*/genomic.gbff"), None)
    if g is None and (D / f"n2-genomes-expansion/{acc}").exists():
        g = next((D / f"n2-genomes-expansion/{acc}").rglob("genomic.gbff"), None)
    return g


def query_proteome(acc, dest: pathlib.Path):
    """The QUERY genome's own non-pseudo CDS translations.

    This is the reverse-search database. Written standalone rather than through
    analysis_reference_sensitivity.build_refseq_reference, so that a concurrent
    full-corpus run cannot race this script on the same cache files.
    """
    g = gbff_path(acc)
    if g is None:
        return None, 0
    n = 0
    with dest.open("w") as fh:
        for rec in SeqIO.parse(str(g), "genbank"):
            for ft in rec.features:
                if ft.type != "CDS":
                    continue
                q = ft.qualifiers
                if "pseudo" in q or "pseudogene" in q:
                    continue
                tr = q.get("translation", [None])[0]
                if not tr:
                    continue
                ident = (q.get("locus_tag") or q.get("protein_id")
                         or [f"{rec.id}_{int(ft.location.start)}"])[0]
                fh.write(f">{ident}\n{tr}\n")
                n += 1
    return dest, n


def sample_partner_proteins(partner, k, dest: pathlib.Path):
    """k partner proteins, as the reverse search's query set.

    In the real analysis these are the forward search's best hits - one protein
    per locus. Here any k of them size the search identically, because DIAMOND's
    cost is driven by query count and length, not by which proteins they are.
    """
    src = REFDIR / f"{partner}.faa"
    if not src.exists():
        return None, 0
    n = 0
    with dest.open("w") as fh:
        for rec in SeqIO.parse(str(src), "fasta"):
            fh.write(f">{rec.id}\n{rec.seq}\n")
            n += 1
            if n >= k:
                break
    return dest, n


def run(cmd):
    t0 = time.time()
    subprocess.run(cmd, check=True)
    return time.time() - t0


# --------------------------------------------------------------- main
def main():
    print("RBH FOOTPRINT PROBE — measures the cost, does not run RBH\n")
    free0 = guard("start")
    if free0 < FLOOR_START_GB:
        sys.exit(f"refusing to start: {free0} GB < {FLOOR_START_GB} GB floor")

    if not REFDIR.exists():
        sys.exit(f"missing partner references at {REFDIR}; run "
                 "analysis_reference_sensitivity.py first")
    if shutil.which("diamond") is None:
        sys.exit("diamond not on PATH")

    OUT.mkdir(parents=True, exist_ok=True)
    # loci per genome, to size the reverse query set realistically
    import csv
    nloci = {}
    with open(P / "manifests/refsens_all.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            nloci[row["accession"]] = int(row["n_dis"]) + int(row["n_con"])

    rows = []
    for acc, partner, label in TARGETS:
        print(f"\n=== {acc} vs {partner}  ({label})")
        work = OUT / acc
        work.mkdir(exist_ok=True)
        faa = work / "query_proteome.faa"

        t0 = time.time()
        _, n_prot = query_proteome(acc, faa)
        t_extract = time.time() - t0
        if not n_prot:
            print("  no proteome; skipped")
            continue
        print(f"  query proteome      {n_prot:6d} proteins  "
              f"{faa.stat().st_size / 1e6:7.2f} MB  {t_extract:5.1f} s")

        dmnd = work / "query_proteome.dmnd"
        t_db = run(["diamond", "makedb", "--in", str(faa), "-d", str(dmnd),
                    "--quiet", "--threads", str(THREADS)])
        print(f"  DIAMOND index       {'':6s}            "
              f"{dmnd.stat().st_size / 1e6:7.2f} MB  {t_db:5.1f} s")

        k = min(nloci.get(acc, 500), 1500)
        qf = work / "reverse_queries.faa"
        _, n_q = sample_partner_proteins(partner, k, qf)
        print(f"  reverse query set   {n_q:6d} proteins  "
              f"{qf.stat().st_size / 1e6:7.2f} MB")

        hits = work / "reverse_hits.tsv"
        t_search = run(["diamond", "blastp", "-q", str(qf), "-d", str(dmnd),
                        "-o", str(hits), "--quiet", "--threads", str(THREADS),
                        "--max-target-seqs", "1", "--evalue", "1e-5",
                        "--outfmt", "6", "qseqid", "sseqid", "pident",
                        "length", "evalue", "bitscore"])
        n_hits = sum(1 for _ in hits.open())
        print(f"  reverse search      {n_hits:6d} hits      "
              f"{hits.stat().st_size / 1e6:7.2f} MB  {t_search:5.1f} s")

        peak = du_bytes(work)
        print(f"  PEAK on disk, this genome: {peak / 1e6:.2f} MB")
        guard(f"after {acc}")

        rows.append({
            "accession": acc, "partner": partner, "label": label,
            "n_query_proteins": n_prot, "n_reverse_queries": n_q,
            "n_reverse_hits": n_hits,
            "bytes_proteome": faa.stat().st_size,
            "bytes_index": dmnd.stat().st_size,
            "bytes_reverse_queries": qf.stat().st_size,
            "bytes_hits": hits.stat().st_size,
            "bytes_peak": peak,
            "sec_extract": round(t_extract, 2),
            "sec_makedb": round(t_db, 2),
            "sec_search": round(t_search, 2),
            "sec_total": round(t_extract + t_db + t_search, 2),
        })
        # keep the measurement, drop the bulk
        for f in (dmnd, qf, hits, faa):
            f.unlink(missing_ok=True)

    if not rows:
        sys.exit("no genome measured")

    # ------------------------------------------------- extrapolation
    n_pairs = len(nloci)
    peak_max = max(r["bytes_peak"] for r in rows)
    per_genome_persist = max(r["bytes_proteome"] + r["bytes_index"] for r in rows)
    sec_max = max(r["sec_total"] for r in rows)
    sec_mean = sum(r["sec_total"] for r in rows) / len(rows)

    bound = {
        "note": "POST HOC; footprint bound only, RBH itself not run",
        "n_eligible_pairs": n_pairs,
        "threads_used": THREADS,
        "measured": rows,
        "serial_peak_bytes": peak_max,
        "serial_peak_mb": round(peak_max / 1e6, 2),
        "cached_all_pairs_bytes": per_genome_persist * n_pairs,
        "cached_all_pairs_gb": round(per_genome_persist * n_pairs / 1e9, 3),
        "runtime_serial_hours_at_%d_threads" % THREADS:
            round(sec_mean * n_pairs / 3600, 2),
        "runtime_worst_genome_sec": sec_max,
        "host_free_gb_start": free0,
        "host_free_gb_end": host_free_gb(),
        "floor_start_gb": FLOOR_START_GB,
        "floor_abort_gb": FLOOR_ABORT_GB,
    }

    print("\n" + "=" * 68)
    print("BOUND")
    print("=" * 68)
    print(f"  eligible (query, partner) pairs            {n_pairs}")
    print(f"  peak transient disk, serial, one genome    "
          f"{bound['serial_peak_mb']:.2f} MB")
    print(f"  if every pair's proteome+index is cached   "
          f"{bound['cached_all_pairs_gb']:.3f} GB")
    print(f"  serial runtime at {THREADS} threads            "
          f"{bound['runtime_serial_hours_at_%d_threads' % THREADS]:.2f} h")
    print(f"  host C: free, start -> end                 "
          f"{free0} -> {bound['host_free_gb_end']} GB")

    headroom = bound["host_free_gb_end"] - FLOOR_ABORT_GB
    need_gb = max(bound["cached_all_pairs_gb"], peak_max / 1e9)
    verdict = ("WITHIN GUARD" if need_gb < headroom else "EXCEEDS GUARD")
    bound["verdict"] = verdict
    bound["headroom_gb_above_abort_floor"] = headroom
    bound["worst_case_need_gb"] = round(need_gb, 3)
    print(f"\n  worst-case need {need_gb:.3f} GB against {headroom} GB of "
          f"headroom above the {FLOOR_ABORT_GB} GB floor  ->  {verdict}")

    (R / "rbh_footprint.json").write_text(json.dumps(bound, indent=1,
                                                     sort_keys=True))
    print(f"\nwrote {R / 'rbh_footprint.json'}")

    for f in ("loci.tsv", "locus.json"):
        print(f"  {f} untouched by this script (never opened for writing)")


if __name__ == "__main__":
    main()
