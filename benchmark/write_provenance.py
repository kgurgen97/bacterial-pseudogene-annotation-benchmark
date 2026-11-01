#!/usr/bin/env python3
"""Provenance for this session: versions, inputs, checksums, commands.

SIMULATION_PLAN.md and the session brief both require that a computation be
reproducible from the record rather than from memory. This writes one manifest:
software versions actually used, the input stores with their sizes, SHA-256 for
every committed result and every frozen plan, and the command that produced
each analysis.

Raw assemblies and caller output are deliberately NOT checksummed file by file -
they are 1.3 GB of regenerable intermediates. Their directories and sizes are
recorded instead, with the command to rebuild them.

    python3 benchmark/write_provenance.py
"""
import hashlib
import json
import pathlib
import platform
import subprocess
import sys

P = pathlib.Path(__file__).resolve().parent.parent
D = pathlib.Path.home() / "research/datasets"

PLANS = ["DNDS_PLAN.md", "FRAME_PLAN.md", "FRAME2_PLAN.md", "LOCUS_PLAN.md",
         "PHYLO_PLAN.md", "LITERATURE_RETEST_PLAN.md",
         "CORPUS_EXPANSION_PLAN.md", "PF_CONFIG_PLAN.md",
         "benchmark/SIMULATION_PLAN.md"]

COMMANDS = {
    "classifier fix verification": "python3 benchmark/check_classifier_fix.py 523",
    "parser audit": "python3 benchmark/audit_fixed_parser.py 523",
    "RNG coupling": "python3 benchmark/check_rng_coupling.py",
    "dN/dS, clustering, locus, dfast, frame": "bash benchmark/rerun_after_fix.sh",
    "omega pool": "bash benchmark/run_pool.sh",
    "seed sweep": "bash benchmark/run_sweep.sh",
    "locus sensitivity": "python3 analysis_locus_sensitivity.py",
    "G3 contig sensitivity": "python3 analysis_g3_contig_sensitivity.py",
    "G3 by lesion class": "python3 analysis_g3_by_lesion.py",
    "caller comparison": "python3 analysis_caller_comparison.py",
    "frame2": "python3 analysis_frame2.py",
    "frame2 confound": "python3 benchmark/check_frame2_confound.py",
    "simulation": "python3 benchmark/simulate.py --stage A --replicates 2 "
                  "--accessions GCF_000006925.2 GCF_000007765.2 GCF_000018105.1",
    "Stage A caller runs": "bash benchmark/run_stage_a.sh",
    "Stage A evaluation": "python3 benchmark/evaluate_stage_a.py",
    "B1 diagnosis": "python3 benchmark/diagnose_b1_failure.py",
    "manuscript audit": "python3 benchmark/audit_manuscript_claims.py",
}


def sha(p):
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def run(cmd):
    try:
        return subprocess.run(cmd, shell=True, capture_output=True,
                              text=True, timeout=60).stdout.strip().splitlines()[0]
    except Exception:
        return "unavailable"


def main() -> None:
    env = {
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "kernel": run("uname -r"),
        "cpu": run("grep -m1 'model name' /proc/cpuinfo | cut -d: -f2"),
        "nproc": run("nproc"),
        "mem_total_gb": run("free -g | awk '/^Mem:/{print $2}'"),
        "host_disk_free": run("df -h /mnt/c | tail -1 | awk '{print $4}'"),
    }
    callers = {
        "bakta": run("~/miniforge3/envs/callers/bin/bakta --version"),
        "dfast": run("~/miniforge3/envs/callers/bin/dfast --version"),
        "diamond": run("~/.local/bin-analysis/diamond --version"),
        "pseudofinder_commit": run(
            "git -C tools/pseudofinder log --oneline -1"),
    }
    pyver = {}
    for m in ("numpy", "pandas", "scipy", "Bio"):
        pyver[m] = run(f"python3 -c \"import {m}; "
                       f"print(getattr({m}, '__version__', '?'))\"")

    plans = {}
    for f in PLANS:
        p = P / f
        if p.exists():
            plans[f] = sha(p)

    results = {}
    for p in sorted((P / "results").glob("*.json")):
        results[f"results/{p.name}"] = sha(p)
    for p in sorted((P / "benchmark/results").glob("*.json")):
        results[f"benchmark/results/{p.name}"] = sha(p)
    tables = {}
    for rel in ("results/loci.tsv", "results/dnds_pairs.tsv",
                "results/frame_pairs.tsv", "results/omega_pool.tsv",
                "results/frame2_genes.tsv",
                "benchmark/results/ground_truth.tsv",
                "benchmark/results/eligible_genes.tsv",
                "benchmark/results/stage_a_predictions.tsv"):
        p = P / rel
        if p.exists():
            tables[rel] = {"sha256": sha(p),
                           "bytes": p.stat().st_size,
                           "rows": sum(1 for _ in p.open()) - 1}

    stores = {}
    for name in ("n2-genomes", "n2-genomes-expansion", "n2-bakta", "n2-dfast",
                 "n2-pseudofinder", "n2-benchmark", "bakta-db"):
        d = D / name
        if d.exists():
            stores[name] = {
                "path": str(d),
                "size": run(f"du -sh {d} 2>/dev/null | cut -f1"),
                "entries": len(list(d.iterdir()))}

    git = {
        "head": run("git rev-parse HEAD"),
        "branch": run("git rev-parse --abbrev-ref HEAD"),
        "remote": run("git config --get remote.origin.url"),
        "dirty": run("git status --porcelain | wc -l"),
    }

    out = {"generated": run("date -u +%Y-%m-%dT%H:%M:%SZ"),
           "environment": env, "callers": callers, "python_packages": pyver,
           "frozen_plan_sha256": plans,
           "result_json_sha256": results,
           "result_tables": tables,
           "data_stores": stores,
           "git": git,
           "commands": COMMANDS,
           "note": "Raw assemblies and caller output are regenerable "
                   "intermediates and are not checksummed individually; their "
                   "stores and the commands that rebuild them are recorded."}
    (P / "benchmark/results/provenance.json").write_text(
        json.dumps(out, indent=2, sort_keys=True) + "\n")

    print(f"environment: {env['platform']}, python {env['python']}, "
          f"{env['nproc']} cpus, {env['mem_total_gb']} GB")
    print(f"host disk free: {env['host_disk_free']}")
    print("callers:")
    for k, v in callers.items():
        print(f"  {k:22s} {v}")
    print("python packages:")
    for k, v in pyver.items():
        print(f"  {k:22s} {v}")
    print(f"frozen plans hashed: {len(plans)}")
    print(f"result json hashed : {len(results)}")
    print("result tables:")
    for k, v in tables.items():
        print(f"  {k:46s} {v['rows']:>8} rows  {v['bytes'] / 1e6:6.1f} MB")
    print(f"git {git['branch']} at {git['head'][:8]}, "
          f"{git['dirty']} uncommitted paths")
    print("\nwrote benchmark/results/provenance.json")


if __name__ == "__main__":
    main()
