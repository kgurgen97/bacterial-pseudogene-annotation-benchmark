# Reproducibility and scope

This is a publication snapshot, not a newly executed biological analysis. Source results and original figure files were preserved from the research repository. No genome-wide annotation analysis was rerun during publication preparation.

## Contents

- `results/` and `benchmark/results/`: committed numerical summaries and available per-locus evidence; the twelve omitted large intermediate tables are itemized in `LARGE_INTERMEDIATE_DATA.md`.
- `supplementary/`: all complete S1–S9 tables and S7b.
- `figures/`: five original figures in PNG, PDF and SVG; `figures/data/` contains figure inputs.
- `analysis_*.py`, `benchmark/`, and original frozen plans: analysis and quality checks.
- `results/corpus523.tsv`: reference assembly accessions.

**Important:** The GitHub snapshot omits 12 large intermediate TSV files, including `results/loci.tsv` and `results/omega_pool.tsv`. See `LARGE_INTERMEDIATE_DATA.md`; these must be obtained from the research archive for scripts that require them. Do not interpret their absence as successful full reproducibility from this checkout.

## Environment

Install baseline Python dependencies using `python -m pip install -r requirements.txt`. Testing and plot regeneration:
```bash
python -m pytest -q benchmark/tests
python figures/make_figures.py
```

Reproducing the entire pipeline from raw sequence requires external reference assemblies, exact original annotation output files, larger private/local intermediates and tools including PGAP annotations, DIAMOND, Pseudofinder, Bakta, and DFAST. Tool versions and experiment parameters are specified in the study materials where available. We do not claim that a fresh workstation can run the entire source pipeline directly from this repository.

## Statistical safeguards

- The 523-genome descriptive corpus is distinct from the 363 trait-annotated subset and the 362-genome direct paired contrast set.
- The 92-genome/28-genus locus analysis comprises 103,224 scored loci, with a separate 36,535-locus reference-sensitivity subset. Do not combine their denominators.
- The genus-level sign test is **two-sided** (25/28, p = 2.744048833847046e-5).
- The 11-gene frameshift-arm magnitude is uncertain; 5/11 leave-one-out estimates cross the planned +0.030 band edge.
- Overlap-threshold sweeps reused existing omega values and are conditional on the previously scored set.
- Differential RBH retention (67.5% disputed, 97.2% consensus-intact) can induce selection bias; failed gates must not be described as passed.
- The 12-genome lesion-evidence pilot is neither independent functional validation nor proof that a given pseudogene is nonfunctional.

Frozen plans are left unchanged to preserve their original SHA-256 provenance. Their failures and deviations remain part of the scientific record. An independent external biological ground truth set was not available for caller accuracy ranking.
