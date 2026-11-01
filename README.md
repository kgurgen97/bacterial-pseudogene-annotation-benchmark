# Bacterial pseudogene annotation: methods and reproducibility

Research code, result tables, source figures and supplementary materials for a comparison of bacterial pseudogene annotation configurations.

**Authors:** Gurgen Karapetyan¹˒² and Astghik Pepoyan¹

¹ Armenian National Agrarian University, Yerevan, Armenia  
² Scientific Center of Zoology and Hydroecology, Yerevan, Armenia

**Corresponding author:** Gurgen Karapetyan  
**Funding:** grant 25AA-4A009 (support details to be verified by the authors before journal submission).

## Study design

We compare four annotation configurations on the same **523 complete RefSeq bacterial genomes**: PGAP, Pseudofinder applied to PGAP-derived annotations, Bakta, and DFAST. Pseudofinder is **not an independent fourth annotation pipeline**. Corpus-total pseudogene counts span **19.1×**, and unrounded median per-genome pseudogene fractions span **77.2×**.

| Analysis | Population | Principal files |
|---|---|---|
| Pseudogene abundance | 523 genomes | `results/corpus523.tsv`, `supplementary/S1.tsv` |
| Genome-size association | 363 annotated; 362 eligible for direct paired contrasts (185 genera) | `results/literature_retest.tsv`, `results/caller_effect_contrast.json` |
| Locus disagreement and dN/dS | 92 genomes, 28 genera, 103,224 scored loci, 4,886 disputed | `results/loci.tsv`, `results/locus.json` |
| Reference and RBH checks | separate 36,535-locus subset | `results/reference_sensitivity_all.json`, `results/rbh_audited_all.json` |
| Sequence-lesion pilot | 12 genomes, 6,360 loci | `results/lesion_extended.json`, `supplementary/S7.tsv` |

## Files

- `manuscript/MANUSCRIPT.md`: scientific manuscript text with limitations and author information. The submitted Word version may be formatted separately to meet journal requirements.
- `figures/figure1_study_design.*` through `figure5_lesion_evidence.*`: original full-panel figures (PDF, PNG and SVG), not redrawn approximations.
- `figures/FIGURE_CAPTIONS.md`: figure legends; `figures/data/` contains plotting data.
- `supplementary/S1.tsv`–`S9.tsv` and `supplementary/S7b.tsv`: complete supplied supplementary tables.
- `analysis_*.py` and `benchmark/`: analysis and testing code, statistical results, frozen research plans, and validation records.
- `REPRODUCIBILITY.md`, `DATA_AVAILABILITY.md`, and `PUBLICATION_CHECKLIST.md`: scope, data sources, and outstanding publication requirements.

## Main findings and limitations

Direct paired comparisons find statistically different associations for two of three independent configuration pairs using **raw counts**; no such pair is distinguishable using fraction or counts/Mb. Failure to detect a difference does **not** demonstrate equivalence.

Disputed loci have higher dN/dS at genus level (25/28 genera; **two-sided** exact sign-test p = 2.744 × 10⁻⁵), but simulated frameshifts mechanically affect this estimator. The +0.0315 median frameshift effect is imprecise (11 genes; 5/11 leave-one-out medians fall below the prespecified +0.030 threshold). RBH retains 67.5% of disputed and 97.2% of concordant-intact loci; this differential retention may bias sensitivity analyses. The overlap sweep is conditional on previously scored loci, and the 12-genome lesion pilot is not independent biological ground truth. We make **no ranking of callers by biological accuracy** and no claim about loss of gene function.

## Large intermediate data

The publication-facing GitHub snapshot omits twelve large intermediate TSV artifacts required by some full reanalysis workflows. The original scientific results are unchanged; see `LARGE_INTERMEDIATE_DATA.md` for exact names, sizes, and reproducibility limitations.

## Reproducibility

```bash
python -m pip install -r requirements.txt
python -m pytest -q benchmark/tests
python figures/make_figures.py
```

Figure regeneration and tests use supplied artifacts. Complete genome-wide annotation reruns require the accessioned assemblies, external tools and larger intermediate data not included in this distribution. See `REPRODUCIBILITY.md`.

## Citation and licence

The article DOI has not yet been assigned; see `CITATION.cff`. Code is MIT licensed; manuscript, figures and derived tables are covered by the terms recorded in `LICENSE`, with third-party inputs retaining their own terms.
