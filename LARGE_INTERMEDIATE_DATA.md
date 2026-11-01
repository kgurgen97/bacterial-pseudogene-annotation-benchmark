# Large intermediate datasets

This repository provides the complete supplementary tables, original publication figures and their plotting data, committed result summaries, analysis scripts and test suite.

The following **12 large intermediate TSV files** were not included in this GitHub snapshot because they could not be transferred through the current GitHub file API. They remain in the original private research archive and are **not** to be mistaken for missing supplementary tables:

| Intermediate result | Source size (bytes) |
|---|---:|
| `results/dnds_pairs.tsv` | 5040702 |
| `results/frame_pairs.tsv` | 5382833 |
| `results/loci.tsv` | 15500196 |
| `results/omega_pool.tsv` | 49254050 |
| `results/overlap_sensitivity/loci_recip_0.30.tsv` | 12560569 |
| `results/overlap_sensitivity/loci_recip_0.40.tsv` | 12611423 |
| `results/overlap_sensitivity/loci_recip_0.50.tsv` | 12624354 |
| `results/overlap_sensitivity/loci_recip_0.60.tsv` | 12570534 |
| `results/overlap_sensitivity/loci_recip_0.70.tsv` | 12518352 |
| `results/rbh3_all/input_snapshot.tsv` | 6379359 |
| `results/reference_sensitivity/per_locus_all.tsv` | 6379359 |
| `results/rbh3_all/per_locus_rbh_audited.tsv` | 14493599 |

**Consequence:** Some end-to-end reanalysis scripts require these files and cannot run from the GitHub checkout alone. The full scientific evidence should be archived with a persistent DOI before journal submission. This GitHub snapshot is a publication-facing distribution, **not a complete raw-data-to-results rerun bundle**.

The manuscript findings are not changed by omission of large intermediates from this distribution. Their source values and analyses remain in the research archive.
