# Expanding the corpus beyond two phyla, with the reduction gradient populated

**Written 2026-10-22, before any genome was selected or downloaded.** Availability was surveyed first
(section 3) because how many genomes exist is a feasibility fact, not an outcome. The selection rules
below are fixed and may not be changed because of a result.

## 1. Why

`LITERATURE_RETEST_RESULTS.md` has two weaknesses that this addresses, and one opportunity.

**Weakness 1 — two phyla.** The corpus is Pseudomonadota 259 and Actinomycetota 104. No claim about
bacteria in general can rest on that.

**Weakness 2 — no obligate endosymbionts.** C2 (host-associated versus free-living) rests on 74 genomes
whose labels come from species ecology, and the host-associated side contains no genuinely
host-restricted, genome-reduced lineage. The canonical claim is really about those lineages.

**The opportunity.** C1 came out *reversed*: pseudogene fraction rose with genome size (PGAP +0.126,
Pseudofinder +0.233). The explanation offered was that the literature's relationship is non-monotonic —
lineages mid-decay carry many pseudogenes, fully reduced endosymbionts carry few because the dead genes
are already deleted — and that this corpus samples only the free-living limb. **That explanation is
currently untested, because the corpus has neither mid-decay nor fully reduced lineages.** Adding both
turns a negative result into a testable prediction:

> **C3 (new, primary for the expansion):** pseudogene fraction is **non-monotonic** in genome reduction.
> It is low in free-living bacteria, **highest in recently host-restricted lineages**, and low again in
> extremely reduced obligate endosymbionts.

## 2. Sampling is stratified by ecology, never by genome size

This matters and is easy to get wrong. If genomes were chosen to fill a range of genome sizes, and then a
correlation between pseudogene fraction and genome size were computed on that sample, the deliberate
range expansion would inflate the correlation by construction. So:

**Strata are defined by lifestyle only. Genome size is never a selection criterion and is allowed to fall
wherever it falls.** C1 on the expanded corpus is reported alongside C1 on the original 363 so that any
change can be attributed.

| stratum | what it is for | target n | expected size |
|---|---|---|---|
| S1 extreme obligate endosymbiont | the reduced endpoint; predicted **few** pseudogenes | 25 | < 1 Mb |
| S2 obligate intracellular | mid-reduction | 35 | 1–2.5 Mb |
| S3 recent host restriction | mid-decay; predicted **many** pseudogenes | 30 | 2–5 Mb |
| S4 free-living breadth, new phyla | removes the two-phyla limit | 70 | unconstrained |
| | | **160** | |

160 is a compute budget, not a statistical one: Bakta runs at roughly 5 minutes per genome on this
machine, so 160 additions cost about 13 hours on top of the 42 hours still outstanding for the original
corpus. Stated so that the number is not mistaken for a power calculation.

**Lifestyle labels for the additions come from the stratum itself**, which is far stronger than the
Madin genus-modal fallback that limited C2: nothing is in doubt about *Buchnera* being an obligate
endosymbiont. The additions therefore also repair C2's label quality.

## 3. Availability (surveyed 2026-10-22, complete RefSeq, atypical excluded)

Abundant in every stratum, so the budget binds, not supply. Counts are genome counts, not species counts:
Buchnera 82, Carsonella 57, *Ca.* Sulcia 95, *Ca.* Portiera 14, *Ca.* Blochmannia 31; Rickettsia 163,
Wolbachia 329, Chlamydia 257, Mycoplasma 186, Coxiella 64; *Yersinia pestis* 234, Shigella 227,
*Bordetella pertussis* 715, *Salmonella* Typhi 143, *Francisella tularensis* 65; and thousands across
Bacillota, Bacteroidota, Campylobacterota and the rest.

*Ca.* Nasuia and *Ca.* Hodgkinia returned 0 complete RefSeq genomes and are dropped.

## 4. Selection rules — deterministic, fixed now

1. Query the NCBI Datasets v2 API: `assembly_level = complete_genome`, `assembly_source = refseq`,
   `exclude_atypical = true`.
2. **One genome per species.** This is the single most important rule: *Bordetella pertussis* alone has
   715 near-clonal assemblies, and without it one species would dominate every stratum. The original 363
   are one-per-species and the expansion must match.
3. Within a species, prefer the assembly flagged **reference**, then **representative**, then the
   lexicographically smallest accession. No randomness, so the set is reproducible from this file.
4. **Cap of 5 genomes per genus**, and for S4 a cap of 5 per phylum with at most one per family, so no
   single clade dominates a stratum.
5. Exclude anything already in the corpus of 363.
6. Require CheckM-style completeness as reported by NCBI where available; drop assemblies flagged as
   anomalous. Record what was dropped and why.

## 5. What is run on the additions

All four callers, with the counting rules already fixed and unchanged: PGAP from the RefSeq annotation,
Pseudofinder (defaults, Swiss-Prot), Bakta 1.12.2 (DB 6.0 full, its own `pseudogenes:` line), DFAST 1.5.3
via `dfast_runtime_fix.py` (disrupted classes only, never `partial`).

## 6. Analyses on the expanded corpus

- **C3 (new primary).** Pseudogene fraction across the four ecological strata, Kruskal–Wallis across
  strata then pairwise Mann–Whitney, per caller. Predicted ordering: S3 > S2 > S4 ≈ S1, i.e. a peak at
  S3 and a return to low values at S1. **Non-monotonicity is declared confirmed only if S3 > S1
  significantly and S3 > S4 significantly for at least three of four callers.**
- **C1 re-run**, reported side by side with the original-363 result. If the pooled correlation changes
  sign or significance, that is attributed to the strata now present, and both numbers are reported.
- **C2 re-run** with S1+S2+S3 as host-restricted and S4 plus the existing free-living labels as
  free-living. This is the properly labelled version of C2.
- **Caller dispersion** recomputed; the prediction from `analysis_caller_dispersion.py` is that
  Pseudofinder's relative dispersion stays low even across a far wider biological range, which would make
  the low-dispersion finding much stronger.

## 7. What this still will not fix

- Still RefSeq complete genomes only, so still biased toward cultured and medically interesting taxa.
- Still one configuration per caller.
- S1–S3 assignments are lineage-level ecology from the literature, not measured for each strain. They are
  far more secure than genus-modal labels but they are still assignments.
- Archaea are not included; the claim stays about bacteria.
- Adding reduced genomes widens the genome-size range, so **C1 on the expanded corpus is not directly
  comparable in magnitude to C1 on the original 363** — only the sign and significance are. This is
  stated because range expansion inflates correlation coefficients.

## 7a. AMENDMENT, 2026-10-22 — before any caller was run on any new genome

The first selection run delivered **140** genomes, not 160, and two strata under-filled:
**S1 = 16** (target 25) and **S3 = 19** (target 30). The cause is the one-genome-per-species rule working
as intended: *Buchnera*, *Carsonella* and *Ca.* Sulcia have many assemblies but few species, and the ten
S3 taxa together contain only 20 species.

Inspecting the selection also showed **S3 is heterogeneous in reduction stage**. It contains the classic
recently restricted pathogens whose genomes are not yet small (*Yersinia pestis*, *Shigella*,
*Bordetella pertussis*, *Salmonella* Typhi, ~4–5 Mb) alongside lineages already far reduced: a
Sodalis-like aphid symbiont at **0.47 Mb** and *Treponema pallidum* at **1.14 Mb**. With n = 19 that
spread would blur precisely the mid-decay peak C3 is meant to detect.

**Amendment: the S1 and S3 taxon lists are extended to reach the targets already declared in section 2.**
No target, threshold, selection rule or endpoint changes. The reason is a shortfall in supply and a
within-stratum heterogeneity visible from assembly metadata — **no pseudogene value has been computed for
any of these 140 genomes**, so this cannot be outcome-driven. The original list is preserved in commit
`17bd7d0`.

- **S3 additions** (recently host-restricted lineages with documented pseudogene accumulation):
  *Mycobacterium ulcerans*, *Bordetella parapertussis*, *Bordetella holmesii*,
  *Salmonella enterica* subsp. *enterica* serovar Paratyphi A, *Arsenophonus*, *Hamiltonella*,
  *Regiella*, *Candidatus* Riesia, *Baumannia*.
- **S1 additions** (further extremely reduced endosymbionts): *Candidatus* Moranella,
  *Candidatus* Profftella, *Candidatus* Evansia, *Candidatus* Vidania, *Candidatus* Purcelliella,
  *Candidatus* Uzinura, *Candidatus* Walczuchella, *Candidatus* Westeberhardia,
  *Candidatus* Nardonella, *Candidatus* Kinetoplastibacterium.

**S3's residual heterogeneity is reported, not engineered away.** Reassigning a genome from S3 to S1 or
S2 because it is small would violate section 2's rule that strata are ecological and size is never a
criterion. So *T. pallidum* and the Sodalis-like symbiont stay in S3 as their lineage ecology places
them, and the per-stratum size ranges are published with the result so a reader can see the overlap.

## 8. Deliverables

`select_expansion_genomes.py` → `manifests/expansion_candidates.tsv` and
`manifests/expansion_selected.tsv`; the downloaded genomes under
`~/research/datasets/n2-genomes-expansion/`; caller counts appended to the existing
`results/*_counts.tsv`; and `CORPUS_EXPANSION_RESULTS.md` reporting C3 first.
