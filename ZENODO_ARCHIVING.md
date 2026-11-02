# DOI archiving: GitHub release and complete intermediate dataset

**Status: DOI pending.** This repository is public and has a clean independent commit history. The release must be archived in Zenodo only after the maintainer has enabled the GitHub repository integration.

## 1. Connect GitHub and Zenodo (maintainer authorisation)

Open https://zenodo.org/ and log in with GitHub, then open https://zenodo.org/account/settings/github/ . Sync the repository list and enable `kgurgen97/bacterial-pseudogene-annotation-benchmark`. GitHub documentation: https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content .

## 2. Publish an archival software release

After Zenodo is enabled for this repository, publish the release `v1.0.0`. A GitHub Actions workflow `.github/workflows/archive-release.yml` can create the release through the Actions tab (Run workflow), or a maintainer can create a release directly at https://github.com/kgurgen97/bacterial-pseudogene-annotation-benchmark/releases/new and use `RELEASE_NOTES.md`.

The workflow can also be triggered by a new commit on `main` with the exact text `[publish-release]` in the commit message. No other push triggers a release. **Do not trigger it until Zenodo integration is enabled.**

Zenodo will ingest the subsequent GitHub release and mint a version DOI.

## 3. Archive the 12 additional large intermediate tables

The public repository omits 12 large TSV intermediate outputs because of the transfer limits, not because the tables are unnecessary. They can be independently deposited as a Zenodo **dataset** using a complete data ZIP, with dataset authors Gurgen Karapetyan and Astghik Pepoyan and the project title. The archive file, manifest, SHA-256 checksums and draft Zenodo metadata have been prepared in the associated delivery package. The 12 files total 165,315,330 uncompressed bytes. See `LARGE_INTERMEDIATE_DATA.md`.

Keep the **software DOI** and **full data DOI** distinct and insert both into the article Data Availability/Code Availability paragraphs after they have actually been minted.

## 4. Human approvals and scientific limits

Coauthors must confirm authorship, grant funding, rights, conflict-of-interest declarations and the manuscript. Genome-wide reruns additionally need publicly accessioned RefSeq source genomes and external tools. Do not assert that the repository alone fully reproduces the original biological pipeline or that software release is a journal acceptance.

Zenodo guidance: https://help.zenodo.org/docs/github/enable-repository/ and https://help.zenodo.org/docs/github/archive-software/github-upload/
