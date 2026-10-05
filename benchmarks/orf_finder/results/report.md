# ORFfinder benchmark results

This report is generated from the saved RefSeq inputs and preserved tool outputs. It distinguishes raw `orffinder` package behavior, Lverage's corrected adapter, and NCBI standalone ORFfinder. No additional ORFs are classified as errors merely because they are not annotated CDS.

## Results per transcript

| Accession | Gene | Raw package | Corrected adapter | NCBI | Raw/adapter translation differences | Adapter/NCBI exact | Adapter/NCBI positive overlaps | Same longest interval? |
|---|---|---:|---:|---:|---:|---:|---:|---|
| NM_000546.6 | TP53 | 16 | 16 | 14 | 16 | 14 | 29 | yes |
| NM_002467.6 | MYC | 22 | 22 | 20 | 22 | 20 | 30 | yes |
| NM_003106.4 | SOX2 | 11 | 11 | 11 | 11 | 11 | 15 | yes |
| NM_000280.5 | PAX6 | 17 | 17 | 16 | 17 | 16 | 30 | yes |
| NM_000125.4 | ESR1 | 43 | 43 | 39 | 40 | 38 | 69 | yes |
| NM_004496.4 | FOXA1 | 17 | 17 | 16 | 16 | 15 | 29 | yes |
| NM_002049.4 | GATA1 | 11 | 11 | 11 | 9 | 9 | 23 | yes |
| NM_000457.5 | HNF4A | 54 | 54 | 52 | 52 | 51 | 94 | yes |
| NM_001754.5 | RUNX1 | 43 | 43 | 42 | 43 | 42 | 65 | yes |
| NM_152739.4 | HOXA9 | 16 | 16 | 14 | 13 | 13 | 31 | no |

## Interpretation

- Called-region differences are represented in `candidates.tsv`, `overlaps.tsv`, `matches.tsv`, and `unmatched.tsv`; exact matches and non-exact overlaps are separate.
- Overlap means a positive intersection on the same strand; `overlaps.tsv` includes exact interval pairs with IoU 1.0 and marks them explicitly. Exact matches are also listed separately in `matches.tsv`. The additional one-to-one summary uses deterministic greedy pairing at IoU ≥ 0.5 after exact matches are reserved.
- Longest ORFs are ranked by normalized nucleotide span. All tied longest locations and proteins are retained in `comparisons.tsv`; the table's same-longest flag means at least one exact location-and-strand intersection between the tied sets.
- Raw package translation differences are counted only for exact same-region raw/adapter candidate pairs in `translation_differences.tsv`; candidates found by only one discovery path remain region differences. The corrected adapter candidate reconstruction is checked against production `OrffinderOrfSearcher.get_orfs()` as an ordered protein list and `Counter`, preserving duplicate counts.
- The archived GenBank CDS annotation is a reference. Candidate/CDS interval IoU, strand, exactness, annotated translation, and complete CDS qualifiers are in the candidate and dataset tables; additional ORFs are not automatically errors.
- Tool settings and exact commands are in `metadata.json`. NCBI outfmt 1 locations are 1-based inclusive; they are normalized to half-open intervals while retaining stop codons. Outfmt 0 supplies proteins and is cross-checked against the matching CDS record.
