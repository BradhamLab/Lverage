# ORFfinder benchmark results

Completed comparison of the three assigned genomic RefSeqGene FASTA records using raw `orffinder`, Lverage's corrected adapter, and NCBI standalone ORFfinder. ORF scanning on genomic sequence does not perform splicing or establish the protein annotated for the named gene.

## Results per assigned input

| Input file | RefSeqGene accession.version | Raw package ORFs | Corrected adapter ORFs | NCBI ORFs | Raw/adapter exact-region translation differences | Adapter/NCBI exact matches | Adapter/NCBI positive overlaps | Same longest interval? |
|---|---|---:|---:|---:|---:|---:|---:|---|
| HoxA13.fa | NG_008181.2 | 90 | 90 | 84 | 90 | 84 | 155 | yes |
| Jun.fa | NG_047027.2 | 62 | 62 | 55 | 62 | 55 | 92 | yes |
| MITF.fa | NG_011631.1 | 1761 | 1761 | 1645 | 1761 | 1645 | 2750 | yes |

## Findings

- Across the three assigned inputs, the corrected adapter returned 1,913 candidates and NCBI returned 1,784; there were 1,784 exact interval-and-strand matches.
- Every NCBI candidate had an exact adapter match. The adapter returned 129 additional candidates.
- All three inputs share at least one exact longest-ORF interval and strand between the adapter and NCBI.
- Raw-package/adapter protein-string differences include terminal stop-symbol removal. Their counts do not indicate incorrect proteins.
- These findings apply only to the assigned genomic sequences and recorded settings. ORF scanning does not perform splicing or establish the annotated protein for the named gene.

## Interpretation

- Called-region differences are represented in `candidates.tsv`, `overlaps.tsv`, `matches.tsv`, and `unmatched.tsv`; exact matches and non-exact overlaps are separate.
- Overlap means a positive intersection on the same strand; `overlaps.tsv` includes exact interval pairs with IoU 1.0 and marks them explicitly. Exact matches are also listed separately in `matches.tsv`. The additional one-to-one summary uses deterministic greedy pairing at IoU ≥ 0.5 after exact matches are reserved.
- Longest ORFs are ranked by normalized nucleotide span. All tied longest locations and proteins are retained in `comparisons.tsv`; same-longest means the tied sets share at least one exact interval-and-strand location.
- Raw package translation differences are counted only for exact same-region raw/adapter candidate pairs in `translation_differences.tsv`; candidates found by only one discovery path remain region differences. The corrected adapter candidate reconstruction is checked against production `OrffinderOrfSearcher.get_orfs()` as an ordered protein list and `Counter`, preserving duplicate counts.
- Tool settings and exact commands are in `metadata.json`. NCBI outfmt 1 locations are 1-based inclusive; they are normalized to half-open intervals while retaining stop codons. Outfmt 0 supplies proteins and is cross-checked against the matching CDS record.
- `metadata.json` records each input filename, original FASTA header, accession version, sequence length, and file/sequence SHA-256 checksums. Synthetic controls are validation fixtures and are not part of these assigned-input results.
