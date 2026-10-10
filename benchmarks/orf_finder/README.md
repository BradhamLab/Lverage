# Assigned-input ORFfinder benchmark

This issue-scoped benchmark compares three ORF-finding paths without changing the production adapter or pipeline:

1. **Raw Python package** — `orffinder==1.8` `getORFProteins(..., return_loci=True)` at its direct minimum-length threshold. Native loci are retained, and protein strings are preserved verbatim, including terminal `*` symbols.
2. **Corrected adapter** — benchmark-side reconstruction of `OrffinderOrfSearcher` locus discovery, length filtering, translation, and ordering. Proteins are checked against the production adapter as an ordered list; candidate coordinates come from locus dictionaries.
3. **NCBI standalone** — ORFfinder 0.4.3 CDS FASTA (outfmt 1) and protein FASTA (outfmt 0), joined by ORF identifiers and cross-checked against each other and the original input interval.

## Assigned inputs and scope

The only benchmark dataset is the three professor-provided files under `data/test_inputs/`:

| File | RefSeqGene record |
|---|---|
| `HoxA13.fa` | `NG_008181.2` |
| `Jun.fa` | `NG_047027.2` |
| `MITF.fa` | `NG_011631.1` |

Each file must contain exactly one versioned RefSeqGene FASTA record with a non-empty DNA sequence. The runner checks the accession in the original header, record identifier, DNA alphabet, uniqueness, and exact filename-to-accession assignment. `data/assigned_input_manifest.tsv` records each original header, accession version, sequence length, SHA-256 of the untouched input file, and SHA-256 of the upper-case sequence; the runner verifies the supplied files against this manifest. A run repeats this provenance in `results/metadata.json` and records a checksum for the combined FASTA passed to NCBI. Inputs are consumed locally; no replacement sequences or GenBank annotations are fetched or required.

These are **genomic RefSeqGene sequences**, not spliced mRNA sequences. ORF scanning does not perform exon joining or splicing and cannot establish the mature transcript, the annotated protein, or whether an ORF is the named gene's protein. Results describe the scan of the supplied genomic strings under the recorded settings only.

`data/controls.fasta` contains synthetic **validation fixtures** for strand coordinates, stop handling, nested filtering, and terminal incomplete ORFs. The controls are run and reported separately; they are never included in assigned-input candidate counts or comparisons.

## Settings and comparison semantics

All three tools use ATG starts, standard genetic code (NCBI `-g 1`), minimum 75 nt, no nested filtering, and both strands (`-strand both`). NCBI has a 30-nt minimum floor; the benchmark threshold is above it. The Python adapter requests loci at 74 nt then enforces the actual 75-nt span because the package reports incomplete terminal loci one nucleotide short. The raw-package path is not threshold-corrected.

Candidate coordinates use zero-based, half-open bounds on the original input with strand stored separately. NCBI CDS output includes the terminal stop in the interval while protein output excludes it; the parser checks the outputs pairwise, verifies translation, and checks the CDS sequence against the corresponding input interval. The adapter retains stop codons in genomic intervals but excludes them from protein strings. The adapter's incomplete-terminal translation uses complete codons only.

Exact interval-and-strand matches are paired one-to-one by stable candidate ID. `overlaps.tsv` retains every positive same-strand overlap, including exact matches, and records IoU. The additional one-to-one summary reserves exact matches, then greedily pairs non-exact edges with IoU ≥ 0.5 in deterministic order; it is descriptive, not an optimal assignment. Longest ORFs are ranked by normalized nucleotide span and all tied locations and proteins are retained.

## Findings

For these three assigned genomic RefSeqGene inputs at the recorded settings, the corrected adapter returned 1,913 candidates and NCBI returned 1,784. There were 1,784 exact interval-and-strand matches in total:

| Input | Adapter candidates | NCBI candidates | Exact matches |
|---|---:|---:|---:|
| HoxA13 (`NG_008181.2`) | 90 | 84 | 84 |
| Jun (`NG_047027.2`) | 62 | 55 | 55 |
| MITF (`NG_011631.1`) | 1,761 | 1,645 | 1,645 |
| **Total** | **1,913** | **1,784** | **1,784** |

Every NCBI candidate had an exact adapter match; the adapter returned 129 additional candidates. Each of the three inputs had at least one exact longest-ORF interval-and-strand match between the adapter and NCBI. Raw-package/adapter protein-string differences include terminal stop-symbol removal; their counts are not counts of incorrect proteins. These conclusions apply only to these supplied genomic sequences and the recorded settings. ORF scanning does not perform splicing or establish the annotated protein for any named gene.

## Results and outputs

On completion, the runner writes:

- `results/candidates.tsv`: candidate-level coordinates and proteins from the three tools.
- `results/overlaps.tsv`: all positive same-strand overlap pairs and IoUs.
- `results/matches.tsv`: exact and deterministic IoU-threshold pairs.
- `results/unmatched.tsv`: candidates left unmatched by those summaries.
- `results/translation_differences.tsv`: raw-package/adapter protein comparisons only for exact same-region pairs.
- `results/comparisons.tsv`: per-input counts, tied longest calls, overlaps, matches, and protein multiplicities.
- `results/metadata.json`: tool versions, settings, commands, input provenance/checksums, and separate synthetic-control validation metadata.
- `results/raw/`: NCBI outputs, combined input FASTA, and raw Python package loci for the assigned inputs and controls.
- `results/report.md`: generated summary of completed assigned-input results and their genomic-sequence limitations.

## Setup and execution

The SCC Miniconda module, existing Conda environment, and installed standalone executable are required. From the repository root, run focused tests on the login node:

```bash
module load miniconda/25.3.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /projectnb/paxlab/thomas/envs/lverage-v2
python -m unittest discover -s benchmarks/orf_finder/tests -v
bash -n benchmarks/orf_finder/run.qsub
git diff --check
```

Run the complete assigned-input benchmark on a compute node using `run.qsub`. SCC SGE requires the confirmed account project ID via `-P`; do not infer it from a Unix group or directory name:

```bash
qsub -P "$SCC_PROJECT" benchmarks/orf_finder/run.qsub
```

Do not submit until `SCC_PROJECT` has been confirmed with the SCC account/project tools or the SCC administrator. The job initializes the configured Conda environment, verifies the executable, and runs only the three assigned inputs plus separate synthetic validation fixtures.
