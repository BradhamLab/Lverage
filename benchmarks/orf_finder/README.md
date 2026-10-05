# Human transcript ORFfinder benchmark

This is an issue-scoped benchmark only; it does not change the production adapter or pipeline. It compares three outputs separately:

1. **Raw Python package** — `orffinder==1.8` `getORFProteins(..., return_loci=True)` at its direct minimum-length threshold. Candidate loci come from package coordinates; protein strings are retained verbatim, including any `*`.
2. **Corrected adapter** — benchmark-only reproduction of `OrffinderOrfSearcher` locus discovery, length filtering, terminal-stop removal, complete-codon translation, and candidate ordering. Candidate proteins are verified against the actual production adapter as an ordered list and `Counter`, so duplicates matter. Candidate coordinates come from the locus dictionaries, never by matching protein strings or list positions.
3. **NCBI standalone** — ORFfinder 0.4.3 CDS FASTA (outfmt 1) and protein FASTA (outfmt 0), joined through ORF identifiers and verified against each other and the original input.

## Dataset and provenance

The agreed data set is ten version-pinned human RefSeq mRNA records:

| Accession | Gene |
|---|---|
| NM_000546.6 | TP53 |
| NM_002467.6 | MYC |
| NM_003106.4 | SOX2 |
| NM_000280.5 | PAX6 |
| NM_000125.4 | ESR1 |
| NM_004496.4 | FOXA1 |
| NM_002049.4 | GATA1 |
| NM_000457.5 | HNF4A |
| NM_001754.5 | RUNX1 |
| NM_152739.4 | HOXA9 |

Before saving, the fetch command requires exact returned accession versions and exact CDS gene qualifiers; a mismatch is an error, not a reason to substitute a record. It fetches once in one Entrez request and saves the raw GenBank response to `data/human_refseq.gb`. The exact sequences used for both tools are in `data/human_refseq.fasta`. `data/dataset.tsv` records sequence SHA-256, accession/version, gene, CDS bounds/strand/length, annotated CDS translation and its checksum, the complete translation-related CDS qualifier map, record definition, and fetch timestamp. The source is NCBI RefSeq nuccore GenBank. `data/controls.fasta` contains synthetic controls separately; these never enter the human results.

The annotated CDS is reference evidence, not a binary correctness label. Its translation and qualifiers are retained; unsupported compound or multiple principal CDS representations stop the run rather than being silently simplified. One annotation caveat was observed for MYC (NM_002467.6): its CDS starts with CTG and the GenBank note specifies non-AUG initiation; the annotated translation begins with methionine while direct standard-code translation begins with leucine. The annotated translation is preserved verbatim and the difference is recorded in `dataset.tsv`; this does not change the benchmark's ATG-only discovery setting.

## Matched settings and semantics

Both tools use ATG starts only, standard genetic code (NCBI `-g 1`), minimum 75 nt, no nested filtering, and both strands (NCBI `-strand both`). NCBI also has a hard 30-nt minimum floor; the 75-nt run is above it. The Python adapter asks the package for loci at 74 nt then enforces the actual span at 75 nt because `orffinder==1.8` reports incomplete terminal loci one nucleotide short. The raw-package series is deliberately not threshold-corrected and reports its direct package result.

For complete ORFs, NCBI's CDS output interval includes the terminal stop, and NCBI's protein output excludes it. The normalized candidate interval retains all CDS nucleotides including the stop. The corrected adapter also retains the stop in its genomic interval but excludes it from translated protein. Python's incomplete-terminal behavior is separately captured: the corrected adapter counts the actual sequence span and trims only residual one/two bases from translation; the raw package output remains raw. Synthetic controls verify plus/minus coordinates, complete stop handling, nested filtering, and incomplete terminal calls at both 30 and 75 nt. NCBI output is parsed only after its output formats are checked; both FASTA outputs are preserved and cross-validated.

## Coordinates, overlap, and pairing

Every candidate interval is zero-based, half-open on the original input, with `+` or `-` stored separately. Candidate identity is tool-specific and stable; it is derived from source record and native locus/ORF identifier, not a protein sequence.

Exact interval-and-strand matches are paired separately, one-to-one, in ascending stable candidate-ID order. `overlaps.tsv` contains **every** same-strand pair with positive intersection, including exact pairs, and its intersection-over-union (IoU); exact pairs are flagged and also listed in `matches.tsv`. The extra IoU-threshold summary is deterministic greedy pairing: reserve exact matches first; sort remaining non-exact edges with IoU ≥ 0.5 by descending IoU, then ascending left candidate ID, then ascending right candidate ID; accept an edge only if neither candidate is already paired. `unmatched.tsv` lists candidates not paired in either exact or threshold summaries. This greedy summary is descriptive, not an optimal assignment claim; the full overlap table retains all positive overlaps and their IoUs.

The longest candidates are ranked by normalized nucleotide span. All tied locations and translated proteins are retained. Raw package translation differences are explicitly distinct from region-discovery differences. Additional ORFs are not presumed erroneous because they are absent from the annotated CDS.

## Findings

Across these ten version-pinned human RefSeq transcripts, the corrected adapter returned 250 candidates and NCBI returned 235. There were 229 exact interval-and-strand matches, alongside additional non-exact same-strand overlaps. Nine transcripts shared an exact longest-ORF interval.

HOXA9 was the exception: the adapter's longest interval was [0,905) on the reverse strand, versus [2,905) for NCBI. Their longest-ORF proteins were identical; the adapter interval includes two untranslated trailing bases. Raw-package versus adapter protein-string differences include terminal stop-symbol removal, so they should not all be characterized as truncation. These results do not establish that either tool is universally more accurate; findings apply only to this dataset and the recorded settings.

## Files and outputs

- `benchmark.py`: dataset validation/fetching, controlled checks, benchmark runner and report writer.
- `tests/test_benchmark.py`: parser, coordinate, overlap/pairing, adapter-association tests.
- `data/`: saved inputs/provenance; the synthetic controls are intentionally separate.
- `results/raw/`: both NCBI human formats and raw synthetic outputs.
- `results/raw/python_package_loci.json`: untouched package-native locus fields and protein strings for human records; control package loci are separately preserved as `controls_python_package_loci.json`.
- `results/candidates.tsv`: one candidate-level row per tool output, including coordinates and proteins.
- `results/overlaps.tsv`: all same-strand positive-overlap pairs, including exact pairs, and IoUs; exact pairs are marked and separately listed in `matches.tsv`.
- `results/matches.tsv`: exact and additional one-to-one IoU-threshold pairs.
- `results/translation_differences.tsv`: protein comparison for exact same-region raw-package/adapter pairs, kept distinct from candidate discovery differences.
- `results/unmatched.tsv`: unpaired candidate IDs by tool.
- `results/comparisons.tsv`: per-transcript counts, longest calls/ties, match counts and protein multiplicities.
- `results/metadata.json`: versions, settings, exact commands, checksums, control results, and normalization/matching semantics.
- `results/report.md`: concise generated answer to region overlap and largest-ORF questions.

## Setup and execution

Use the SCC Miniconda module, the existing Conda environment, and the installed executable. From the repository root:

```bash
module load miniconda/25.3.1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /projectnb/paxlab/thomas/envs/lverage-v2
python -m benchmarks.orf_finder.benchmark fetch --email YOUR_CONTACT_EMAIL
python -m unittest discover -s benchmarks/orf_finder/tests -v
```

Run the focused unit tests on the login node after setup:

```bash
python -m unittest discover -s benchmarks/orf_finder/tests -v
```

The fetch step is run once; the benchmark run consumes only saved files. Run the complete human benchmark on a compute node through `run.qsub`, which initializes the same explicit Conda environment and sets the NCBI executable path. SCC SGE project directives use `-P`; this checkout does not verify the account's scheduler project ID. Provide it explicitly at submission (do not infer it from the Unix group or directory):

```bash
qsub -P "$SCC_PROJECT" benchmarks/orf_finder/run.qsub
```

Do not submit until `SCC_PROJECT` has been confirmed with the SCC account/project tools or the SCC administrator.
