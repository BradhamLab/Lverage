# ORFfinder benchmark

This benchmark compares three ORF-finding paths without changing the production adapter or pipeline:

1. **Raw Python package** — `orffinder==1.8` `getORFProteins(..., return_loci=True)` at its direct minimum-length threshold. Native loci are retained, and protein strings are preserved verbatim, including terminal `*` symbols.
2. **Corrected adapter** — benchmark-side reconstruction of `OrffinderOrfSearcher` locus discovery, length filtering, translation, and ordering. Proteins are checked against the production adapter as an ordered list; candidate coordinates come from locus dictionaries.
3. **NCBI standalone** — ORFfinder 0.4.3 CDS FASTA (outfmt 1) and protein FASTA (outfmt 0), joined by ORF identifiers and cross-checked against each other and the original input interval.

## Test inputs and scope

The only benchmark dataset is three test files:

| File | RefSeqGene record |
|---|---|
| `HoxA13.fa` | `NG_008181.2` |
| `Jun.fa` | `NG_047027.2` |
| `MITF.fa` | `NG_011631.1` |

The FASTA files themselves are **not versioned** (the repository holds code, not data). Place a local copy of the three files in a directory and pass it with `--input-dir` or the `LVERAGE_ORF_BENCHMARK_INPUTS` environment variable; `data/test_inputs/` is the git-ignored default location. Each file must contain exactly one versioned RefSeqGene FASTA record with a non-empty DNA sequence. The runner checks the accession in the original header, record identifier, DNA alphabet, uniqueness, and exact filename-to-accession assignment. `data/test_input_manifest.tsv` records each original header, accession version, sequence length, SHA-256 of the untouched input file, and SHA-256 of the upper-case sequence; whenever that manifest is present the runner verifies the files against it, so a run on different bytes fails loudly. A run repeats this provenance in `results/metadata.json` and records a checksum for the combined FASTA passed to NCBI. Inputs are consumed locally; no replacement sequences or GenBank annotations are fetched or required.

This runner is deliberately **not a general ORF benchmark**: it accepts only these three filenames and accessions and refuses anything else. That is so the recorded numbers can only ever be reproduced against the same bytes. Generalizing it to arbitrary inputs would mean dropping the manifest and accession checks, which is a different tool.

These are **genomic RefSeqGene sequences**, not spliced mRNA sequences. ORF scanning does not perform exon joining or splicing and cannot establish the mature transcript, the annotated protein, or whether an ORF is the named gene's protein. Results describe the scan of the genomic strings under the recorded settings only.

`data/controls.fasta` contains synthetic **validation fixtures** for strand coordinates, stop handling, nested filtering, and terminal incomplete ORFs. The controls are run and reported separately; they are never included in test-input candidate counts or comparisons.

## Settings and comparison semantics

The preserved baseline uses ATG starts, standard genetic code (NCBI `-g 1`), no nested filtering, and both strands (`-strand both`). Its raw Python package and NCBI settings request 75 nt; the adapter requests loci at 74 nt, then enforces an actual 75-nt span because the package reports incomplete terminal loci one nucleotide short. The raw-package path is not threshold-corrected. NCBI has a 30-nt minimum floor.

Equal numeric thresholds do not imply equivalent length filtering across tools. In this baseline, NCBI ORFfinder 0.4.3 at `-ml 75` returned 1,784 candidates, while the adapter's 75-nt filter returned 1,913; 129 adapter candidates had no exact NCBI interval-and-strand match. All 129 were complete 75-nt ORFs whose intervals included a terminal stop codon. A separate validation reran NCBI ORFfinder with `-ml 72` and all other baseline options unchanged. Against the adapter at 75 nt, this produced 1,913 exact interval-and-strand matches across these inputs, preserving candidate multiplicities; all paired protein strings were identical, and the longest intervals matched for all three inputs. This follow-up does not change the benchmark runner's baseline settings or default.

Candidate coordinates use zero-based, half-open bounds on the original input with strand stored separately. NCBI CDS output includes the terminal stop in the interval while protein output excludes it; the parser checks the outputs pairwise, verifies translation, and checks the CDS sequence against the corresponding input interval. The adapter retains stop codons in genomic intervals but excludes them from protein strings. The adapter's incomplete-terminal translation uses complete codons only.

Exact interval-and-strand matches are paired one-to-one by stable candidate ID. `overlaps.tsv` retains every positive same-strand overlap, including exact matches, and records IoU. The additional one-to-one summary reserves exact matches, then greedily pairs non-exact edges with IoU ≥ 0.5 in deterministic order; it is descriptive, not an optimal assignment. Longest ORFs are ranked by normalized nucleotide span and all tied locations and proteins are retained.

## Findings

For these three genomic RefSeqGene test inputs, the preserved `-ml 75` baseline comparison returned 1,913 corrected-adapter candidates and 1,784 NCBI candidates, with 1,784 exact interval-and-strand matches:

| Input | Adapter candidates | NCBI candidates | Exact matches |
|---|---:|---:|---:|
| HoxA13 (`NG_008181.2`) | 90 | 84 | 84 |
| Jun (`NG_047027.2`) | 62 | 55 | 55 |
| MITF (`NG_011631.1`) | 1,761 | 1,645 | 1,645 |
| **Total** | **1,913** | **1,784** | **1,784** |

Every baseline NCBI candidate had an exact adapter match; the adapter returned 129 additional candidates, all complete 75-nt ORFs with an interval-terminal stop codon. Separately, the NCBI ORFfinder 0.4.3 `-ml 72` validation matched all 1,913 adapter interval-and-strand occurrences, with identical protein strings and matching longest intervals for each input. The Counter-based comparison preserved multiplicities. Raw-package/adapter protein-string differences include terminal stop-symbol removal; their counts are not counts of incorrect proteins. These conclusions apply only to these three genomic sequences and the stated settings. ORF scanning does not perform splicing or establish the annotated protein for any named gene.

The separate lower-threshold check can be reproduced by running NCBI ORFfinder 0.4.3 on the same combined FASTA used in the baseline run, changing only `-ml` from 75 to 72, and saving both output formats:

```bash
ORFfinder -in /path/to/combined_input.fasta -s 0 -g 1 -ml 72 -n false -strand both -outfmt 0 -out /path/to/results/ncbi_min72_proteins.fasta
ORFfinder -in /path/to/combined_input.fasta -s 0 -g 1 -ml 72 -n false -strand both -outfmt 1 -out /path/to/results/ncbi_min72_cds.fasta
```

These commands describe a follow-up validation only; they do not change the benchmark runner's default `-ml 75` comparison.

## Results and outputs

Generated outputs are written under `results/` (override with `--results-dir`). That directory is git-ignored; the summary of the recorded run is the Findings section above. The runner writes:

- `results/candidates.tsv`: candidate-level coordinates and proteins from the three tools.
- `results/overlaps.tsv`: all positive same-strand overlap pairs and IoUs.
- `results/matches.tsv`: exact and deterministic IoU-threshold pairs.
- `results/unmatched.tsv`: candidates left unmatched by those summaries, one row per pairwise comparison (`comparison` and `unmatched_against` columns say which), so a candidate may legitimately appear up to twice.
- `results/translation_differences.tsv`: raw-package/adapter protein comparisons only for exact same-region pairs.
- `results/comparisons.tsv`: per-input counts, tied longest calls, overlaps, matches, and protein multiplicities.
- `results/metadata.json`: tool versions, settings, commands, input provenance/checksums, and separate synthetic-control validation metadata.
- `results/raw/`: NCBI outputs, combined input FASTA, and raw Python package loci for the test inputs and controls.
- `results/report.md`: generated summary of completed test-input results and their genomic-sequence limitations.

## Setup and execution

Nothing in this directory refers to a particular user's environment, clone, or tool install. You need:

1. A Python 3.11+ environment with this checkout installed, so that the adapter under test is the one that gets imported. From the repository root: `pip install -e .` (installs `orffinder==1.8` and `biopython`). If the environment already has another Lverage clone installed in editable mode, either reinstall from this checkout or prefix commands with `PYTHONPATH=src` as below.
2. The NCBI standalone ORFfinder binary (0.4.3 was used for the recorded run). It is not installed by Conda or pip; download it from NCBI. Tell the benchmark where it is with `--executable`, the `ORFFINDER_EXECUTABLE` environment variable, or by putting an `ORFfinder` command on `PATH`.
3. A local directory containing `HoxA13.fa`, `Jun.fa`, and `MITF.fa` (see above).

Run the focused tests from the repository root:

```bash
PYTHONPATH=src python -m unittest discover -s benchmarks/orf_finder/tests -v
bash -n benchmarks/orf_finder/run.qsub
```

Run the benchmark directly:

```bash
PYTHONPATH=src python -m benchmarks.orf_finder.benchmark run \
    --executable /path/to/ORFfinder \
    --input-dir /path/to/combined_fastas
```

Or on a cluster that uses Grid Engine (SGE/UGE), submit `run.qsub` **from the repository root** of the checkout you want to benchmark. The job script reads its configuration from environment variables (see the header of `run.qsub`); add whatever scheduler options your site requires, such as `-P <project>`:

```bash
qsub \
    -v ORFFINDER_EXECUTABLE=/path/to/ORFfinder \
    -v ORF_BENCHMARK_INPUTS=/path/to/test_fastas \
    -v CONDA_ENV=my-lverage-env \
    benchmarks/orf_finder/run.qsub
```

The job verifies it was submitted from a repository root, verifies the executable, and runs only the three test inputs plus separate synthetic validation fixtures. On a workstation, just run the `python -m` command above instead.
