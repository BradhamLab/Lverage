# Singularity Setup Requirements

This document tracks all software, external tools, databases, environment configuration, and setup steps required to run Lverage.

## Python Environment

- Python 3.11
- Biopython
- orffinder 1.8 (Python ORF finder)
- Requests

Developer/documentation dependencies:
- pytest
- Sphinx
- numpydoc
- pydata-sphinx-theme

## BLAST

Required executables:
- `blastp`
- `blastdbcmd`
- `makeblastdb`

Required data:
- Local protein BLAST databases

Setup notes:
- BLAST+ must be installed and available on `PATH`.
- Protein databases should be created with `makeblastdb`.
- Use `-parse_seqids` so sequences can be retrieved later with `blastdbcmd`.

## ORFfinder comparison benchmark

- The benchmark compares Python `orffinder==1.8` with NCBI standalone ORFfinder.
- NCBI ORFfinder 0.4.3 executable used on the SCC:
	`/projectnb/paxlab/thomas/tools/orffinder/ORFfinder`.
- The production Lverage adapter continues to use the Python package; the NCBI
	executable is required only for the optional benchmark.
- The benchmark uses only the three professor-provided genomic RefSeqGene FASTA
	files under `benchmarks/orf_finder/data/test_inputs/`; it does not fetch
	replacement sequences or require GenBank records/CDS annotations.
- These are genomic sequences: ORF scanning does not splice exons and cannot
	establish the named gene's mature transcript or annotated protein. Synthetic
	sequences under `data/controls.fasta` are validation fixtures kept separate
	from assigned-input results.
- The runner, parser/coordinate tests, SCC job script, input provenance/checksums,
	and generated outputs are under `benchmarks/orf_finder/`. The assigned-input
	benchmark completed on SCC, and its saved results are under
	`benchmarks/orf_finder/results/`.
- The standalone binary is not installed by Conda or pip. Confirm its version
	and executable path before running. The benchmark job script initializes
	`/projectnb/paxlab/thomas/envs/lverage-v2` explicitly.
- SCC SGE requires an account scheduler project via `qsub -P`; this repository
	does not infer that identifier from a Unix group or directory name. Set the
	confirmed project at submission time as documented in the benchmark README.