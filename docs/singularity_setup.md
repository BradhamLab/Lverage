# Singularity Setup Requirements

This document tracks all software, external tools, databases, environment configuration, and setup steps required to run Lverage.

## Python Environment

- Python 3.11
- Biopython
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