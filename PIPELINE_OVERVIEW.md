# Lverage Pipeline Overview

This document explains what **Lverage** does from start to finish.

Lverage is a motif-discovery pipeline for transcription factors. It starts with gene DNA sequences, translates them into candidate proteins, finds DNA-binding domains (DBDs), searches for orthologous proteins in other species, checks whether the DBD is conserved, and then uses those orthologs to infer candidate binding motifs from a motif database such as **JASPAR**.

## Big Picture

Lverage answers this question:

> Given a gene from my species of interest, can I infer likely transcription factor binding motifs by looking at orthologous proteins in better-annotated species?

The core idea is:

1. Find the protein encoded by the input DNA.
2. Identify which part of that protein is the DNA-binding domain.
3. Search for orthologous proteins in selected species.
4. Keep orthologs whose DBD is similar enough to the original gene's DBD.
5. Use those orthologs to query a motif database and collect candidate motifs.

## Inputs

Lverage expects:

- A folder of FASTA files, usually one file per gene
- One or more ortholog species to search against
- An email address for external services
- A list of valid DBD PFAM accessions
- Optional local BLAST configuration
- Thresholds such as DBD identity and e-value cutoffs

Each FASTA file may contain:

- A single DNA sequence
- Multiple scaffold sequences for the same gene

## End-to-End Pipeline

## 1. Read the input gene sequence

For each FASTA file, Lverage reads all DNA records and treats them as belonging to one gene.

Why this matters:

- Some genes may be split across scaffolds
- The pipeline tries multiple candidate translations before committing to one protein

## 2. Search for candidate open reading frames

Lverage uses `orffinder` to translate the gene DNA into possible protein products.

What happens:

- Every input scaffold is translated into possible ORFs
- ORFs are deduplicated
- ORFs are sorted by length
- Only the top ORFs are examined further

Why this matters:

- The longest ORF is not always the one with the biologically relevant DBD
- Lverage looks for an ORF that contains a valid DNA-binding domain

## 3. Scan ORFs for DNA-binding domains

Each candidate ORF is sent to the Pfam/EMBL domain scanning step.

What happens:

- Lverage queries the DBD scanner
- It collects domain hits with:
  - DBD name
  - Start position
  - End position
  - PFAM accession
- It compares those PFAM accessions against the user-provided valid DBD list

Outcome:

- If no valid DBD is found, Lverage moves to the next ORF
- If one or more valid DBDs are found, that ORF becomes the working protein sequence

This is the point where Lverage decides:

> "This looks like a transcription factor protein worth following."

## 4. Search for orthologous proteins

Once a protein sequence is selected, Lverage runs `blastp` to find similar proteins in the requested ortholog species.

This can be done through:

- Remote NCBI BLAST
- Local BLAST, if configured

What Lverage keeps from each BLAST hit:

- Description
- Species
- Accession
- E-value
- Percent identity
- Query coverage
- Protein sequence

Lverage also removes hits that are likely poor annotations, such as records labeled with terms like:

- `hypothetical`
- `unnamed`
- `uncharacterized`
- `partial`
- `isoform`

## 5. Remove orthologs that cannot be mapped cleanly

Lverage checks whether the ortholog species can be mapped to an NCBI taxonomic identifier.

Why this matters:

- Later motif-database queries need a taxonomic context
- Species labels in BLAST output are not always clean enough to use directly

Only orthologs with valid species mapping are allowed to continue.

## 6. Scan ortholog proteins for matching DBDs

Each ortholog protein is scanned for domains again.

Lverage then tries to match ortholog domains to the gene's DBDs by **PFAM accession**.

This is important:

- It does not just ask whether the ortholog has any domain
- It asks whether the ortholog has the **same DBD family** as the gene of interest

For each ortholog, Lverage builds an ordered list of domains aligned to the gene's DBD list:

- Matching DBD -> stored as a DBD object
- Missing DBD -> stored as `None`

## 7. Compare DBD conservation

Now Lverage measures whether the gene DBD and ortholog DBD are similar enough to trust motif transfer.

What happens:

- The gene DBD amino-acid segment is extracted
- The ortholog DBD amino-acid segment is extracted
- The two are globally aligned
- A similarity score is calculated

This score is compared to the configured identity threshold, for example `0.7`.

Meaning:

- If the DBD is conserved enough, the ortholog is considered a valid motif source
- If not, the ortholog is skipped for that DBD

Conceptually, this is the key biological filter in the pipeline:

> "Only use orthologs whose DNA-binding machinery still looks close enough to the original gene."

## 8. Query the motif database

For orthologs that pass the DBD conservation threshold, Lverage asks the motif database for candidate motifs.

In the current code, this is done through **JASPAR**.

What gets sent:

- The ortholog protein sequence
- The ortholog species taxonomic ID
- The main gene protein sequence as context
- The gene DBD
- The matching ortholog DBD

If the ortholog sequence is very long:

- Lverage trims the ortholog protein to a window no longer than 2000 aa
- That window is centered around the **ortholog's matching DBD**
- The full DBD is guaranteed to remain inside the window

Why the ortholog is used instead of the gene:

- The ortholog is the annotated species-specific evidence source
- JASPAR inference is meant to leverage that known ortholog context

## 9. Collect motif hits

Each successful motif hit becomes a `LverageRecord`.

That record combines:

- Gene DBD information
- Ortholog BLAST metadata
- Ortholog-to-gene DBD similarity
- Motif database output

The goal is not just to say:

> "Here is a motif."

It is also to say:

> "Here is which ortholog supported it, how similar that ortholog was, and which DBD relationship justified the transfer."

## 10. Write output

The final output is a tab-separated table.

Each row represents:

- One gene
- One supporting ortholog
- One inferred motif

Typical fields include:

- Gene DBD name
- Gene DBD accession
- Ortholog description
- Ortholog species
- Ortholog e-value
- Ortholog percent identity
- Ortholog query coverage
- Gene-vs-ortholog DBD similarity
- Motif ID
- Motif name
- Motif matrix/PFM
- Motif class
- Motif logo link

## Decision Points and Failure Modes

Lverage can stop early for a gene at several points:

## No usable ORF

- No translated ORF contains a valid DBD

## No usable orthologs

- BLAST finds nothing acceptable
- Orthologs cannot be mapped to species cleanly

## No matching ortholog DBDs

- Ortholog proteins do not contain the same DBD family as the query gene

## DBD not conserved enough

- Ortholog DBD exists, but similarity is below threshold

## No motif database hits

- Ortholog passed previous filters, but no motif could be inferred

This means an empty result does **not** always mean the gene is not a transcription factor. It may simply mean the evidence chain was not strong enough for motif transfer.

## Why the Pipeline Is Structured This Way

Lverage is designed to be conservative.

It does not transfer motifs based only on whole-protein similarity. Instead, it focuses on the part of the protein most responsible for DNA recognition: the DNA-binding domain.

That design has three major advantages:

- It is more biologically meaningful for transcription factors
- It reduces false transfers from unrelated protein regions
- It gives a traceable explanation for every reported motif

## In One Sentence

Lverage takes a gene sequence, finds a likely transcription factor protein, verifies its DNA-binding domains, finds orthologs with conserved DBDs, and uses those orthologs to infer likely DNA-binding motifs.

