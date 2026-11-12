"""Local JASPAR profile inference adapted for Lverage's motif adapter.

This module adapts the inference functions from ``infer_profile.py`` in
https://github.com/amathelier/JASPAR-inference-tool, source revision
ac4024148f1f140c7370d7b13a70556b4d13f7a3 (version 2026.1.0). The original
copyright notice and MIT license terms are preserved in
``jaspar_inference_license.txt``. Upstream attribution names Anthony Mathelier,
Oriol Fornes, Xi Zhang, and the JASPAR Consortium.

The inference calculations and decisions below retain the upstream procedure:
HMMER Pfam domain detection/alignment; BLASTP against each configured taxon DB;
Rost-curve identity filtering; lower-case insertion removal; the best positional
identity/BLOSUM62 alignment; Cis-BP HSim scoring; and one matched DBD's JASPAR
profiles per accepted BLAST target. External programs are called directly with
argument vectors, not through a shell.
"""

from pathlib import Path
import json
import math
import re
import subprocess

from Bio import SearchIO, SeqIO
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np

from .jaspar_models import load_cisbp_models, score_alignment_result


BLOSUM62 = substitution_matrices.load("BLOSUM62")
JASPAR_TAXA = ("fungi", "insects", "nematodes", "plants", "urochordates", "vertebrates")
BLAST_OUTFMT = "sseqid pident length qstart qend sstart send evalue bitscore ppos qlen slen"


def load_jaspar_data(data_dir, taxons=JASPAR_TAXA):
    """Load the fork's per-taxon Pfam/profile/UniProt mappings."""
    data_dir = Path(data_dir)
    jaspar = {}
    pfams = {}
    profiles = {}
    uniprots = {}

    for taxon in taxons:
        for suffix, destination in (
            ("pfam", pfams), ("profiles", profiles), ("uniprot", uniprots),
        ):
            path = data_dir / f"{taxon}.{suffix}.json"
            with path.open(encoding="utf-8") as handle:
                for key, values in json.load(handle).items():
                    destination.setdefault(key, values)

    for uniprot, domains in pfams.items():
        entry = jaspar.setdefault(uniprot, {})
        entry.setdefault("pfam", domains)
        entry.setdefault("profiles", [])
        for profile in uniprots[uniprot][0]:
            entry["profiles"].append([profile, profiles[profile]])
    return jaspar


def _checked_run(arguments, *, stage, stdout=None):
    try:
        return subprocess.run(
            arguments,
            check=True,
            stdout=subprocess.PIPE if stdout is None else stdout,
            stderr=subprocess.PIPE,
            text=True,
        )
    except FileNotFoundError as error:
        raise RuntimeError(f"{stage} executable was not found: {arguments[0]}") from error
    except subprocess.CalledProcessError as error:
        diagnostic = (error.stderr or error.stdout or "").strip()
        if not diagnostic:
            diagnostic = "no diagnostic output"
        raise RuntimeError(
            f"{stage} failed with status {error.returncode}: {diagnostic[-4000:]}"
        ) from error


def _write_fasta(record, path):
    with Path(path).open("w", encoding="utf-8") as handle:
        SeqIO.write(record, handle, "fasta")


def read_hmmscan_domains(path):
    """Apply upstream model (1e-5) and conditional-domain (0.01) cutoffs."""
    domains = []
    for result in SearchIO.parse(str(path), "hmmscan3-domtab"):
        for model in result.iterhits():
            if model.evalue > 1e-5:
                continue
            for domain in model.hsps:
                if domain.evalue_cond > 0.01:
                    continue
                domains.append((model.id, domain.query_start, domain.query_end, domain.evalue_cond))
    return domains


def non_overlapping_domains(domains):
    """Retain the lowest conditional-E-value domains with upstream overlap rule."""
    retained = []
    for domain in sorted(domains, key=lambda item: item[-1]):
        if any(domain[1] < previous[2] and domain[2] > previous[1]
               for previous in retained):
            continue
        retained.append(domain)
    return retained


def scan_and_align_domains(record, data_dir, work_dir, *, hmmscan_path, hmmalign_path):
    """Run upstream Pfam scan, thresholding, overlap selection, and HMM alignment."""
    data_dir = Path(data_dir)
    work_dir = Path(work_dir)
    sequence_path = work_dir / "query.fasta"
    _write_fasta(record, sequence_path)
    domtbl_path = work_dir / "query.domtbl"
    _checked_run(
        [hmmscan_path, "--domtblout", str(domtbl_path),
         str(data_dir / "pfam" / "All.hmm"), str(sequence_path)],
        stage="hmmscan",
        stdout=subprocess.DEVNULL,
    )
    domains = non_overlapping_domains(read_hmmscan_domains(domtbl_path))
    alignments = []
    domain_path = work_dir / "domain.fasta"
    for pfam_id, start, end, evalue in sorted(domains, key=lambda item: item[1]):
        domain_record = SeqRecord(
            record.seq[start:end], id=record.id, name=record.name,
            description=record.description,
        )
        _write_fasta(domain_record, domain_path)
        process = _checked_run(
            [hmmalign_path, "--outformat", "PSIBLAST",
             str(data_dir / "pfam" / f"{pfam_id}.hmm"), str(domain_path)],
            stage="hmmalign",
        )
        aligned = read_psiblast_alignment(process.stdout)
        alignments.append((pfam_id, aligned, start + 1, end, evalue))
    return alignments


def read_psiblast_alignment(text):
    """Concatenate non-whitespace trailing alignment blocks as upstream does."""
    aligned = ""
    for chunk in text.split("\n"):
        match = re.search(r"\s+(\S+)$", chunk)
        if match:
            aligned += match.group(1)
    return aligned


def _rost_identity_cutoff(alignment_length, n=5):
    """Upstream Rost pairwise sequence-identity curve."""
    return n + 480 * pow(
        alignment_length,
        -0.32 * (1 + math.exp(float(-alignment_length) / 1000)),
    )


def _is_over_rost_curve(percent_identity, alignment_length, n=5):
    return percent_identity >= _rost_identity_cutoff(alignment_length, n)


def blast_search(record, data_dir, query_fasta, *, blastp_path, taxons=JASPAR_TAXA, n=5):
    """Run upstream BLASTP outfmt and Rost identity/coverage filtering."""
    data_dir = Path(data_dir)
    blast_results = set()
    for taxon in taxons:
        process = _checked_run(
            [blastp_path, "-db", str(data_dir / f"{taxon}.fa"),
             "-outfmt", f"6 {BLAST_OUTFMT}", "-query", str(query_fasta)],
            stage=f"blastp ({taxon})",
        )
        for line in process.stdout.splitlines():
            fields = line.split("\t")
            if len(fields) != 12:
                continue
            target_id = fields[0]
            percent_identity = float(fields[1])
            alignment_length = int(fields[2])
            query_start_end = f"{fields[3]}-{fields[4]}"
            target_start_end = f"{fields[5]}-{fields[6]}"
            evalue = float(fields[7])
            bit_score = float(fields[8])
            percent_similarity = float(fields[9])
            query_coverage = (int(fields[4]) - int(fields[3]) + 1) * 100 / float(fields[10])
            target_coverage = (int(fields[6]) - int(fields[5]) + 1) * 100 / float(fields[11])
            joint_coverage = math.sqrt(query_coverage * target_coverage)
            blast_results.add((
                record.id, target_id, query_start_end, target_start_end,
                evalue, bit_score, percent_identity, alignment_length,
                percent_similarity, joint_coverage,
            ))

    ordered = sorted(blast_results, key=lambda item: item[-1], reverse=True)
    return [result for result in ordered
            if _is_over_rost_curve(result[6], result[7], n)]


def _target_pfam_alignments(blast_results, jaspar):
    pfam_alignments = {}
    for result in blast_results:
        target_id = result[1]
        pfam_alignments.setdefault(target_id, {})
        for alignment in jaspar[target_id]["pfam"]:
            pfam_alignments[target_id].setdefault(alignment[0], []).append(alignment[1])
    return pfam_alignments


def _remove_insertions(alignment):
    return alignment.translate(str.maketrans("", "", "abcdefghijklmnopqrstuvwxyz"))


def _score_residue(first, second, similarity):
    if similarity == "identity":
        return 1 if first == second else 0
    if similarity == "blosum62":
        if first == second == "-":
            return 1
        if first == "-" or second == "-":
            return -4
        try:
            return int(BLOSUM62[first][second])
        except (KeyError, IndexError):
            return 0
    return 0


def positional_alignment_score(first_alignments, second_alignments, similarity="identity"):
    """Upstream __get_X: maximize summed positional score over contiguous DBD alignments."""
    if not first_alignments or not second_alignments:
        raise ValueError("at least one aligned domain is required for similarity regression")
    if len(first_alignments) < len(second_alignments):
        first_alignments, second_alignments = second_alignments, first_alignments

    scores = []
    for offset in range(len(first_alignments) - len(second_alignments) + 1):
        row = [0] * len(first_alignments[offset])
        for second_offset, second_alignment in enumerate(second_alignments):
            first_alignment = first_alignments[offset + second_offset]
            if len(first_alignment) != len(second_alignment):
                raise ValueError("upstream HMM alignment rows have unequal lengths")
            for position, (first, second) in enumerate(zip(first_alignment, second_alignment)):
                row[position] += _score_residue(first, second, similarity)
        scores.append(np.asarray(row) / len(second_alignments))
    scores.sort(key=lambda values: sum(values), reverse=True)
    return scores[0]


def infer_sequence(sequence, *, data_dir, cisbp_models, jaspar_data, work_dir,
                   blastp_path, hmmscan_path, hmmalign_path, query_id="lverage_query",
                   n=5, taxons=JASPAR_TAXA):
    """Return upstream inference rows for one protein sequence, without CLI/TSV wrappers."""
    record = SeqRecord(Seq(sequence), id=query_id, name=query_id, description=query_id)
    query_alignments = scan_and_align_domains(
        record, data_dir, work_dir,
        hmmscan_path=hmmscan_path, hmmalign_path=hmmalign_path,
    )
    if not query_alignments:
        return []

    query_by_domain = {}
    for pfam_id, aligned, _start, _end, _evalue in query_alignments:
        query_by_domain.setdefault(pfam_id, []).append(aligned)

    blast_results = blast_search(
        record, data_dir, Path(work_dir) / "query.fasta",
        blastp_path=blastp_path, taxons=taxons, n=n,
    )
    target_by_protein = _target_pfam_alignments(blast_results, jaspar_data)
    predictions = []
    for result in blast_results:
        target_id = result[1]
        # Upstream directly indexes this result: a BLAST target absent from
        # reference data is an execution/data mismatch, not a biological no-hit.
        target_domains = target_by_protein[target_id]
        for pfam_id, query_domain_alignments in query_by_domain.items():
            if pfam_id not in target_domains:
                continue
            query_rows = [_remove_insertions(row) for row in query_domain_alignments]
            target_rows = [_remove_insertions(row) for row in target_domains[pfam_id]]
            percent_identity_by_position = positional_alignment_score(
                list(query_rows), list(target_rows), "identity",
            )
            average_blosum62_by_position = positional_alignment_score(
                list(query_rows), list(target_rows), "blosum62",
            )
            percent_identity = sum(percent_identity_by_position) / len(percent_identity_by_position)
            alignment_features = {
                "ByPos.PctID": percent_identity_by_position,
                "ByPos.AvgB62": average_blosum62_by_position,
                "PctID_L": percent_identity,
            }
            model = cisbp_models.get(pfam_id, cisbp_models[None])
            _score, classification = score_alignment_result(alignment_features, model)
            if classification != "HSim":
                continue

            for matrix_id, tf_name in jaspar_data[target_id]["profiles"]:
                predictions.append({
                    "query": result[0],
                    "tf_name": tf_name,
                    "matrix_id": matrix_id,
                    "evalue": result[4],
                    "query_region": result[2],
                    "tf_region": result[3],
                    "dbd_identity": round(percent_identity, 3),
                })
            # Upstream stops after the first matching/accepted query DBD per target.
            break

    predictions.sort(key=lambda row: (
        row["evalue"], row["tf_name"], float(row["matrix_id"][2:]),
    ))
    return predictions
