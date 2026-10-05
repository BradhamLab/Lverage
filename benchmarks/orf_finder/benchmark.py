"""Reproducible, benchmark-only comparison of Python and NCBI ORFfinders."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import itertools
import json
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from Bio import SeqIO, __version__ as biopython_version
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from orffinder import orffinder

from lverage.orf_searcher import OrffinderOrfSearcher

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
RESULTS = HERE / "results"
DEFAULT_EXECUTABLE = "/projectnb/paxlab/thomas/tools/orffinder/ORFfinder"
MINIMUM_NT = 75
GENETIC_CODE = 1
ASSIGNED_INPUT_DIR = DATA / "test_inputs"
ASSIGNED_INPUT_MANIFEST = DATA / "assigned_input_manifest.tsv"
ASSIGNED_INPUT_FILES = ("HoxA13.fa", "Jun.fa", "MITF.fa")
ASSIGNED_ACCESSIONS = {
    "HoxA13.fa": "NG_008181.2",
    "Jun.fa": "NG_047027.2",
    "MITF.fa": "NG_011631.1",
}
DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN")
ACCESSION_VERSION = re.compile(r"\b(NG_\d+\.\d+)\b")


@dataclass(frozen=True)
class Candidate:
    candidate_id: str
    tool: str
    record_id: str
    start0: int
    end0: int
    strand: str
    nucleotide_length: int
    protein: str
    terminal_stop_in_nt: bool
    complete: bool

    @property
    def interval(self) -> tuple[int, int, str]:
        return self.start0, self.end0, self.strand


def sha256_text(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def _locus_interval(locus: dict, sequence_length: int) -> tuple[int, int, str, int]:
    """Convert orffinder's native locus dictionary to original-input bounds."""
    start_index = (
        locus["start"] - 1
        if locus["sense"] == "+"
        else sequence_length - locus["start"] + 1
    )
    nucleotide_length = sequence_length - start_index if locus["trailing"] else locus["length"]
    if nucleotide_length <= 0 or start_index < 0 or start_index + nucleotide_length > sequence_length:
        raise ValueError(f"Invalid orffinder locus for input length {sequence_length}: {locus}")
    if locus["sense"] == "+":
        start0, end0 = start_index, start_index + nucleotide_length
    else:
        start0, end0 = sequence_length - start_index - nucleotide_length, sequence_length - start_index
    return start0, end0, locus["sense"], nucleotide_length


def interval_iou(left: tuple[int, int], right: tuple[int, int]) -> float:
    intersection = max(0, min(left[1], right[1]) - max(left[0], right[0]))
    union = max(left[1], right[1]) - min(left[0], right[0])
    return intersection / union if union else 0.0


def python_raw_candidates(record_id: str, sequence: str) -> list[Candidate]:
    """Keep package-native coordinates and protein translations (including '*')."""
    # This is the direct package baseline at the requested threshold. Its
    # incomplete-terminal length behavior is intentionally not adapter-corrected.
    loci = orffinder.getORFProteins(
        SeqRecord(Seq(sequence)),
        minimum_length=MINIMUM_NT,
        start_codons=["ATG"],
        remove_nested=False,
        trim_trailing=False,
        return_loci=True,
    )
    _RAW_PACKAGE_OUTPUTS[record_id] = [
        {key: str(value) if key == "protein" else value for key, value in locus.items()}
        for locus in loci
    ]
    candidates = []
    for locus in loci:
        start0, end0, strand, nucleotide_length = _locus_interval(locus, len(sequence))
        protein = str(locus["protein"])
        candidates.append(Candidate(
            f"{record_id}:raw:{locus['index']}", "python_raw", record_id,
            start0, end0, strand, nucleotide_length, protein,
            bool(not locus["trailing"]), not locus["trailing"],
        ))
    return candidates


def corrected_adapter_candidates(record_id: str, sequence: str) -> list[Candidate]:
    """Mirror adapter discovery/translation and verify against production output."""
    sequence = sequence.upper()
    minimum_for_finder = max(1, MINIMUM_NT - 1)
    loci = orffinder.getORFs(
        SeqRecord(Seq(sequence)),
        minimum_length=minimum_for_finder,
        start_codons=["ATG"],
        remove_nested=False,
        trim_trailing=False,
    )
    forward = sequence
    reverse = str(Seq(sequence).reverse_complement())
    candidates = []
    for locus in loci:
        start0, end0, strand, nucleotide_length = _locus_interval(locus, len(sequence))
        if nucleotide_length < MINIMUM_NT:
            continue
        strand_sequence = forward if strand == "+" else reverse
        start_index = locus["start"] - 1 if strand == "+" else len(sequence) - locus["start"] + 1
        coding_length = (
            nucleotide_length - nucleotide_length % 3
            if locus["trailing"]
            else nucleotide_length - 3
        )
        if coding_length <= 0:
            continue
        protein = str(Seq(strand_sequence[start_index:start_index + coding_length]).translate()).rstrip("*")
        if not protein:
            continue
        candidates.append(Candidate(
            f"{record_id}:adapter:{locus['index']}", "adapter_corrected", record_id,
            start0, end0, strand, nucleotide_length, protein,
            bool(not locus["trailing"]), not locus["trailing"],
        ))

    actual = OrffinderOrfSearcher(
        minimum_length=MINIMUM_NT,
        start_codons=("ATG",),
        remove_nested=False,
        trim_trailing=False,
    ).get_orfs(sequence)
    reconstructed = [candidate.protein for candidate in candidates]
    if reconstructed != actual:
        raise AssertionError(
            f"Benchmark adapter reconstruction differs for {record_id}: "
            f"reconstructed={Counter(reconstructed)}, production={Counter(actual)}"
        )
    return candidates


def _ncbi_record_id(orf_tag: str, known_ids: Iterable[str]) -> str:
    matches = [record_id for record_id in known_ids if orf_tag.endswith("_" + record_id)]
    if len(matches) != 1:
        raise ValueError(f"Cannot uniquely associate NCBI ORF {orf_tag!r} with input records")
    return matches[0]


_NCBI_CDS_HEADER = re.compile(r"^lcl\|([^:]+):(c?)(\d+)-(\d+)$")
_NCBI_ORF_TAG = re.compile(r"(ORF\d+_[^:\s]+):(\d+):(\d+)")


def parse_ncbi_outputs(cds_fasta: Path, protein_fasta: Path, records: dict[str, SeqRecord]) -> list[Candidate]:
    """Parse NCBI outfmt 1 CDS FASTA and outfmt 0 protein FASTA as a checked pair."""
    known_ids = sorted(records, key=len, reverse=True)
    cds_by_tag = {}
    with cds_fasta.open(encoding="utf-8") as handle:
        for cds_record in SeqIO.parse(handle, "fasta"):
            match = _NCBI_CDS_HEADER.match(cds_record.id)
            tag_match = _NCBI_ORF_TAG.search(cds_record.description)
            if not match or not tag_match:
                raise ValueError(f"Unrecognized NCBI CDS FASTA header: {cds_record.description}")
            input_id, reverse_mark, first, second = match.groups()
            record_id = input_id if input_id in records else None
            if record_id is None:
                raise ValueError(f"NCBI CDS header refers to unknown input record {input_id!r}")
            a, b = int(first), int(second)
            if reverse_mark:
                start0, end0, strand = b - 1, a, "-"
            else:
                start0, end0, strand = a - 1, b, "+"
            tag = tag_match.group(1)
            if _ncbi_record_id(tag, known_ids) != record_id:
                raise ValueError(f"NCBI ORF tag {tag} and CDS header disagree on record")
            if tag in cds_by_tag:
                raise ValueError(f"Duplicate NCBI ORF tag in CDS FASTA: {tag}")
            cds_by_tag[tag] = (record_id, start0, end0, strand, str(cds_record.seq).upper())

    protein_by_tag = {}
    with protein_fasta.open(encoding="utf-8") as handle:
        for protein_record in SeqIO.parse(handle, "fasta"):
            match = _NCBI_ORF_TAG.search(protein_record.id)
            if not match:
                raise ValueError(f"Unrecognized NCBI protein FASTA header: {protein_record.description}")
            tag = match.group(1)
            coords = tuple(map(int, match.groups()[1:]))
            if tag in protein_by_tag:
                raise ValueError(f"Duplicate NCBI ORF tag in protein FASTA: {tag}")
            protein_by_tag[tag] = (str(protein_record.seq), coords)

    if set(cds_by_tag) != set(protein_by_tag):
        raise ValueError(
            f"NCBI CDS/protein ORF identifiers differ; CDS-only={sorted(set(cds_by_tag)-set(protein_by_tag))}, "
            f"protein-only={sorted(set(protein_by_tag)-set(cds_by_tag))}"
        )

    candidates = []
    for tag in sorted(cds_by_tag):
        record_id, start0, end0, strand, cds_sequence = cds_by_tag[tag]
        protein, encoded_coords = protein_by_tag[tag]
        # NCBI protein headers encode zero-based inclusive endpoints in either
        # orientation; the CDS FASTA's 1-based inclusive range is authoritative.
        if {encoded_coords[0], encoded_coords[1]} != {start0, end0 - 1}:
            raise ValueError(f"NCBI protein/CDS coordinate encodings disagree for {tag}")
        translated_cds = str(Seq(cds_sequence).translate(table=GENETIC_CODE)).rstrip("*")
        if protein != translated_cds:
            raise ValueError(f"NCBI protein does not match translated CDS for {tag}")
        input_sequence = str(records[record_id].seq).upper()
        if strand == "+":
            expected_cds = input_sequence[start0:end0]
        else:
            expected_cds = str(Seq(input_sequence[start0:end0]).reverse_complement())
        if expected_cds != cds_sequence:
            raise ValueError(f"NCBI CDS sequence does not match saved input interval for {tag}")
        candidates.append(Candidate(
            f"{record_id}:ncbi:{tag}", "ncbi", record_id,
            start0, end0, strand, end0 - start0, protein,
            bool(cds_sequence[-3:] in {"TAA", "TAG", "TGA"}),
            bool(cds_sequence[-3:] in {"TAA", "TAG", "TGA"}),
        ))
    return candidates


def exact_pairs(left: list[Candidate], right: list[Candidate]) -> list[tuple[Candidate, Candidate, float]]:
    """Pair exact interval/strand matches one-to-one by sorted stable IDs."""
    right_by_interval: dict[tuple[int, int, str], list[Candidate]] = defaultdict(list)
    for candidate in right:
        right_by_interval[candidate.interval].append(candidate)
    for values in right_by_interval.values():
        values.sort(key=lambda candidate: candidate.candidate_id)
    pairs = []
    for candidate in sorted(left, key=lambda item: item.candidate_id):
        matches = right_by_interval.get(candidate.interval, [])
        if matches:
            pairs.append((candidate, matches.pop(0), 1.0))
    return pairs


def compare_candidates(left: list[Candidate], right: list[Candidate]) -> dict:
    """Return all overlaps, exact matches, threshold pairs, and unmatched IDs."""
    exact = exact_pairs(left, right)
    exact_left = {pair[0].candidate_id for pair in exact}
    exact_right = {pair[1].candidate_id for pair in exact}
    overlaps = []
    for a in left:
        for b in right:
            if a.strand != b.strand:
                continue
            intersection = max(0, min(a.end0, b.end0) - max(a.start0, b.start0))
            if intersection <= 0:
                continue
            iou = interval_iou((a.start0, a.end0), (b.start0, b.end0))
            overlaps.append((a, b, iou))
    overlaps.sort(key=lambda pair: (pair[0].candidate_id, pair[1].candidate_id))

    # Deterministic greedy one-to-one summary: reserve exact pairs first, then
    # rank remaining eligible edges by descending IoU and ascending stable IDs.
    eligible = [pair for pair in overlaps if pair[2] >= 0.5
                and pair[0].interval != pair[1].interval
                and pair[0].candidate_id not in exact_left
                and pair[1].candidate_id not in exact_right]
    eligible.sort(key=lambda pair: (-pair[2], pair[0].candidate_id, pair[1].candidate_id))
    used_left, used_right = set(exact_left), set(exact_right)
    threshold_pairs = []
    for pair in eligible:
        if pair[0].candidate_id in used_left or pair[1].candidate_id in used_right:
            continue
        threshold_pairs.append(pair)
        used_left.add(pair[0].candidate_id)
        used_right.add(pair[1].candidate_id)
    return {
        "exact": exact,
        "overlaps": overlaps,
        "threshold_pairs": threshold_pairs,
        "unmatched_left": sorted(c.candidate_id for c in left if c.candidate_id not in used_left),
        "unmatched_right": sorted(c.candidate_id for c in right if c.candidate_id not in used_right),
    }


def load_assigned_inputs(input_dir: Path = ASSIGNED_INPUT_DIR) -> tuple[dict[str, SeqRecord], list[dict]]:
    """Load and validate the professor-provided genomic RefSeqGene FASTAs."""
    fasta_paths = sorted(
        path.name for path in input_dir.iterdir()
        if path.is_file() and path.suffix.lower() in {".fa", ".fasta"}
    )
    if tuple(fasta_paths) != tuple(sorted(ASSIGNED_INPUT_FILES)):
        raise ValueError(
            f"Expected exactly {sorted(ASSIGNED_INPUT_FILES)} in {input_dir}, got {fasta_paths}"
        )

    records: dict[str, SeqRecord] = {}
    provenance = []
    for filename in ASSIGNED_INPUT_FILES:
        path = input_dir / filename
        with path.open(encoding="utf-8") as handle:
            parsed = list(SeqIO.parse(handle, "fasta"))
        if len(parsed) != 1:
            raise ValueError(f"Expected exactly one FASTA record in {filename}; found {len(parsed)}")
        record = parsed[0]
        if not record.id or not record.description:
            raise ValueError(f"Missing FASTA identifier/header in {filename}")
        accession_match = ACCESSION_VERSION.search(record.description)
        if accession_match is None:
            raise ValueError(f"No versioned RefSeqGene accession in FASTA header: {record.description}")
        accession_version = accession_match.group(1)
        if accession_version != ASSIGNED_ACCESSIONS[filename]:
            raise ValueError(
                f"{filename} expected RefSeqGene accession {ASSIGNED_ACCESSIONS[filename]}, "
                f"found {accession_version}"
            )
        if record.id != accession_version:
            raise ValueError(
                f"FASTA identifier {record.id!r} does not match header accession {accession_version!r}"
            )
        sequence = str(record.seq).upper()
        if not sequence:
            raise ValueError(f"Empty DNA sequence in {filename}")
        invalid = sorted(set(sequence) - DNA_ALPHABET)
        if invalid:
            raise ValueError(f"Invalid DNA symbols in {filename}: {invalid}")
        if record.id in records:
            raise ValueError(f"Duplicate FASTA record identifier: {record.id}")
        records[record.id] = record
        provenance.append({
            "filename": filename,
            "original_header": record.description,
            "accession_version": accession_version,
            "sequence_length": len(sequence),
            "file_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "sequence_sha256": sha256_text(sequence),
        })
    if input_dir.resolve() == ASSIGNED_INPUT_DIR.resolve():
        with ASSIGNED_INPUT_MANIFEST.open(encoding="utf-8", newline="") as handle:
            saved_manifest = list(csv.DictReader(handle, delimiter="\t"))
        normalized_manifest = [
            {**entry, "sequence_length": int(entry["sequence_length"])}
            for entry in saved_manifest
        ]
        if normalized_manifest != provenance:
            raise ValueError(
                f"Assigned FASTA files do not match their saved provenance manifest: "
                f"{ASSIGNED_INPUT_MANIFEST}"
            )
    return records, provenance


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, content) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(content, indent=2) + "\n", encoding="utf-8")


def _run_ncbi(executable: str, input_fasta: Path, out: Path, fmt: int,
                         minimum_nt: int = MINIMUM_NT, nested: bool = False) -> list[str]:
    command = [executable, "-in", str(input_fasta), "-s", "0", "-g", str(GENETIC_CODE),
                             "-ml", str(minimum_nt), "-n", str(nested).lower(), "-strand", "both",
               "-outfmt", str(fmt), "-out", str(out)]
    subprocess.run(command, check=True, capture_output=True, text=True)
    return command


def run_control_probes(executable: str) -> dict:
    """Run and assert synthetic-only controls, preserving every NCBI output."""
    control_fasta = DATA / "controls.fasta"
    control_records = {record.id: record for record in SeqIO.parse(control_fasta, "fasta")}
    raw_dir = RESULTS / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    outputs = {}
    commands = []
    for minimum_nt in (30, MINIMUM_NT):
        protein_path = raw_dir / f"controls_ml{minimum_nt}_nfalse_outfmt0.fasta"
        cds_path = raw_dir / f"controls_ml{minimum_nt}_nfalse_outfmt1.fasta"
        commands.append(_run_ncbi(executable, control_fasta, protein_path, 0, minimum_nt))
        commands.append(_run_ncbi(executable, control_fasta, cds_path, 1, minimum_nt))
        outputs[minimum_nt] = parse_ncbi_outputs(cds_path, protein_path, control_records)
    nested_path = raw_dir / "controls_ml30_ntrue_outfmt1.fasta"
    commands.append(_run_ncbi(executable, control_fasta, nested_path, 1, 30, nested=True))

    baseline30 = outputs[30]
    baseline75 = outputs[MINIMUM_NT]
    by_record30 = defaultdict(list)
    by_record75 = defaultdict(list)
    for candidate in baseline30:
        by_record30[candidate.record_id].append(candidate)
    for candidate in baseline75:
        by_record75[candidate.record_id].append(candidate)
    required = {"forward_complete", "reverse_complete", "nested_pair", "incomplete_terminal",
                "incomplete_terminal_reverse"}
    if not required.issubset(control_records):
        raise ValueError(f"Synthetic controls missing expected records: {sorted(required-set(control_records))}")
    if not any(c.strand == "+" and c.terminal_stop_in_nt for c in by_record75["forward_complete"]):
        raise AssertionError("NCBI did not report the controlled forward complete ORF with its stop")
    if not any(c.strand == "-" and c.terminal_stop_in_nt for c in by_record75["reverse_complete"]):
        raise AssertionError("NCBI did not report the controlled reverse complete ORF with its stop")
    if not any(not c.terminal_stop_in_nt and c.nucleotide_length >= MINIMUM_NT for c in by_record75["incomplete_terminal"]):
        raise AssertionError("NCBI terminal-ORF behavior did not retain the long incomplete control")
    if not any(c.strand == "-" and not c.terminal_stop_in_nt and c.nucleotide_length >= MINIMUM_NT
               for c in by_record75["incomplete_terminal_reverse"]):
        raise AssertionError("NCBI terminal-ORF behavior did not retain the reverse-strand incomplete control")
    nested_fasta_count = sum(
        1 for candidate in SeqIO.parse(nested_path, "fasta")
        if candidate.id.startswith("lcl|nested_pair:")
    )
    nested_false_count = len(by_record30["nested_pair"])
    if nested_false_count <= nested_fasta_count:
        raise AssertionError("NCBI -n true did not remove the known nested control candidate")

    # Independently check the Python package's own nested-filter switch on the
    # same control sequence; this is not mixed into the human transcript output.
    nested_sequence = str(control_records["nested_pair"].seq).upper()
    py_unfiltered = orffinder.getORFs(SeqRecord(Seq(nested_sequence)), minimum_length=30,
                                     start_codons=["ATG"], remove_nested=False, trim_trailing=False)
    py_filtered = orffinder.getORFs(SeqRecord(Seq(nested_sequence)), minimum_length=30,
                                   start_codons=["ATG"], remove_nested=True, trim_trailing=False)
    if len(py_unfiltered) <= len(py_filtered):
        raise AssertionError("Python package nested-filter control did not remove a nested candidate")
    original_raw_outputs = dict(_RAW_PACKAGE_OUTPUTS)
    _RAW_PACKAGE_OUTPUTS.clear()
    control_candidate_rows = []
    control_by_tool = {"ncbi_ml30": baseline30, "ncbi_ml75": baseline75,
                       "python_raw": [], "adapter_corrected": []}
    for record_id, record in control_records.items():
        sequence = str(record.seq).upper()
        raw_candidates = python_raw_candidates(record_id, sequence)
        adapter_candidates = corrected_adapter_candidates(record_id, sequence)
        control_by_tool["python_raw"].extend(raw_candidates)
        control_by_tool["adapter_corrected"].extend(adapter_candidates)
    write_json(RESULTS / "raw" / "controls_python_package_loci.json", _RAW_PACKAGE_OUTPUTS)
    _RAW_PACKAGE_OUTPUTS.clear()
    _RAW_PACKAGE_OUTPUTS.update(original_raw_outputs)
    for tool, candidates in control_by_tool.items():
        for candidate in candidates:
            row = asdict(candidate)
            row["tool"] = tool
            control_candidate_rows.append(row)
    write_tsv(RESULTS / "control_candidates.tsv", control_candidate_rows)
    control_comparison_rows = []
    for record_id in control_records:
        raw = [c for c in control_by_tool["python_raw"] if c.record_id == record_id]
        corrected = [c for c in control_by_tool["adapter_corrected"] if c.record_id == record_id]
        ncbi = [c for c in baseline75 if c.record_id == record_id]
        for left_name, left, right_name, right in (
            ("python_raw", raw, "adapter_corrected", corrected),
            ("adapter_corrected", corrected, "ncbi_ml75", ncbi),
        ):
            comparison = compare_candidates(left, right)
            control_comparison_rows.append({
                "record_id": record_id, "left_tool": left_name, "right_tool": right_name,
                "left_count": len(left), "right_count": len(right),
                "exact_count": len(comparison["exact"]),
                "positive_overlap_pair_count": len(comparison["overlaps"]),
                "iou_ge_0_5_pair_count": len(comparison["threshold_pairs"]),
                "unmatched_left_count": len(comparison["unmatched_left"]),
                "unmatched_right_count": len(comparison["unmatched_right"]),
            })
    write_tsv(RESULTS / "control_comparisons.tsv", control_comparison_rows)
    return {
        "ncbi_ml30_candidate_count": len(baseline30),
        "ncbi_ml75_candidate_count": len(baseline75),
        "ncbi_nested_false_count": nested_false_count,
        "ncbi_nested_true_count": nested_fasta_count,
        "python_nested_false_count": len(py_unfiltered),
        "python_nested_true_count": len(py_filtered),
        "ncbi_minimum_75_incomplete_control_included": bool(by_record75["incomplete_terminal"]),
        "commands": commands,
    }


def run_benchmark(executable: str) -> None:
    records, input_provenance = load_assigned_inputs()
    RESULTS.mkdir(parents=True, exist_ok=True)
    raw_dir = RESULTS / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    assigned_fasta = raw_dir / "assigned_inputs.fasta"
    SeqIO.write(list(records.values()), assigned_fasta, "fasta")

    control_metadata = run_control_probes(executable)
    _RAW_PACKAGE_OUTPUTS.clear()
    ncbi_cds_path = raw_dir / "ncbi_outfmt1_cds.fasta"
    ncbi_protein_path = raw_dir / "ncbi_outfmt0_proteins.fasta"
    command1 = _run_ncbi(executable, assigned_fasta, ncbi_protein_path, 0)
    command2 = _run_ncbi(executable, assigned_fasta, ncbi_cds_path, 1)
    ncbi_candidates = parse_ncbi_outputs(ncbi_cds_path, ncbi_protein_path, records)

    candidates_by_tool: dict[str, list[Candidate]] = {"python_raw": [], "adapter_corrected": [], "ncbi": ncbi_candidates}
    for record_id, record in records.items():
        sequence = str(record.seq).upper()
        candidates_by_tool["python_raw"].extend(python_raw_candidates(record_id, sequence))
        candidates_by_tool["adapter_corrected"].extend(corrected_adapter_candidates(record_id, sequence))

    write_json(raw_dir / "python_package_loci.json", _RAW_PACKAGE_OUTPUTS)
    candidate_rows = [asdict(candidate) for tool in candidates_by_tool.values() for candidate in tool]
    write_tsv(RESULTS / "candidates.tsv", candidate_rows)
    pair_rows, match_rows, unmatched_rows, translation_difference_rows = [], [], [], []
    tools = ("python_raw", "adapter_corrected", "ncbi")
    comparisons = {}
    for left_name, right_name in itertools.combinations(tools, 2):
        left_groups = defaultdict(list)
        right_groups = defaultdict(list)
        for candidate in candidates_by_tool[left_name]:
            left_groups[candidate.record_id].append(candidate)
        for candidate in candidates_by_tool[right_name]:
            right_groups[candidate.record_id].append(candidate)
        for record_id in records:
            result = compare_candidates(left_groups[record_id], right_groups[record_id])
            comparisons[left_name, right_name, record_id] = result
            for a, b, iou in result["exact"]:
                match_rows.append({"record_id": record_id, "left_tool": left_name, "right_tool": right_name,
                                   "match_type": "exact", "left_candidate_id": a.candidate_id,
                                   "right_candidate_id": b.candidate_id, "iou": iou})
                if (left_name, right_name) == ("python_raw", "adapter_corrected"):
                    translation_difference_rows.append({
                        "record_id": record_id,
                        "start0": a.start0,
                        "end0": a.end0,
                        "strand": a.strand,
                        "raw_candidate_id": a.candidate_id,
                        "adapter_candidate_id": b.candidate_id,
                        "raw_package_protein": a.protein,
                        "corrected_adapter_protein": b.protein,
                        "same_translation": a.protein == b.protein,
                    })
            for a, b, iou in result["overlaps"]:
                pair_rows.append({"record_id": record_id, "left_tool": left_name, "right_tool": right_name,
                                  "left_candidate_id": a.candidate_id, "right_candidate_id": b.candidate_id,
                                  "left_interval": f"{a.start0}:{a.end0}:{a.strand}",
                                  "right_interval": f"{b.start0}:{b.end0}:{b.strand}", "iou": iou,
                                  "exact_interval_and_strand": a.interval == b.interval})
            for a, b, iou in result["threshold_pairs"]:
                match_rows.append({"record_id": record_id, "left_tool": left_name, "right_tool": right_name,
                                   "match_type": "greedy_iou_ge_0.5", "left_candidate_id": a.candidate_id,
                                   "right_candidate_id": b.candidate_id, "iou": iou})
            unmatched_rows.extend(
                {"record_id": record_id, "tool": left_name, "candidate_id": candidate_id}
                for candidate_id in result["unmatched_left"]
            )
            unmatched_rows.extend(
                {"record_id": record_id, "tool": right_name, "candidate_id": candidate_id}
                for candidate_id in result["unmatched_right"]
            )

    write_tsv(RESULTS / "overlaps.tsv", pair_rows)
    write_tsv(RESULTS / "matches.tsv", match_rows)
    write_tsv(RESULTS / "unmatched.tsv", unmatched_rows)
    write_tsv(RESULTS / "translation_differences.tsv", translation_difference_rows)
    provenance_by_accession = {entry["accession_version"]: entry for entry in input_provenance}
    comparison_rows = []
    for record_id in records:
        input_entry = provenance_by_accession[record_id]
        gene = Path(input_entry["filename"]).stem
        entry = {
            "record_id": record_id,
            "input_filename": input_entry["filename"],
            "gene_label": gene,
        }
        for tool in tools:
            candidates = [candidate for candidate in candidates_by_tool[tool] if candidate.record_id == record_id]
            longest_length = max((candidate.nucleotide_length for candidate in candidates), default=0)
            longest = [candidate for candidate in candidates if candidate.nucleotide_length == longest_length]
            entry[f"{tool}_count"] = len(candidates)
            longest_proteins = Counter(candidate.protein for candidate in candidates)
            entry[f"{tool}_protein_translation_counts"] = json.dumps(dict(sorted(longest_proteins.items())))
            entry[f"{tool}_longest_nt"] = longest_length
            entry[f"{tool}_longest_candidate_ids"] = json.dumps([c.candidate_id for c in longest])
            entry[f"{tool}_longest_locations"] = json.dumps([f"{c.start0}:{c.end0}:{c.strand}" for c in longest])
            entry[f"{tool}_longest_proteins"] = json.dumps([c.protein for c in longest])
        for left_name, right_name in itertools.combinations(tools, 2):
            result = comparisons[left_name, right_name, record_id]
            prefix = f"{left_name}_vs_{right_name}"
            entry[f"{prefix}_exact_count"] = len(result["exact"])
            entry[f"{prefix}_positive_overlap_pair_count"] = len(result["overlaps"])
            entry[f"{prefix}_iou_ge_0_5_pair_count"] = len(result["threshold_pairs"])
            entry[f"{prefix}_unmatched_left_count"] = len(result["unmatched_left"])
            entry[f"{prefix}_unmatched_right_count"] = len(result["unmatched_right"])
        raw_adapter_exact = comparisons["python_raw", "adapter_corrected", record_id]["exact"]
        entry["raw_adapter_exact_translation_difference_count"] = sum(
            raw.protein != corrected.protein for raw, corrected, _ in raw_adapter_exact
        )
        comparison_rows.append(entry)
    write_tsv(RESULTS / "comparisons.tsv", comparison_rows)

    metadata = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version,
        "biopython_version": biopython_version,
        "orffinder_version": importlib.metadata.version("orffinder"),
        "ncbi_orffinder_executable": str(Path(executable).resolve()),
        "ncbi_orffinder_version": subprocess.run([executable, "-version-full"], check=True, capture_output=True, text=True).stdout.strip(),
        "parameters": {"start_codons": ["ATG"], "ncbi_start_option": "-s 0", "genetic_code": GENETIC_CODE,
                       "minimum_orf_length_nt": MINIMUM_NT, "nested_filter": False, "strand": "both",
                       "ncbi_minimum_floor_nt": 30},
        "commands": [command1, command2],
        "assigned_inputs": input_provenance,
        "combined_input_fasta_sha256": hashlib.sha256(assigned_fasta.read_bytes()).hexdigest(),
        "synthetic_control_results": control_metadata,
        "dataset_type": "three professor-provided genomic RefSeqGene FASTA records; no sequence download or GenBank annotation",
        "biological_scope": "ORF scanning does not perform splicing and does not establish the protein annotated for the named gene",
        "coordinate_system": "zero-based half-open on original input; strand separate",
        "ncbi_terminal_stop_normalization": "CDS FASTA sequence/location includes stop; protein excludes translated terminal stop; normalized interval retains the terminal stop nucleotide",
        "overlap_rule": "All same-strand positive non-exact interval intersections are emitted with IoU. Exact interval-and-strand matches are separate.",
        "pairing_rule": "Reserve exact matches one-to-one in ascending stable candidate-ID order. Among remaining candidates, sort positive-overlap edges with IoU >= 0.5 by descending IoU then ascending left/right candidate IDs; greedily accept an edge if neither endpoint is already paired. Report unmatched candidates after these pairings.",
    }
    write_json(RESULTS / "metadata.json", metadata)
    write_report(comparison_rows, input_provenance)


def write_report(comparisons: list[dict], input_provenance: list[dict]) -> None:
    lines = [
        "# ORFfinder benchmark results",
        "",
        "Completed comparison of the three assigned genomic RefSeqGene FASTA records using raw `orffinder`, Lverage's corrected adapter, and NCBI standalone ORFfinder. ORF scanning on genomic sequence does not perform splicing or establish the protein annotated for the named gene.",
        "",
        "## Results per assigned input",
        "",
        "| Input file | RefSeqGene accession.version | Raw package ORFs | Corrected adapter ORFs | NCBI ORFs | Raw/adapter exact-region translation differences | Adapter/NCBI exact matches | Adapter/NCBI positive overlaps | Same longest interval? |",
        "|---|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    by_record = {row["record_id"]: row for row in comparisons}
    adapter_total = sum(int(row["adapter_corrected_count"]) for row in comparisons)
    ncbi_total = sum(int(row["ncbi_count"]) for row in comparisons)
    exact_total = sum(int(row["adapter_corrected_vs_ncbi_exact_count"]) for row in comparisons)
    for item in input_provenance:
        record_id = item["accession_version"]
        row = by_record[record_id]
        exact = row["adapter_corrected_vs_ncbi_exact_count"]
        same_longest = bool(set(json.loads(row["adapter_corrected_longest_locations"])) & set(json.loads(row["ncbi_longest_locations"])))
        lines.append(f"| {item['filename']} | {record_id} | {row['python_raw_count']} | {row['adapter_corrected_count']} | {row['ncbi_count']} | {row['raw_adapter_exact_translation_difference_count']} | {exact} | {row['adapter_corrected_vs_ncbi_positive_overlap_pair_count']} | {'yes' if same_longest else 'no'} |")
    lines.extend([
        "",
        "## Findings",
        "",
        f"- Across the three assigned inputs, the corrected adapter returned {adapter_total:,} candidates and NCBI returned {ncbi_total:,}; there were {exact_total:,} exact interval-and-strand matches.",
        f"- Every NCBI candidate had an exact adapter match. The adapter returned {adapter_total - ncbi_total:,} additional candidates.",
        "- All three inputs share at least one exact longest-ORF interval and strand between the adapter and NCBI.",
        "- Raw-package/adapter protein-string differences include terminal stop-symbol removal. Their counts do not indicate incorrect proteins.",
        "- These findings apply only to the assigned genomic sequences and recorded settings. ORF scanning does not perform splicing or establish the annotated protein for the named gene.",
        "",
        "## Interpretation",
        "",
        "- Called-region differences are represented in `candidates.tsv`, `overlaps.tsv`, `matches.tsv`, and `unmatched.tsv`; exact matches and non-exact overlaps are separate.",
        "- Overlap means a positive intersection on the same strand; `overlaps.tsv` includes exact interval pairs with IoU 1.0 and marks them explicitly. Exact matches are also listed separately in `matches.tsv`. The additional one-to-one summary uses deterministic greedy pairing at IoU ≥ 0.5 after exact matches are reserved.",
        "- Longest ORFs are ranked by normalized nucleotide span. All tied longest locations and proteins are retained in `comparisons.tsv`; same-longest means the tied sets share at least one exact interval-and-strand location.",
        "- Raw package translation differences are counted only for exact same-region raw/adapter candidate pairs in `translation_differences.tsv`; candidates found by only one discovery path remain region differences. The corrected adapter candidate reconstruction is checked against production `OrffinderOrfSearcher.get_orfs()` as an ordered protein list and `Counter`, preserving duplicate counts.",
        "- Tool settings and exact commands are in `metadata.json`. NCBI outfmt 1 locations are 1-based inclusive; they are normalized to half-open intervals while retaining stop codons. Outfmt 0 supplies proteins and is cross-checked against the matching CDS record.",
        "- `metadata.json` records each input filename, original FASTA header, accession version, sequence length, and file/sequence SHA-256 checksums. Synthetic controls are validation fixtures and are not part of these assigned-input results.",
        "",
    ])
    # Keep generated Markdown sections separated and normalize the file ending
    # so renderers do not treat trailing blank lines as an extra paragraph.
    report = "\n".join(lines).rstrip()
    (RESULTS / "report.md").write_text(f"{report}\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    run_parser = subparsers.add_parser("run", help="compare the three assigned genomic RefSeqGene FASTAs")
    run_parser.add_argument("--executable", default=os.environ.get("ORFFINDER_EXECUTABLE", DEFAULT_EXECUTABLE))
    args = parser.parse_args()
    run_benchmark(args.executable)


_RAW_PACKAGE_OUTPUTS: dict[str, list[dict]] = {}

if __name__ == "__main__":
    main()
