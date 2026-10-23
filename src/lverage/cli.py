"""Command-line interface for the V2 evidence-preserving pipeline."""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys


DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN")
FASTA_SUFFIXES = frozenset({".fa", ".fasta", ".fna"})
RESULT_HEADERS = [
    "Row Type", "Input Gene", "Query ORF ID", "Query ORF Selected",
    "Stage", "Status", "Error", "Outcome Detail", "Attempted Hit Identifier",
    "Ortholog Hit ID", "Ortholog Accession",
    "Ortholog Description", "Ortholog Species Name", "Ortholog Species Tax ID",
    "BLAST E-value", "Ortholog Identity Ratio", "Ortholog Identity Percent",
    "Query Coverage Ratio", "Query Coverage Percent",
    "Query Domain Name", "Query Domain ID", "Query Domain Start 0-based",
    "Query Domain End 0-based-exclusive", "Ortholog Domain Name", "Ortholog Domain ID",
    "Ortholog Domain Start 0-based", "Ortholog Domain End 0-based-exclusive",
    "Domain Pair Identity Ratio", "Domain Pair Identity Percent", "Threshold Ratio",
    "Threshold Passed", "Motif Database", "Matrix ID", "Motif Name", "Motif Class",
    "Motif Inference E-value", "Motif URL", "PFM JSON",
]


class InputValidationError(Exception):
    """A user-provided input or setting failed preflight validation."""


def _positive_float(value):
    try:
        parsed = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be a number") from error
    if not math.isfinite(parsed) or parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than zero")
    return parsed


def _identity_ratio(value):
    try:
        parsed = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be a number") from error
    if not math.isfinite(parsed) or not 0 <= parsed <= 1:
        raise argparse.ArgumentTypeError("must be between zero and one")
    return parsed


def _positive_int(value):
    try:
        parsed = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be an integer") from error
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than zero")
    return parsed


def build_parser():
    """Build argument parser without loading any runtime adapters."""
    parser = argparse.ArgumentParser(
        prog="lverage",
        description=(
            "Run V2 on one nucleotide FASTA file per gene. --output-dir creates "
            "a results.tsv, orfs.fasta, and metadata.json evidence bundle."
        ),
        epilog=(
            "Rows are explicitly typed: input, ORF, query_domain, ortholog, "
            "ortholog_retrieval, ortholog_domain, domain_comparison, motif, stage_status, or gene_error. "
            "Earlier-stage rows are not ortholog/motif hits. --stop-after alignment "
            "is an intentional pre-JASPAR stop. Exit 0 means completed or intentionally "
            "stopped (including biological no-hits), 2 is preflight failure, and 1 "
            "means one or more execution failures with partial evidence preserved."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input-dir", required=True, metavar="DIR", help="directory of per-gene nucleotide FASTA files")
    parser.add_argument("--blast-db", required=True, metavar="PREFIX", help="local protein BLAST database prefix")
    parser.add_argument("--species-map", required=True, metavar="JSON", help="JSON mapping species names to NCBI tax IDs")
    parser.add_argument("--pfam-db", required=True, metavar="DIR", help="local PfamScan database directory")
    parser.add_argument("--output-dir", required=True, metavar="DIR", help="new directory for results.tsv, orfs.fasta, and metadata.json")
    parser.add_argument("--blastp", default="blastp", metavar="EXE", help="blastp executable name or path")
    parser.add_argument("--blastdbcmd", default="blastdbcmd", metavar="EXE", help="blastdbcmd executable name or path")
    parser.add_argument("--pfamscan", default="pfam_scan.pl", metavar="EXE", help="pfam_scan.pl executable name or path")
    parser.add_argument("--identity-threshold", type=_identity_ratio, default=0.7, metavar="RATIO", help="minimum domain identity ratio")
    parser.add_argument("--valid-pfam", action="append", default=None, metavar="ACCESSION", help="accepted Pfam accession; repeatable")
    parser.add_argument("--jaspar-hits", type=_positive_int, default=10, metavar="N", help="maximum JASPAR motif hits")
    parser.add_argument("--jaspar-evalue", type=_positive_float, default=1e-6, metavar="VALUE", help="maximum JASPAR inference E-value")
    parser.add_argument("--stop-after", choices=("alignment",), help="intentionally stop after domain alignment, without contacting JASPAR")
    return parser


def _read_species_map(path):
    def reject_duplicate_keys(pairs):
        mapping = {}
        for key, value in pairs:
            if key in mapping:
                raise ValueError(f"duplicate JSON key: {key}")
            mapping[key] = value
        return mapping

    try:
        with path.open(encoding="utf-8") as handle:
            mapping = json.load(handle, object_pairs_hook=reject_duplicate_keys)
    except (OSError, UnicodeError, ValueError) as error:
        raise InputValidationError(f"cannot read species mapping {path}: {error}") from error
    if not isinstance(mapping, dict) or not mapping:
        raise InputValidationError("species mapping must be a nonempty JSON object")
    normalized = set()
    for name, tax_id in mapping.items():
        if not isinstance(name, str) or not name.strip() or name != name.strip():
            raise InputValidationError("species mapping keys must be nonempty, trimmed names")
        if name.casefold() in normalized:
            raise InputValidationError(f"duplicate species name ignoring case: {name}")
        normalized.add(name.casefold())
        if isinstance(tax_id, bool) or not isinstance(tax_id, int) or tax_id <= 0:
            raise InputValidationError(f"species mapping tax ID for {name!r} must be a positive integer")
    return mapping


def _read_fasta(path):
    from Bio import SeqIO
    try:
        with path.open(encoding="utf-8") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
    except Exception as error:
        raise InputValidationError(f"cannot parse FASTA {path}: {error}") from error
    if not records:
        raise InputValidationError(f"FASTA file contains no records: {path}")
    fragments = []
    for index, record in enumerate(records, 1):
        sequence = str(record.seq)
        if not sequence:
            raise InputValidationError(f"empty FASTA record {index} in {path}")
        invalid = sorted(set(sequence.upper()) - DNA_ALPHABET)
        if invalid:
            raise InputValidationError(f"invalid DNA symbols in FASTA record {index} of {path}: {''.join(invalid)}")
        fragments.append(sequence)
    return fragments


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _validate_inputs(args):
    input_dir = Path(args.input_dir)
    if not input_dir.is_dir():
        raise InputValidationError(f"input directory does not exist: {input_dir}")
    fasta_files = sorted(
        (item for item in input_dir.iterdir() if item.is_file() and item.suffix.casefold() in FASTA_SUFFIXES),
        key=lambda item: (item.name.casefold(), item.name),
    )
    if not fasta_files:
        raise InputValidationError(f"no .fa, .fasta, or .fna files found in {input_dir}")
    identifiers = set()
    inputs = []
    for fasta_path in fasta_files:
        gene = fasta_path.stem
        if gene.casefold() in identifiers:
            raise InputValidationError(f"duplicate gene identifier from filename stem: {gene}")
        identifiers.add(gene.casefold())
        inputs.append({"gene": gene, "path": str(fasta_path), "fragments": _read_fasta(fasta_path), "sha256": _sha256(fasta_path)})
    species_path = Path(args.species_map)
    species_map = _read_species_map(species_path)
    pfam_dir = Path(args.pfam_db)
    if not pfam_dir.is_dir():
        raise InputValidationError(f"PfamScan database directory does not exist: {pfam_dir}")
    if not args.blast_db.strip():
        raise InputValidationError("--blast-db must be a nonempty database prefix")
    if not args.blastp.strip() or not args.blastdbcmd.strip() or not args.pfamscan.strip():
        raise InputValidationError("executable paths must not be empty")
    if args.valid_pfam is not None and any(not value.strip() for value in args.valid_pfam):
        raise InputValidationError("--valid-pfam values must not be empty")
    output_dir = Path(args.output_dir)
    if os.path.lexists(output_dir):
        raise InputValidationError(f"output directory already exists: {output_dir}")
    if not output_dir.parent.is_dir():
        raise InputValidationError(f"output parent directory does not exist: {output_dir.parent}")
    return inputs, species_map, species_path


def _build_pipeline(args, species_map):
    from .blast import LocalBlastSearcher
    from .jaspar import Jaspar2024MotifDB
    from .orf_searcher import OrffinderOrfSearcher
    from .pfam import LocalPfamScanner
    from .pipeline import Lverage
    ortholog_searcher = LocalBlastSearcher(
        database_path=args.blast_db, species_map=species_map,
        blastp_path=args.blastp, blastdbcmd_path=args.blastdbcmd,
    )
    return Lverage(
        motif_database_list=[Jaspar2024MotifDB(n_hits=args.jaspar_hits, escore_threshold=args.jaspar_evalue)],
        orf_searcher=OrffinderOrfSearcher(),
        domain_scanner=LocalPfamScanner(database_path=args.pfam_db, pfamscan_path=args.pfamscan),
        ortholog_searcher=ortholog_searcher,
        valid_pfam_list=args.valid_pfam or [], dbd_identity_thresh=args.identity_threshold,
    )


def _timestamp():
    return datetime.now(timezone.utc).isoformat()


def _redact_error(error):
    text = re.sub(r"https?://[^\s]+", "<url-redacted>", str(error))
    text = re.sub(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}", "<email-redacted>", text)
    text = re.sub(r"(?i)(token|password|secret|api[_-]?key)=([^\s&]+)", r"\1=<redacted>", text)
    text = re.sub(r"(?i)bearer\s+[A-Za-z0-9._~+/-]+=*", "Bearer <redacted>", text)
    return text


def _orf_id(sequence):
    return "orf_" + hashlib.sha256(sequence.encode("ascii")).hexdigest()


def _ratio_and_percent(value):
    return value, value * 100 if value is not None else None


def _event_row(event, gene, sequence_ids, selection_status, selection_details):
    row = {header: "" for header in RESULT_HEADERS}
    row.update({"Input Gene": gene, "Stage": event.get("stage", ""), "Status": event.get("status", ""), "Error": event.get("error", ""), "Outcome Detail": event.get("detail", "")})
    kind = event["event"]
    if kind == "gene_start":
        row["Row Type"] = "input"
    elif kind == "orf_candidate":
        sequence = event["orf_sequence"]
        row["Row Type"] = "orf"
        row["Query ORF ID"] = sequence_ids[sequence]
        row["Query ORF Selected"] = "true" if selection_status.get((gene, sequence)) == "selected" else "false"
        if kind == "orf_candidate":
            row["Stage"] = "orf"
            row["Status"] = selection_status.get((gene, sequence), "not-reached")
            row["Outcome Detail"] = selection_details.get((gene, sequence), "")
    elif kind == "query_domain":
        domain = event["domain"]
        row["Row Type"] = "query_domain"
        row["Query ORF ID"] = sequence_ids[event["orf_sequence"]]
        row["Query ORF Selected"] = "true" if selection_status.get((gene, event["orf_sequence"])) == "selected" else "false"
        row.update({"Query Domain Name": domain.name, "Query Domain ID": domain.accession, "Query Domain Start 0-based": domain.start, "Query Domain End 0-based-exclusive": domain.end})
    elif kind == "ortholog":
        ortholog = event["ortholog"]
        row["Row Type"] = "ortholog"
        row["Query ORF ID"] = sequence_ids[event["orf_sequence"]]
        row.update(_ortholog_fields(ortholog))
    elif kind == "ortholog_retrieval":
        row["Row Type"] = "ortholog_retrieval"
        row["Stage"] = "blast"
        if event.get("orf_sequence"):
            row["Query ORF ID"] = sequence_ids[event["orf_sequence"]]
        row["Attempted Hit Identifier"] = event.get("attempted_identifier", "")
        row["Ortholog Hit ID"] = event.get("hit_id") or ""
        row["Ortholog Accession"] = event.get("accession") or ""
    elif kind in {"blast_search", "blast_search_summary"}:
        row["Row Type"] = "stage_status"
        row["Stage"] = "blast"
        if event.get("orf_sequence"):
            row["Query ORF ID"] = sequence_ids[event["orf_sequence"]]
        details = []
        for key in ("hit_count", "retained_hit_count", "retrieved_count", "skipped_count"):
            if key in event:
                details.append(f"{key}={event[key]}")
        row["Outcome Detail"] = "; ".join(details)
    elif kind == "ortholog_domain":
        domain = event["domain"]
        row["Row Type"] = "ortholog_domain"
        row["Query ORF ID"] = sequence_ids[event["orf_sequence"]]
        row.update(_ortholog_fields(event["ortholog"]))
        row.update({"Ortholog Domain Name": domain.name, "Ortholog Domain ID": domain.accession, "Ortholog Domain Start 0-based": domain.start, "Ortholog Domain End 0-based-exclusive": domain.end})
    elif kind == "domain_comparison":
        qdomain, odomain = event["query_domain"], event["ortholog_domain"]
        row["Row Type"] = "domain_comparison"
        row["Query ORF ID"] = sequence_ids[event["orf_sequence"]]
        row.update(_ortholog_fields(event["ortholog"]))
        row.update({"Query Domain Name": qdomain.name, "Query Domain ID": qdomain.accession, "Query Domain Start 0-based": qdomain.start, "Query Domain End 0-based-exclusive": qdomain.end, "Ortholog Domain Name": odomain.name, "Ortholog Domain ID": odomain.accession, "Ortholog Domain Start 0-based": odomain.start, "Ortholog Domain End 0-based-exclusive": odomain.end, "Domain Pair Identity Ratio": event["domain_identity"], "Domain Pair Identity Percent": event["domain_identity"] * 100, "Threshold Ratio": event.get("threshold"), "Threshold Passed": "true" if event["status"] == "accepted" else "false"})
    elif kind == "motif":
        row["Row Type"] = "motif"
        row["Query ORF ID"] = sequence_ids[event["orf_sequence"]]
        row.update(_ortholog_fields(event["ortholog"]))
        qdomain, odomain = event["query_domain"], event["ortholog_domain"]
        row.update({"Query Domain Name": qdomain.name, "Query Domain ID": qdomain.accession, "Query Domain Start 0-based": qdomain.start, "Query Domain End 0-based-exclusive": qdomain.end, "Ortholog Domain Name": odomain.name, "Ortholog Domain ID": odomain.accession, "Ortholog Domain Start 0-based": odomain.start, "Ortholog Domain End 0-based-exclusive": odomain.end, "Domain Pair Identity Ratio": event["domain_identity"], "Domain Pair Identity Percent": event["domain_identity"] * 100, "Threshold Ratio": event.get("threshold"), "Threshold Passed": "true", "Motif Database": event["motif_database_name"]})
        motif = event["motif_record"]
        motif_values = dict(zip(motif.get_headers(), motif.get_values()))
        mapping = {"Matrix ID": "Matrix ID", "Motif Name": "Motif Name", "Motif Class": "Motif Class", "Motif Inference E-value": "Inference E-value", "Motif URL": "Motif URL"}
        for target, source in mapping.items():
            row[target] = motif_values.get(source, "")
        pfm = motif_values.get("PFM", "")
        row["PFM JSON"] = json.dumps(pfm, ensure_ascii=False, separators=(",", ":")) if pfm != "" else ""
    elif kind == "gene_outcome":
        row["Row Type"] = "stage_status"
    elif kind == "gene_failure":
        row["Row Type"] = "gene_error"
    elif kind == "stage_status":
        row["Row Type"] = "stage_status"
    else:
        return None
    return row


def _ortholog_fields(ortholog):
    identity, identity_percent = _ratio_and_percent(ortholog.identity)
    coverage, coverage_percent = _ratio_and_percent(ortholog.query_coverage)
    return {
        "Ortholog Hit ID": ortholog.hit_id or "", "Ortholog Accession": ortholog.accession,
        "Ortholog Description": ortholog.description, "Ortholog Species Name": ortholog.species_name,
        "Ortholog Species Tax ID": ortholog.species_tax_id, "BLAST E-value": ortholog.evalue,
        "Ortholog Identity Ratio": identity, "Ortholog Identity Percent": identity_percent,
        "Query Coverage Ratio": coverage, "Query Coverage Percent": coverage_percent,
    }


class _Bundle:
    def __init__(self, destination, metadata):
        self.path = destination
        self.metadata = metadata
        self.events = []
        self.sequence_by_id = {}
        self.sequence_ids = {}
        self._write_results()
        self._write_orfs()
        self.flush_metadata()

    def add_event(self, gene, event):
        self.events.append((gene, event))
        sequence = event.get("orf_sequence")
        if sequence:
            identifier = _orf_id(sequence)
            self.sequence_ids[sequence] = identifier
            self.sequence_by_id[identifier] = sequence
        self.metadata["current_gene"] = gene
        self.metadata["current_stage"] = event.get("stage", "unknown")
        self.metadata["updated_at"] = _timestamp()
        self._write_results()
        self._write_orfs()
        self.flush_metadata()

    def _write_results(self):
        selection = {}
        selection_details = {}
        for gene, event in self.events:
            if event["event"] == "orf_selection":
                selection[(gene, event["orf_sequence"])] = event["status"]
                selection_details[(gene, event["orf_sequence"])] = event.get("detail", "")
        selected_by_gene = {
            gene: any(status == "selected" for (event_gene, _), status in selection.items() if event_gene == gene)
            for gene, _ in self.events
        }
        for gene, event in self.events:
            if event["event"] == "orf_candidate" and (gene, event["orf_sequence"]) not in selection and selected_by_gene.get(gene):
                selection[(gene, event["orf_sequence"])] = "not-selected"
        path = self.path / "results.tsv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=RESULT_HEADERS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            for gene, event in self.events:
                row = _event_row(event, gene, self.sequence_ids, selection, selection_details)
                if row is not None:
                    writer.writerow(row)
            handle.flush()
            os.fsync(handle.fileno())

    def _write_orfs(self):
        with (self.path / "orfs.fasta").open("w", encoding="utf-8", newline="\n") as handle:
            for identifier in sorted(self.sequence_by_id):
                sequence = self.sequence_by_id[identifier]
                handle.write(f">{identifier}\n")
                for start in range(0, len(sequence), 80):
                    handle.write(sequence[start:start + 80] + "\n")
            handle.flush()
            os.fsync(handle.fileno())

    def flush_metadata(self):
        temporary = self.path / ".metadata.json.tmp"
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(self.metadata, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, self.path / "metadata.json")


def _metadata(args, inputs, species_path):
    try:
        package_version = importlib.metadata.version("lverage")
    except importlib.metadata.PackageNotFoundError:
        package_version = "unknown"
    return {
        "schema_version": 1,
        "started_at": _timestamp(), "updated_at": _timestamp(), "finished_at": None,
        "completion": "in-progress", "current_gene": None, "current_stage": "preflight-complete",
        "summary": {"genes_total": len(inputs), "genes_completed": 0, "genes_failed": 0, "genes_intentionally_stopped": 0},
        "errors": [],
        "tools": {
            "lverage": {"version": package_version},
            "blastp": {"configured_executable": args.blastp, "resolved_executable": shutil.which(args.blastp), "version": "unknown"},
            "blastdbcmd": {"configured_executable": args.blastdbcmd, "resolved_executable": shutil.which(args.blastdbcmd), "version": "unknown"},
            "pfamscan": {"configured_executable": args.pfamscan, "resolved_executable": shutil.which(args.pfamscan), "version": "unknown"},
            "orffinder": {"version": _distribution_version("orffinder")},
            "biopython": {"version": _distribution_version("biopython")},
            "jaspar": {"release": "JASPAR2024", "service_version": "unknown"},
        },
        "settings": {
            "identity_threshold_ratio": args.identity_threshold,
            "valid_pfam": args.valid_pfam or [], "jaspar_hits": args.jaspar_hits,
            "jaspar_evalue": args.jaspar_evalue, "stop_after": args.stop_after or "full-pipeline",
            "orf_minimum_length_nt": 75, "blast_evalue": 1e-6, "blast_top_n": 20,
        },
        "references": {
            "input_files": [{"gene": item["gene"], "path": item["path"], "sha256": item["sha256"]} for item in inputs],
            "species_map": {"path": str(species_path), "sha256": _sha256(species_path)},
            "blast_database_prefix": args.blast_db, "pfam_database_directory": args.pfam_db,
        },
        "sequence_id_convention": "orf_<full lowercase SHA-256 of ASCII protein sequence>; identical sequences share one FASTA record; table rows retain per-gene associations",
        "limitations": ["Evidence is flushed as pipeline events complete; this is not checkpoint/resume or arbitrary-crash recovery."],
    }


def _distribution_version(name):
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _stage_completion_events(final_status, last_stage):
    stages = ("orf", "query_pfam", "blast", "ortholog_pfam", "alignment", "motif")
    last_index = stages.index(last_stage) if last_stage in stages else -1
    events = []
    if final_status == "STOPPED_AFTER_ALIGNMENT":
        events.append({"event": "stage_status", "stage": "motif", "status": "not-reached", "detail": "intentional stop after alignment"})
        return events
    for stage in stages[last_index + 1:]:
        events.append({"event": "stage_status", "stage": stage, "status": "not-reached", "detail": "upstream stage did not produce required evidence"})
    return events


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        inputs, species_map, species_path = _validate_inputs(args)
    except (InputValidationError, OSError, ValueError) as error:
        print(f"lverage: invalid input/configuration: {_redact_error(error)}", file=sys.stderr)
        return 2

    output_dir = Path(args.output_dir)
    try:
        output_dir.mkdir()
    except OSError as error:
        print(f"lverage: could not create output directory: {_redact_error(error)}", file=sys.stderr)
        return 2

    metadata = _metadata(args, inputs, species_path)
    bundle = _Bundle(output_dir, metadata)
    failures = 0
    try:
        pipeline = _build_pipeline(args, species_map)
    except Exception as error:
        failures = len(inputs)
        safe_error = _redact_error(error)
        for item in inputs:
            gene = item["gene"]
            metadata["errors"].append({"gene": gene, "stage": "initialization", "error": safe_error})
            bundle.add_event(gene, {"event": "gene_failure", "stage": "initialization", "status": "failed", "error": safe_error})
        metadata["summary"]["genes_failed"] = failures
        metadata["completion"] = "failed"
        metadata["finished_at"] = _timestamp()
        bundle.flush_metadata()
        print(f"lverage: pipeline initialization failed: {safe_error}", file=sys.stderr)
        return 1

    from .evidence import EvidenceCollector
    for item in inputs:
        gene = item["gene"]
        bundle.add_event(gene, {"event": "gene_start", "stage": "input", "status": "in-progress", "error": ""})
        collector = EvidenceCollector(callback=lambda event, gene=gene: bundle.add_event(gene, event))
        try:
            pipeline.run(item["fragments"], evidence_collector=collector, stop_after=args.stop_after)
            status = pipeline.lverage_code.name
            for event in _stage_completion_events(status, collector.stage):
                bundle.add_event(gene, event)
            retrieval_failures = [
                event for event in collector.events
                if event.get("event") == "ortholog_retrieval"
                and event.get("status") == "skipped"
            ]
            if retrieval_failures:
                failures += 1
                metadata["summary"]["genes_failed"] = failures
                metadata["errors"].append({
                    "gene": gene, "stage": "blast",
                    "error": f"{len(retrieval_failures)} BLAST hit(s) could not be retrieved",
                })
                print(
                    f"lverage: {len(retrieval_failures)} BLAST hit(s) could not be retrieved for gene {gene}",
                    file=sys.stderr,
                )
            if status == "STOPPED_AFTER_ALIGNMENT":
                metadata["summary"]["genes_intentionally_stopped"] += 1
            elif not retrieval_failures:
                metadata["summary"]["genes_completed"] += 1
            if not any(event["event"] == "gene_outcome" for event in collector.events):
                bundle.add_event(gene, {"event": "gene_outcome", "stage": collector.stage, "status": status, "error": ""})
        except Exception as error:
            failures += 1
            safe_error = _redact_error(error)
            metadata["errors"].append({"gene": gene, "stage": collector.stage, "error": safe_error})
            metadata["summary"]["genes_failed"] = failures
            bundle.add_event(gene, {"event": "gene_failure", "stage": collector.stage, "status": "failed", "error": safe_error})
            stages = ("orf", "query_pfam", "blast", "ortholog_pfam", "alignment", "motif")
            last_index = stages.index(collector.stage) if collector.stage in stages else -1
            for stage in stages[last_index + 1:]:
                bundle.add_event(gene, {"event": "stage_status", "stage": stage, "status": "not-reached", "detail": "execution failed in upstream stage"})
            print(f"lverage: execution failed for gene {gene} at {collector.stage}: {safe_error}", file=sys.stderr)
        metadata["current_gene"] = gene
        metadata["current_stage"] = "complete"
        metadata["updated_at"] = _timestamp()
        bundle.flush_metadata()

    metadata["finished_at"] = _timestamp()
    metadata["completion"] = "failed" if failures else ("intentionally-stopped" if metadata["summary"]["genes_intentionally_stopped"] else "completed")
    metadata["current_gene"] = None
    metadata["current_stage"] = "finished"
    bundle.flush_metadata()
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
