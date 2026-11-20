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
    "Gene", "Ortholog Species", "BLAST Description", "Query Pfam Name",
    "Hit Pfam Name", "Global Alignment Identity Percent", "Motif Name",
    "JASPAR Tool DBD Identity Ratio", "Motif PFM JSON", "ORF ID",
]
ERROR_HEADERS = [
    "Gene", "ORF ID", "Stage", "BLAST Hit ID", "BLAST Accession",
    "Query Pfam Name", "Hit Pfam Name", "Attempted Identifier", "Error",
]


class InputValidationError(Exception):
    """A user-provided input or setting failed preflight validation."""


def _require_external_path(path, description):
    repository_root = Path(__file__).resolve().parents[2]
    resolved = Path(path).expanduser().resolve()
    if resolved == repository_root or repository_root in resolved.parents:
        raise InputValidationError(
            f"{description} must be outside the Lverage repository: {resolved}"
        )


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
            "results.csv, errors.csv, orfs.fasta, and metadata.json."
        ),
        epilog=(
            "results.csv has ten columns ordered by gene and supporting evidence. "
            "errors.csv holds execution diagnostics. Every row "
            "repeats its applicable supporting fields; blank cells mean "
            "unavailable or not reached. --stop-after alignment "
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
    parser.add_argument("--output-dir", required=True, metavar="DIR", help="new directory for results.csv, errors.csv, orfs.fasta, and metadata.json")
    parser.add_argument("--blastp", default="blastp", metavar="EXE", help="blastp executable name or path")
    parser.add_argument("--blastdbcmd", default="blastdbcmd", metavar="EXE", help="blastdbcmd executable name or path")
    parser.add_argument("--pfamscan", default="pfam_scan.pl", metavar="EXE", help="pfam_scan.pl executable name or path")
    parser.add_argument("--identity-threshold", type=_identity_ratio, default=0.7, metavar="RATIO", help="minimum domain identity ratio")
    parser.add_argument("--valid-pfam", action="append", default=None, metavar="ACCESSION", help="accepted Pfam accession; repeatable")
    parser.add_argument("--jaspar-hits", type=_positive_int, default=10, metavar="N", help="maximum JASPAR motif hits")
    parser.add_argument("--jaspar-evalue", type=_positive_float, default=1e-6, metavar="VALUE", help="maximum JASPAR inference-tool BLAST E-value")
    parser.add_argument("--jaspar-data-dir", default=os.environ.get("LVERAGE_JASPAR_DATA_DIR"), metavar="DIR", help="external local JASPAR inference data directory (data only; no source checkout required)")
    parser.add_argument("--jaspar-sqlite-db", default=os.environ.get("LVERAGE_JASPAR_SQLITE_DB"), metavar="FILE", help="local JASPAR2026.sqlite matrix and annotation database")
    parser.add_argument("--hmmscan", default="hmmscan", metavar="EXE", help="hmmscan executable name or path for local JASPAR inference")
    parser.add_argument("--hmmalign", default="hmmalign", metavar="EXE", help="hmmalign executable name or path for local JASPAR inference")
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
    _require_external_path(pfam_dir, "PfamScan database")
    if not args.blast_db.strip():
        raise InputValidationError("--blast-db must be a nonempty database prefix")
    _require_external_path(args.blast_db, "protein BLAST database")
    if any(not executable.strip() for executable in (
            args.blastp, args.blastdbcmd, args.pfamscan,
            args.hmmscan, args.hmmalign)):
        raise InputValidationError("executable paths must not be empty")
    if args.valid_pfam is not None and any(not value.strip() for value in args.valid_pfam):
        raise InputValidationError("--valid-pfam values must not be empty")
    output_dir = Path(args.output_dir)
    if os.path.lexists(output_dir):
        raise InputValidationError(f"output directory already exists: {output_dir}")
    if not output_dir.parent.is_dir():
        raise InputValidationError(f"output parent directory does not exist: {output_dir.parent}")
    _require_external_path(output_dir, "output directory")
    return inputs, species_map, species_path


def _build_pipeline(args, species_map):
    from .blast import LocalBlastSearcher
    from .jaspar import Jaspar2026MotifDB
    from .orf_searcher import OrffinderOrfSearcher
    from .pfam import LocalPfamScanner
    from .pipeline import Lverage
    ortholog_searcher = LocalBlastSearcher(
        database_path=args.blast_db, species_map=species_map,
        blastp_path=args.blastp, blastdbcmd_path=args.blastdbcmd,
    )
    motif_database = Jaspar2026MotifDB(
        n_hits=args.jaspar_hits,
        escore_threshold=args.jaspar_evalue,
        data_dir=args.jaspar_data_dir,
        sqlite_db=args.jaspar_sqlite_db,
        blastp_path=args.blastp,
        hmmscan_path=args.hmmscan,
        hmmalign_path=args.hmmalign,
    )
    return Lverage(
        motif_database_list=[motif_database],
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


def _grouped_result_rows(gene, events, sequence_ids):
    """Order completed evidence without creating diagnostic-only result rows.

    Match domain occurrences by their original objects, rather than their names.
    Repeated motif records and distinct comparisons remain separate rows even
    when the selected export columns have identical values.
    """
    candidates = {}
    queries = {}
    hits = {}
    domains = {}
    comparisons = {}
    motifs = {}
    retrievals = {}

    def pair_key(event):
        return tuple(id(event.get(key)) for key in
                     ("ortholog", "query_domain", "ortholog_domain"))

    for event in events:
        sequence = event.get("orf_sequence")
        kind = event["event"]
        if sequence:
            candidates.setdefault(sequence, None)
        if kind == "query_domain":
            queries.setdefault(sequence, []).append(event)
        elif kind == "ortholog":
            hits.setdefault(sequence, []).append(event["ortholog"])
        elif kind == "ortholog_domain":
            domains.setdefault(id(event["ortholog"]), []).append(event["domain"])
        elif kind == "domain_comparison":
            comparisons[pair_key(event)] = event
        elif kind == "motif":
            motifs.setdefault(pair_key(event), []).append(event)
        elif kind == "ortholog_retrieval" and event.get("error"):
            key = (sequence, event.get("hit_group") or
                   event.get("hit_id") or event.get("accession"))
            retrievals[key] = event

    def row(sequence, query=None, hit=None, domain=None, comparison=None,
            motif=None, retrieval=None):
        result = {header: "" for header in RESULT_HEADERS}
        result.update({"Gene": gene, "ORF ID": sequence_ids.get(sequence, "")})
        if query is not None:
            result["Query Pfam Name"] = query.name
        if hit is not None:
            result["BLAST Description"] = hit.description
            result["Ortholog Species"] = hit.species_name
        if retrieval is not None:
            result["BLAST Description"] = retrieval.get("description", "")
            result["Ortholog Species"] = retrieval.get("species_name", "")
        if domain is not None:
            result["Hit Pfam Name"] = domain.name
        if comparison is not None:
            result["Global Alignment Identity Percent"] = comparison["domain_identity"] * 100
        if motif is not None:
            record = motif["motif_record"]
            values = dict(zip(record.get_headers(), record.get_values()))
            result["Motif Name"] = values.get("Motif Name", "")
            result["JASPAR Tool DBD Identity Ratio"] = values.get(
                "JASPAR Tool DBD Identity Ratio", "")
            pfm = values.get("PFM")
            if pfm is not None:
                result["Motif PFM JSON"] = json.dumps(
                    pfm, ensure_ascii=False, separators=(",", ":"))
        return result

    rows = []
    for sequence in sorted(candidates, key=len, reverse=True):
        query_events = queries.get(sequence, [])
        orthologs = hits.get(sequence, [])
        sequence_rows = []
        paired_domains = set()
        for query_event in query_events:
            query = query_event["domain"]
            query_rows = []
            if query_event.get("status") == "accepted":
                for hit in orthologs:
                    hit_domains = domains.get(id(hit), [])
                    matching = [domain for domain in hit_domains
                                if domain.accession.split(".", 1)[0] ==
                                query.accession.split(".", 1)[0]]
                    for domain in matching:
                        paired_domains.add((id(hit), id(domain)))
                        key = (id(hit), id(query), id(domain))
                        comparison = comparisons.get(key)
                        predictions = motifs.get(key, [])
                        for prediction in predictions or [None]:
                            query_rows.append(row(sequence, query, hit, domain,
                                                  comparison, prediction))
                    if not hit_domains:
                        query_rows.append(row(sequence, query, hit))
            sequence_rows.extend(query_rows or [row(sequence, query)])

        # A domain without a matching query is evidence of a hit-domain scan,
        # not a query/hit-domain comparison. Leave its query field blank.
        for hit in orthologs:
            hit_domains = domains.get(id(hit), [])
            for domain in hit_domains:
                if (id(hit), id(domain)) not in paired_domains:
                    sequence_rows.append(row(sequence, hit=hit, domain=domain))
            if not query_events and not hit_domains:
                sequence_rows.append(row(sequence, hit=hit))

        for (orf_sequence, _), event in retrievals.items():
            if orf_sequence != sequence:
                continue
            recovered = any(
                (event.get("hit_id") and event["hit_id"] == hit.hit_id)
                or (event.get("accession") and event["accession"] == hit.accession)
                for hit in orthologs
            )
            if not recovered:
                sequence_rows.append(row(sequence, retrieval=event))
        rows.extend(sequence_rows or [row(sequence)])
    rows = rows or [row(None)]
    # Rank all retained motifs for the gene, across supporting comparisons.
    # Stable ties preserve the existing evidence order; zero is a real score.
    def motif_rank(result):
        identity = result["JASPAR Tool DBD Identity Ratio"]
        return (identity == "", -float(identity) if identity != "" else 0)

    rows.sort(key=motif_rank)
    return rows


def _error_rows(events, sequence_ids):
    """Keep execution diagnostics out of the scientific results table."""
    for gene, event in events:
        if not event.get("error"):
            continue
        hit = event.get("ortholog")
        query_domain = event.get("query_domain")
        hit_domain = event.get("ortholog_domain")
        yield {
            "Gene": gene,
            "ORF ID": sequence_ids.get(event.get("orf_sequence"), ""),
            "Stage": event.get("stage_name") or event.get("stage", ""),
            "BLAST Hit ID": (hit.hit_id or "") if hit is not None else event.get("hit_id", ""),
            "BLAST Accession": hit.accession if hit is not None else event.get("accession", ""),
            "Query Pfam Name": query_domain.name if query_domain is not None else "",
            "Hit Pfam Name": hit_domain.name if hit_domain is not None else "",
            "Attempted Identifier": event.get("attempted_identifier", ""),
            "Error": _redact_error(event["error"]),
        }


class _Bundle:
    def __init__(self, destination, metadata):
        self.path = destination
        self.metadata = metadata
        self.events = []
        self.sequence_by_id = {}
        self.sequence_ids = {}
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
        self.flush_metadata()

    def finish(self):
        """Render ordered in-memory evidence and publish the final output files."""
        grouped_events = {}
        for gene, event in self.events:
            grouped_events.setdefault(gene, []).append(event)
        rows = []
        for gene, events in grouped_events.items():
            rows.extend(_grouped_result_rows(gene, events, self.sequence_ids))
        errors = list(_error_rows(self.events, self.sequence_ids))
        self._write_csv("results.csv", RESULT_HEADERS, rows)
        self._write_csv("errors.csv", ERROR_HEADERS, errors)
        self._write_orfs()
        self.flush_metadata()

    def _write_csv(self, filename, headers, rows):
        path = self.path / filename
        temporary = self.path / ("." + filename + ".tmp")
        try:
            with temporary.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=headers, lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, path)
        finally:
            if temporary.exists():
                temporary.unlink()

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
        "schema_version": 3,
        "results_format": "compact-csv",
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
            "jaspar": {
                "release": "JASPAR2026",
                "release_verification_method": "JASPAR2026.sqlite schema and exact per-taxonomy CORE inventory comparison",
                "inference_algorithm_provenance": {
                    "repository": "https://github.com/amathelier/JASPAR-inference-tool",
                    "source_revision": "ac4024148f1f140c7370d7b13a70556b4d13f7a3",
                    "source_version": "2026.1.0",
                    "runtime_checkout_required": False,
                    "runtime_python_required": False,
                },
                "data": {
                    "release": "JASPAR2026",
                    "directory": args.jaspar_data_dir,
                    "annotations_and_matrices_database": args.jaspar_sqlite_db,
                    "database_sha256": None,
                },
                "executables": {
                    "blastp": args.blastp,
                    "hmmscan": args.hmmscan,
                    "hmmalign": args.hmmalign,
                },
                "online_requests": False,
                "validation": "not-required" if args.stop_after == "alignment" else "pending",
            },
        },
        "settings": {
            "identity_threshold_ratio": args.identity_threshold,
            "valid_pfam": args.valid_pfam or [], "jaspar_hits": args.jaspar_hits,
            "jaspar_inference_blast_evalue_threshold": args.jaspar_evalue,
            "stop_after": args.stop_after or "full-pipeline",
            "orf_minimum_length_nt": 75, "blast_evalue": 1e-6, "blast_top_n": 20, "blast_seg": "default",
            "blast_description_excluded_terms": None,
        },
        "references": {
            "input_files": [{"gene": item["gene"], "path": item["path"], "sha256": item["sha256"]} for item in inputs],
            "species_map": {"path": str(species_path), "sha256": _sha256(species_path)},
            "blast_database_prefix": args.blast_db, "pfam_database_directory": args.pfam_db,
        },
        "sequence_id_convention": "orf_<full lowercase SHA-256 of ASCII protein sequence>; identical sequences share one FASTA record; table rows retain per-gene associations",
        "result_columns": RESULT_HEADERS.copy(),
        "result_row_order": "gene filename ascending; JASPAR Tool DBD Identity Ratio descending within each gene; blank ratios last; stable ties preserve evidence order",
        "error_columns": ERROR_HEADERS.copy(),
        "limitations": [
            "Results and errors are collected in memory and written after processing, including handled gene failures. An abrupt process termination may leave no final CSV or FASTA.",
            "The compact CSV omits matrix IDs and domain coordinates; distinct versions and occurrences can produce identical-looking rows. Rows are not deduplicated.",
        ],
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
        jaspar_database = pipeline.motif_database_list[0]
        metadata_factory = getattr(jaspar_database, "configuration_metadata", None)
        metadata["tools"]["jaspar"] = (
            metadata_factory() if callable(metadata_factory)
            else {"release": getattr(jaspar_database, "name", "unknown")}
        )
        metadata["tools"]["jaspar"]["validation"] = (
            "not-required" if args.stop_after == "alignment" else "not-run"
        )
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
        bundle.finish()
        print(f"lverage: pipeline initialization failed: {safe_error}", file=sys.stderr)
        return 1

    excluded_terms = getattr(pipeline.ortholog_searcher, "excluded_terms", None)
    metadata["settings"]["blast_description_excluded_terms"] = (
        list(excluded_terms) if excluded_terms is not None else None
    )

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
            failure_event = {
                "event": "gene_failure", "stage": collector.stage,
                "status": "failed", "error": safe_error,
            }
            attempt = next((event for event in reversed(collector.events)
                            if event.get("event") == "stage_attempt"
                            and event.get("stage_name") == collector.stage), None)
            if attempt is not None:
                for key in (
                    "orf_sequence", "ortholog", "query_domain",
                    "ortholog_domain", "motif_database_name",
                ):
                    if key in attempt:
                        failure_event[key] = attempt[key]
            bundle.add_event(gene, failure_event)
            stages = ("orf", "query_pfam", "blast", "ortholog_pfam", "alignment", "motif")
            last_index = stages.index(collector.stage) if collector.stage in stages else -1
            for stage in stages[last_index + 1:]:
                bundle.add_event(gene, {"event": "stage_status", "stage": stage, "status": "not-reached", "detail": "execution failed in upstream stage"})
            print(f"lverage: execution failed for gene {gene} at {collector.stage}: {safe_error}", file=sys.stderr)
        metadata_factory = getattr(jaspar_database, "configuration_metadata", None)
        if callable(metadata_factory):
            metadata["tools"]["jaspar"] = metadata_factory()
        motif_attempted = any(
            event.get("event") == "stage_attempt"
            and event.get("stage_name") == "motif"
            for event in collector.events
        )
        if args.stop_after == "alignment":
            metadata["tools"]["jaspar"]["validation"] = "not-required"
        elif getattr(jaspar_database, "_validated", False):
            metadata["tools"]["jaspar"]["validation"] = "verified"
        elif motif_attempted:
            metadata["tools"]["jaspar"]["validation"] = "failed"
        metadata["current_gene"] = gene
        metadata["current_stage"] = "complete"
        metadata["updated_at"] = _timestamp()
        bundle.flush_metadata()

    metadata["finished_at"] = _timestamp()
    metadata["completion"] = "failed" if failures else ("intentionally-stopped" if metadata["summary"]["genes_intentionally_stopped"] else "completed")
    metadata["current_gene"] = None
    metadata["current_stage"] = "finished"
    bundle.finish()
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
