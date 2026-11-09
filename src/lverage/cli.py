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
    "Gene", "Gene Status", "Gene Error", "Gene Outcome",
    "ORF Candidate Order", "ORF ID", "ORF Length (aa)",
    "ORF Examination Status", "ORF Selection Status", "ORF Outcome",
    "Query Pfam Occurrence", "Query Pfam ID", "Query Pfam Name",
    "Query Pfam Start (0-based aa)", "Query Pfam End (0-based-exclusive aa)",
    "Query Pfam Status", "Query Pfam Error", "Query Pfam Outcome",
    "BLAST Hit Group", "BLAST Hit ID", "BLAST Accession", "BLAST Description",
    "BLAST Species Name", "BLAST Species Tax ID", "BLAST E-value",
    "BLAST Identity Ratio", "BLAST Identity Percent", "Query Coverage Ratio",
    "Query Coverage Percent", "BLAST Status", "BLAST Error", "BLAST Outcome",
    "BLAST Retrieval Attempted Identifier",
    "Hit Pfam Occurrence", "Hit Pfam ID", "Hit Pfam Name",
    "Hit Pfam Start (0-based aa)", "Hit Pfam End (0-based-exclusive aa)",
    "Hit Pfam Status", "Hit Pfam Error", "Hit Pfam Outcome",
    "Global Alignment ID", "Global Alignment Identity Ratio",
    "Global Alignment Identity Percent", "Global Alignment Threshold Ratio",
    "Global Alignment Threshold Passed", "Global Alignment Status",
    "Global Alignment Error",
    "Motif Database", "Motif Matrix ID", "Motif Name", "Motif Class",
    "JASPAR Tool TF Name", "JASPAR Tool BLAST E-value",
    "JASPAR Tool Query Region", "JASPAR Tool TF Region",
    "JASPAR Tool DBD Identity Ratio", "Motif URL", "Motif PFM JSON",
    "Motif Status", "Motif Error", "Motif Outcome",
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
            "a grouped results.csv, orfs.fasta, and metadata.json evidence bundle."
        ),
        epilog=(
            "results.csv groups evidence by gene, ORF candidate, query Pfam, "
            "BLAST hit, hit Pfam, global alignment, and motif. Every row "
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
    parser.add_argument("--output-dir", required=True, metavar="DIR", help="new directory for results.csv, orfs.fasta, and metadata.json")
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


def _ratio_and_percent(value):
    return value, value * 100 if value is not None else None


def _ortholog_fields(ortholog):
    identity, identity_percent = _ratio_and_percent(ortholog.identity)
    coverage, coverage_percent = _ratio_and_percent(ortholog.query_coverage)
    return {
        "BLAST Hit ID": ortholog.hit_id or "", "BLAST Accession": ortholog.accession,
        "BLAST Description": ortholog.description, "BLAST Species Name": ortholog.species_name,
        "BLAST Species Tax ID": ortholog.species_tax_id, "BLAST E-value": ortholog.evalue,
        "BLAST Identity Ratio": identity, "BLAST Identity Percent": identity_percent,
        "Query Coverage Ratio": coverage, "Query Coverage Percent": coverage_percent,
    }


def _grouped_result_rows(gene, events, sequence_ids):
    """Render one gene as contiguous, self-contained result rows.

    Domain associations use the original objects retained by EvidenceCollector,
    not names or accessions: two occurrences can have identical annotations.
    Parent evidence is copied onto each descendant row; associations are built
    from the event hierarchy, never inferred by filling blank CSV cells.
    """
    def blank():
        return {header: "" for header in RESULT_HEADERS}

    def domain_fields(prefix, event, occurrence):
        domain = event["domain"]
        return {
            f"{prefix} Occurrence": occurrence,
            f"{prefix} ID": domain.accession,
            f"{prefix} Name": domain.name,
            f"{prefix} Start (0-based aa)": domain.start,
            f"{prefix} End (0-based-exclusive aa)": domain.end,
            f"{prefix} Status": event.get("status", ""),
        }

    def prepend(fields, children):
        # Repeat supporting fields on each result row without changing its
        # event-specific association.
        if not children:
            children = [blank()]
        for child in children:
            child.update(fields)
        return children

    def pair_key(event):
        return tuple(id(event.get(key)) for key in
                     ("ortholog", "query_domain", "ortholog_domain"))

    def motif_fields(event):
        values = dict(zip(event["motif_record"].get_headers(),
                          event["motif_record"].get_values()))
        row = blank()
        row.update({
            "Motif Database": event["motif_database_name"],
            "Motif Matrix ID": values.get("Matrix ID", ""),
            "Motif Name": values.get("Motif Name", ""),
            "Motif Class": values.get("Motif Class", ""),
            "JASPAR Tool TF Name": values.get("JASPAR Tool TF Name", ""),
            "JASPAR Tool BLAST E-value": values.get("JASPAR Tool BLAST E-value", ""),
            "JASPAR Tool Query Region": values.get("JASPAR Tool Query Region", ""),
            "JASPAR Tool TF Region": values.get("JASPAR Tool TF Region", ""),
            "JASPAR Tool DBD Identity Ratio": values.get("JASPAR Tool DBD Identity Ratio", ""),
            "Motif URL": values.get("Motif URL", ""),
            "Motif Status": event.get("status", "retained"),
        })
        pfm = values.get("PFM", "")
        if pfm != "":
            row["Motif PFM JSON"] = json.dumps(
                pfm, ensure_ascii=False, separators=(",", ":")
            )
        return row

    def diagnostic(prefix, event):
        row = blank()
        row[f"{prefix} Status"] = event.get("status", "failed")
        row[f"{prefix} Error"] = event.get("error", "")
        if f"{prefix} Outcome" in row:
            row[f"{prefix} Outcome"] = event.get("detail", "")
        return row

    def blast_outcome(summary):
        return "; ".join(
            [summary.get("status", "")] +
            [f"{name}={summary[name]}" for name in
             ("hit_count", "retrieved_count", "skipped_count") if name in summary]
        )

    failures = [event for event in events if event["event"] == "gene_failure"]
    retrieval_failures = [event for event in events
                          if event["event"] == "ortholog_retrieval"
                          and event.get("status") == "skipped"]
    handled_failures = set()
    outcomes = [event for event in events if event["event"] == "gene_outcome"]
    final = outcomes[-1] if outcomes else {}
    gene_fields = {
        "Gene": gene,
        "Gene Status": "failed" if failures or any(
            event["event"] == "ortholog_retrieval" and event.get("status") == "skipped"
            for event in events) else final.get("status", "in-progress"),
        "Gene Error": (
            failures[-1].get("error", "") if failures else
            "; ".join(event.get("error", "") for event in retrieval_failures)
        ),
        "Gene Outcome": (
            f"{failures[-1].get('stage', 'unknown')}: execution failed" if failures
            else "blast: partial-retrieval-failure" if retrieval_failures
            else f"{final.get('stage', 'pipeline')}: {final.get('detail') or final.get('status', '')}"
            if final else ""
        ),
    }
    candidates = {}
    for event in events:
        if event["event"] == "orf_candidate":
            candidates.setdefault(event["orf_sequence"], event)
    ordered_candidates = sorted(candidates, key=lambda sequence: -len(sequence))
    selections = {event["orf_sequence"]: event for event in events
                  if event["event"] == "orf_selection"}
    comparisons = {}
    for event in events:
        if event["event"] == "domain_comparison":
            comparisons.setdefault(pair_key(event), []).append(event)
    alignment_ids = {id(event): f"A{index}" for index, event in enumerate(
        (event for event in events if event["event"] == "domain_comparison"), 1)}

    def comparison_rows(comparison):
        key = pair_key(comparison)
        children = []
        # Preserve completion order, including a database failure after motifs.
        for event in events:
            if pair_key(event) != key:
                continue
            kind = event["event"]
            if kind == "motif":
                children.append(motif_fields(event))
            elif (kind == "stage_result" and event.get("stage_name") == "motif"
                  and event.get("motif_count") == 0):
                row = diagnostic("Motif", event)
                row["Motif Database"] = event.get("motif_database_name", "")
                row["Motif Outcome"] = "species unavailable" if event.get("status") == "species-unavailable" else "motif_count=0"
                children.append(row)
            elif kind == "gene_failure" and event.get("stage") == "motif":
                row = diagnostic("Motif", event)
                row["Motif Database"] = event.get("motif_database_name", "")
                children.append(row)
                handled_failures.add(id(event))
        return prepend({
            "Global Alignment ID": alignment_ids[id(comparison)],
            "Global Alignment Identity Ratio": comparison["domain_identity"],
            "Global Alignment Identity Percent": comparison["domain_identity"] * 100,
            "Global Alignment Threshold Ratio": comparison.get("threshold", ""),
            "Global Alignment Threshold Passed": str(comparison.get("status") == "accepted").lower(),
            "Global Alignment Status": comparison.get("status", ""),
        }, children)

    rows = []
    for order, sequence in enumerate(ordered_candidates, 1):
        selection = selections.get(sequence, {})
        selected = selection.get("status") == "selected"
        queries = [event for event in events if event["event"] == "query_domain"
                   and event.get("orf_sequence") == sequence]
        examined = bool(selection or queries) or any(
            event["event"] == "stage_attempt" and event.get("stage_name") == "query_pfam"
            and event.get("orf_sequence") == sequence for event in events)
        query_failures = [event for event in failures if event.get("stage") == "query_pfam"
                          and event.get("orf_sequence") == sequence]
        candidate_fields = {
            "ORF Candidate Order": order,
            "ORF ID": sequence_ids.get(sequence, _orf_id(sequence)),
            "ORF Length (aa)": len(sequence),
        }
        if examined:
            candidate_fields["ORF Examination Status"] = "examined"
            candidate_fields["ORF Selection Status"] = (
                "failed" if any(event.get("orf_sequence") == sequence for event in query_failures)
                else selection.get("status", "not-selected")
            )
        if selection.get("detail"):
            candidate_fields["ORF Outcome"] = selection["detail"]
        children = []
        if query_failures:
            candidate_fields["ORF Examination Status"] = "failed"
            for event in query_failures:
                children.append(diagnostic("Query Pfam", event))
                handled_failures.add(id(event))
        if not queries and selection.get("status") == "biological-no-hit":
            row = blank()
            row["Query Pfam Status"] = "no accepted domain"
            row["Query Pfam Outcome"] = selection.get("detail", "no accepted query domain")
            children.append(row)

        hit_events = [event for event in events if event["event"] == "ortholog"
                      and event.get("orf_sequence") == sequence]
        retrievals = [event for event in events if event["event"] == "ortholog_retrieval"
                      and event.get("orf_sequence") == sequence]
        blast_summaries = [event for event in events if event["event"] in
                           {"blast_search", "blast_search_summary"}
                           and event.get("orf_sequence") == sequence]
        accepted = [event for event in queries if event.get("status") == "accepted"]
        first_query = accepted[0] if accepted else None

        for q_index, query in enumerate(queries, 1):
            qdomain = query["domain"]
            query_children = []
            if selected and query.get("status") == "accepted":
                for h_index, hit_event in enumerate(hit_events, 1):
                    ortholog = hit_event["ortholog"]
                    domains = [event for event in events if event["event"] == "ortholog_domain"
                               and event.get("ortholog") is ortholog]
                    hit_children = []
                    for d_index, domain in enumerate(domains, 1):
                        hdomain = domain["domain"]
                        key = (id(ortholog), id(qdomain), id(hdomain))
                        matching = comparisons.get(key, [])
                        alignment_failures = [event for event in failures
                                              if event.get("stage") == "alignment"
                                              and pair_key(event) == key]
                        matched_any_query = any(
                            other["domain"].accession.split(".", 1)[0] == hdomain.accession.split(".", 1)[0]
                            for other in accepted)
                        # Display unpaired hit domains once, under the first query.
                        # Matched but not yet evaluated pairs stay visible with blank results.
                        same_family = qdomain.accession.split(".", 1)[0] == hdomain.accession.split(".", 1)[0]
                        if not same_family and (matched_any_query or query is not first_query):
                            continue
                        domain_children = []
                        for comparison in matching:
                            domain_children.extend(comparison_rows(comparison))
                        for event in alignment_failures:
                            domain_children.append(diagnostic("Global Alignment", event))
                            handled_failures.add(id(event))
                        fields = domain_fields("Hit Pfam", domain, domain.get("domain_index", d_index))
                        if not matched_any_query:
                            fields["Hit Pfam Outcome"] = "no matching query Pfam"
                        hit_children.extend(prepend(fields, domain_children))
                    scans = [event for event in events if event["event"] == "stage_result"
                             and event.get("stage_name") == "ortholog_pfam"
                             and event.get("ortholog") is ortholog]
                    scan_failures = [event for event in failures
                                     if event.get("stage") == "ortholog_pfam"
                                     and event.get("ortholog") is ortholog]
                    if query is first_query:
                        for event in scan_failures:
                            hit_children.append(diagnostic("Hit Pfam", event))
                            handled_failures.add(id(event))
                        if not domains and scans:
                            row = diagnostic("Hit Pfam", scans[-1])
                            row["Hit Pfam Outcome"] = f"domain_count={scans[-1].get('domain_count', 0)}"
                            hit_children.append(row)
                    # Each hit is visible even before it is scanned. On later
                    # query parents, only display hits with relevant child evidence.
                    if query is first_query or hit_children:
                        fields = _ortholog_fields(ortholog)
                        summary = blast_summaries[-1] if blast_summaries else {}
                        fields["BLAST Outcome"] = blast_outcome(summary)
                        fields.update({"BLAST Hit Group": h_index,
                                       "BLAST Status": hit_event.get("status", "retained")})
                        query_children.extend(prepend(fields, hit_children))
                if query is first_query:
                    # Retrieval attempt diagnostics stay beside their native hit.
                    for event in retrievals:
                        row = diagnostic("BLAST", event)
                        if blast_summaries:
                            row["BLAST Outcome"] = blast_outcome(blast_summaries[-1])
                        row.update({"BLAST Hit ID": event.get("hit_id", ""),
                                    "BLAST Accession": event.get("accession", ""),
                                    "BLAST Hit Group": event.get("hit_group", ""),
                                    "BLAST Description": event.get("description", ""),
                                    "BLAST Species Name": event.get("species_name", ""),
                                    "BLAST Species Tax ID": event.get("species_tax_id", ""),
                                    "BLAST E-value": event.get("evalue", ""),
                                    "BLAST Identity Ratio": event.get("identity_ratio", ""),
                                    "BLAST Identity Percent": event.get("identity_percent", ""),
                                    "Query Coverage Ratio": event.get("query_coverage_ratio", ""),
                                    "Query Coverage Percent": event.get("query_coverage_percent", ""),
                                    "BLAST Retrieval Attempted Identifier": event.get("attempted_identifier", "")})
                        query_children.append(row)
                    blast_failures = [event for event in failures if event.get("stage") == "blast"]
                    for event in blast_failures:
                        query_children.append(diagnostic("BLAST", event))
                        handled_failures.add(id(event))
                    if not hit_events and not retrievals and not blast_failures:
                        if blast_summaries:
                            summary = blast_summaries[-1]
                            row = diagnostic("BLAST", summary)
                            row["BLAST Outcome"] = blast_outcome(summary)
                        elif final.get("stage") == "blast":
                            row = diagnostic("BLAST", final)
                        else:
                            row = None
                        if row:
                            query_children.append(row)
            children.extend(prepend(domain_fields("Query Pfam", query, q_index), query_children))
        rows.extend(prepend(candidate_fields, children))

    # A failure without sufficiently specific evidence must still be visible.
    prefixes = {"query_pfam": "Query Pfam", "blast": "BLAST",
                "ortholog_pfam": "Hit Pfam", "alignment": "Global Alignment", "motif": "Motif"}
    for event in failures:
        if id(event) not in handled_failures and event.get("stage") in prefixes:
            row = diagnostic(prefixes[event["stage"]], event)
            row["Gene"] = gene  # Restart context: do not attach to the last ORF.
            row["Gene Outcome"] = f"unassociated failure at {event['stage']}"
            rows.append(row)
    return prepend(gene_fields, rows)

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
        grouped_events = {}
        gene_order = []
        for gene, event in self.events:
            if gene not in grouped_events:
                grouped_events[gene] = []
                gene_order.append(gene)
            grouped_events[gene].append(event)
        path = self.path / "results.csv"
        temporary = self.path / ".results.csv.tmp"
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=RESULT_HEADERS, delimiter=",", lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            for gene in gene_order:
                for row in _grouped_result_rows(
                    gene, grouped_events[gene], self.sequence_ids
                ):
                    writer.writerow(row)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)

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
        "schema_version": 2,
        "results_format": "grouped-csv",
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
        bundle.flush_metadata()
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
    bundle.flush_metadata()
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
