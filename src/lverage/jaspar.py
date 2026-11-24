"""Local JASPAR 2026 motif database adapter for Lverage."""

from contextlib import closing
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sqlite3
import tempfile

from .jaspar_inference import JASPAR_TAXA, infer_sequence, load_jaspar_data
from .jaspar_models import load_cisbp_models
from .motif_database import MotifDBRecordTemplate, MotifDBTemplate, MotifSearchRequest


# Algorithm provenance only. Runtime does not import from or inspect this checkout.
UPSTREAM_REPOSITORY = "https://github.com/amathelier/JASPAR-inference-tool"
UPSTREAM_REVISION = "ac4024148f1f140c7370d7b13a70556b4d13f7a3"
UPSTREAM_VERSION = "2026.1.0"
JASPAR_RELEASE = "JASPAR2026"
MATRIX_ID_PATTERN = re.compile(r"^(MA\d{4,5})\.(\d+)$")


@dataclass
class JasparRecord(MotifDBRecordTemplate):
    """One exact-version JASPAR profile and its inference evidence.

    Parameters
    ----------
    matrix_id : str
        Exact versioned profile identifier, for example ``MA0593.1``.
    motif_name : str
        Profile name from the local matrix database.
    pfm : dict
        Position frequency matrix with A, C, G and T count rows.
    motif_url : str
        Convenience link to the hosted profile; not used for inference.
    motif_class : str
        Transcription factor class annotated for the profile.
    inference_evalue : float
        BLAST E-value reported by local motif inference.
    inference_tf_name : str, optional
        Reference transcription factor name reported by inference.
    inference_query_region : str, optional
        BLAST query endpoints, one-based and inclusive, in the submitted
        ortholog protein or its window.
    inference_tf_region : str, optional
        One-based inclusive BLAST endpoints in the reference protein.
    inference_dbd_identity : float, optional
        Inference-tool domain identity ratio, distinct from Lverage's
        query-domain versus ortholog-domain alignment identity.
    """

    matrix_id: str
    motif_name: str
    pfm: dict
    motif_url: str
    motif_class: str
    inference_evalue: float
    inference_tf_name: str = ""
    inference_query_region: str = ""
    inference_tf_region: str = ""
    inference_dbd_identity: float = 0.0

    headers = (
        "Matrix ID", "Motif Name", "PFM", "Motif URL", "Motif Class",
        "JASPAR Tool TF Name", "JASPAR Tool BLAST E-value",
        "JASPAR Tool Query Region", "JASPAR Tool TF Region",
        "JASPAR Tool DBD Identity Ratio",
    )

    def get_values(self) -> list:
        """Return the record fields in export-header order.

        Returns
        -------
        list
            Field values in the order returned by ``get_headers()``.
        """
        return [
            self.matrix_id, self.motif_name, self.pfm, self.motif_url,
            self.motif_class, self.inference_tf_name, self.inference_evalue,
            self.inference_query_region, self.inference_tf_region,
            self.inference_dbd_identity,
        ]


class Jaspar2026MotifDB(MotifDBTemplate):
    """Run the adapted upstream algorithm against local JASPAR 2026 data.

    Only ``data_dir`` (the fork-generated ``files`` tree), ``sqlite_db``, and
    external BLASTP/HMMER executable locations are configurable. The upstream
    repository revision below records algorithm provenance, not a runtime
    checkout or interpreter requirement.

    Parameters
    ----------
    n_hits : int, optional
        Maximum retained predictions per request; default is 10.
    escore_threshold : float, optional
        Maximum inference BLAST E-value; default is ``1e-6``. Equivalence
        to the former REST endpoint's E-value has not been established.
    data_dir : str or pathlib.Path, optional
        Local inference reference-data tree. Defaults to the
        ``LVERAGE_JASPAR_DATA_DIR`` environment variable.
    sqlite_db : str or pathlib.Path, optional
        Local ``JASPAR2026.sqlite`` matrix and annotation database.
        Defaults to ``LVERAGE_JASPAR_SQLITE_DB``.
    blastp_path : str, optional
        BLASTP executable name or path; default is ``blastp``.
    hmmscan_path : str, optional
        HMMER domain-scan executable name or path; default is ``hmmscan``.
    hmmalign_path : str, optional
        HMMER alignment executable name or path; default is ``hmmalign``.
    """

    name = "JASPAR2026-local"
    metainfo = {
        "version": JASPAR_RELEASE,
        "collection": "JASPAR CORE",
        "inference_source": UPSTREAM_REPOSITORY,
        "inference_source_revision": UPSTREAM_REVISION,
        "inference_source_version": UPSTREAM_VERSION,
        "online_requests": False,
    }

    def __init__(self, n_hits=10, escore_threshold=1e-6, *, data_dir=None,
                 sqlite_db=None, blastp_path="blastp", hmmscan_path="hmmscan",
                 hmmalign_path="hmmalign"):
        if isinstance(n_hits, bool) or not isinstance(n_hits, int) or n_hits <= 0:
            raise ValueError("n_hits must be a positive integer")
        if (isinstance(escore_threshold, bool)
                or not isinstance(escore_threshold, (int, float))
                or not math.isfinite(escore_threshold) or escore_threshold <= 0):
            raise ValueError("escore_threshold must be a positive finite number")

        self.n_hits = n_hits
        self.escore_threshold = float(escore_threshold)
        self.data_dir = self._configured_path(data_dir, "LVERAGE_JASPAR_DATA_DIR")
        self.sqlite_db = self._configured_path(sqlite_db, "LVERAGE_JASPAR_SQLITE_DB")
        self.blastp_path = blastp_path
        self.hmmscan_path = hmmscan_path
        self.hmmalign_path = hmmalign_path
        self._validated = False
        self._species_ids = None
        self._inference_data = None
        self._database_sha256_value = None

    @staticmethod
    def _configured_path(value, environment_name):
        configured = value if value is not None else os.environ.get(environment_name)
        if configured == "":
            raise ValueError(f"{environment_name} must not be empty")
        return Path(configured).expanduser() if configured is not None else None

    def configuration_metadata(self):
        """Return code provenance separately from the required data release.

        Returns
        -------
        dict
            Source revision, data locations and release, configured executable
            paths, and database checksum once configuration is validated.
        """
        return {
            "release": JASPAR_RELEASE,
            "release_verification_method": (
                "JASPAR2026.sqlite schema and exact per-taxonomy CORE inventory comparison"
            ),
            "inference_algorithm_provenance": {
                "repository": UPSTREAM_REPOSITORY,
                "source_revision": UPSTREAM_REVISION,
                "source_version": UPSTREAM_VERSION,
                "runtime_checkout_required": False,
                "runtime_python_required": False,
            },
            "data": {
                "release": JASPAR_RELEASE,
                "directory": str(self.data_dir) if self.data_dir else None,
                "annotations_and_matrices_database": (
                    str(self.sqlite_db) if self.sqlite_db else None
                ),
                "database_sha256": self._database_sha256_value,
            },
            "executables": {
                "blastp": self.blastp_path,
                "hmmscan": self.hmmscan_path,
                "hmmalign": self.hmmalign_path,
            },
            "online_requests": False,
        }

    def _validate_configuration(self):
        if self._validated:
            return
        missing = [
            name for name, value in (
                ("--jaspar-data-dir / LVERAGE_JASPAR_DATA_DIR", self.data_dir),
                ("--jaspar-sqlite-db / LVERAGE_JASPAR_SQLITE_DB", self.sqlite_db),
            ) if value is None
        ]
        if missing:
            raise RuntimeError(
                "local JASPAR inference is not configured; provide " + ", ".join(missing)
            )
        for name, executable in (
            ("blastp", self.blastp_path),
            ("hmmscan", self.hmmscan_path),
            ("hmmalign", self.hmmalign_path),
        ):
            if not isinstance(executable, str) or not executable.strip():
                raise RuntimeError(f"configured {name} executable must not be empty")

        data_dir = self.data_dir.resolve()
        sqlite_db = self.sqlite_db.resolve()
        repository_root = Path(__file__).resolve().parents[2]
        for label, path in (("JASPAR inference data", data_dir),
                            ("JASPAR SQLite database", sqlite_db)):
            if path == repository_root or repository_root in path.parents:
                raise RuntimeError(f"{label} must be outside the Lverage repository: {path}")
        if not data_dir.is_dir():
            raise RuntimeError(f"JASPAR inference data directory does not exist: {data_dir}")
        if sqlite_db.name != "JASPAR2026.sqlite" or not sqlite_db.is_file():
            raise RuntimeError(
                "JASPAR annotation database must be a local JASPAR2026.sqlite file"
            )

        self._verify_inference_data(data_dir)
        self._verify_annotation_database(sqlite_db, data_dir)
        self.data_dir, self.sqlite_db = data_dir, sqlite_db
        self._database_sha256_value = self._database_sha256()
        self._validated = True

    @staticmethod
    def _verify_inference_data(data_dir):
        hmm_prefix = data_dir / "pfam" / "All.hmm"
        required = [data_dir / "cisbp", hmm_prefix]
        required.extend(Path(str(hmm_prefix) + suffix)
                        for suffix in (".h3f", ".h3i", ".h3m", ".h3p"))
        for taxon in JASPAR_TAXA:
            required.extend(data_dir / f"{taxon}.{suffix}.json"
                            for suffix in ("pfam", "profiles", "uniprot"))
            fasta_prefix = data_dir / f"{taxon}.fa"
            required.append(fasta_prefix)
            required.extend(Path(str(fasta_prefix) + suffix)
                            for suffix in (".phr", ".pin", ".psq"))
        missing = [str(path) for path in required if not path.exists()]
        if not (data_dir / "cisbp").is_dir():
            missing.append(str(data_dir / "cisbp/*.json"))
        elif not any((data_dir / "cisbp").glob("*.json")):
            missing.append(str(data_dir / "cisbp/*.json"))
        if missing:
            raise RuntimeError(
                "JASPAR inference data directory is incomplete; missing: "
                + ", ".join(missing[:12])
            )

    @staticmethod
    def _readonly_connection(path):
        return sqlite3.connect(f"{path.as_uri()}?mode=ro", uri=True)

    def _verify_annotation_database(self, path, data_dir):
        try:
            with closing(self._readonly_connection(path)) as connection:
                tables = {row[0] for row in connection.execute(
                    "SELECT name FROM sqlite_master WHERE type='table'"
                )}
                required = {"MATRIX", "MATRIX_DATA", "MATRIX_SPECIES", "MATRIX_ANNOTATION"}
                if not required.issubset(tables):
                    raise RuntimeError(
                        "JASPAR 2026 SQLite database lacks required motif tables: "
                        + ", ".join(sorted(required - tables))
                    )
                if connection.execute("SELECT 1 FROM MATRIX LIMIT 1").fetchone() is None:
                    raise RuntimeError("JASPAR 2026 SQLite database contains no matrices")
        except sqlite3.Error as error:
            raise RuntimeError(f"cannot read local JASPAR 2026 SQLite database: {error}") from error

        for taxon in JASPAR_TAXA:
            profile_path = data_dir / f"{taxon}.profiles.json"
            try:
                with profile_path.open(encoding="utf-8") as handle:
                    profiles = json.load(handle)
            except (OSError, UnicodeError, ValueError) as error:
                raise RuntimeError(f"cannot read local profile annotations {profile_path}: {error}") from error
            if not isinstance(profiles, dict) or not profiles:
                raise RuntimeError(f"local JASPAR profile annotations are empty or malformed: {profile_path}")
            with closing(self._readonly_connection(path)) as connection:
                rows = connection.execute(
                    "SELECT m.BASE_ID, m.VERSION, m.NAME FROM MATRIX m "
                    "JOIN MATRIX_ANNOTATION a ON a.ID = m.ID "
                    "WHERE upper(m.COLLECTION) = 'CORE' "
                    "AND a.TAG = 'tax_group' AND lower(a.VAL) = ?",
                    (taxon,),
                )
                database_profiles = {
                    f"{base_id}.{version}": name
                    for base_id, version, name in rows
                }
            if profiles != database_profiles:
                raise RuntimeError(
                    f"{profile_path} does not match the {taxon} CORE profiles "
                    "in the configured JASPAR2026 database"
                )

    def _database_sha256(self):
        digest = hashlib.sha256()
        with self.sqlite_db.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()

    def _load_inference_data(self):
        if self._inference_data is None:
            self._inference_data = (
                load_cisbp_models(self.data_dir), load_jaspar_data(self.data_dir),
            )
        return self._inference_data

    def check_species_validity(self, species_tax_id: int) -> bool:
        """Check whether local motif annotations include a species.

        Parameters
        ----------
        species_tax_id : int
            Positive NCBI taxonomy identifier to match exactly. Blank and
            null database annotations are treated as missing.

        Returns
        -------
        bool
            Whether the identifier appears in the local species annotations.
        """
        if isinstance(species_tax_id, bool) or not isinstance(species_tax_id, int) or species_tax_id <= 0:
            raise ValueError("species_tax_id must be a positive integer")
        self._validate_configuration()
        if self._species_ids is None:
            try:
                with closing(self._readonly_connection(self.sqlite_db)) as connection:
                    values = connection.execute("SELECT DISTINCT TAX_ID FROM MATRIX_SPECIES")
                    self._species_ids = {
                        int(tax_id)
                        for (tax_id,) in values
                        if tax_id is not None and str(tax_id).strip()
                    }
            except (sqlite3.Error, TypeError, ValueError) as error:
                raise RuntimeError(f"cannot read local JASPAR species annotations: {error}") from error
        return species_tax_id in self._species_ids

    def search(self, request: MotifSearchRequest) -> list[JasparRecord]:
        """Infer local motifs from the supporting ortholog protein window.

        Parameters
        ----------
        request : MotifSearchRequest
            Query and ortholog sequences, supporting domains, and ortholog
            species identifier. Inference uses the ortholog sequence or an
            at-most-2,000-residue window retaining the complete domain.

        Returns
        -------
        list of JasparRecord
            Exact-version predictions passing the inference E-value and
            species filters, limited to ``n_hits``. Empty if none survive.
        """
        if not isinstance(request, MotifSearchRequest):
            raise TypeError("request must be a MotifSearchRequest")
        if (isinstance(request.ortholog_species_tax_id, bool)
                or not isinstance(request.ortholog_species_tax_id, int)
                or request.ortholog_species_tax_id <= 0):
            raise ValueError("ortholog_species_tax_id must be a positive integer")

        sequence = self._window_ortholog_sequence(
            request.ortholog_sequence, request.ortholog_domain,
        )
        self._validate_configuration()
        cisbp_models, jaspar_data = self._load_inference_data()
        temporary_root = Path(tempfile.gettempdir()).resolve()
        repository_root = Path(__file__).resolve().parents[2]
        if temporary_root == repository_root or repository_root in temporary_root.parents:
            for candidate in (Path("/tmp"), Path("/var/tmp")):
                if candidate.is_dir() and candidate.resolve() != repository_root \
                        and repository_root not in candidate.resolve().parents:
                    temporary_root = candidate.resolve()
                    break
            else:
                raise RuntimeError("no temporary directory is available outside the Lverage repository")

        with tempfile.TemporaryDirectory(prefix="lverage-jaspar-", dir=temporary_root) as work_dir:
            predictions = infer_sequence(
                sequence,
                data_dir=self.data_dir,
                cisbp_models=cisbp_models,
                jaspar_data=jaspar_data,
                work_dir=work_dir,
                blastp_path=self.blastp_path,
                hmmscan_path=self.hmmscan_path,
                hmmalign_path=self.hmmalign_path,
                query_id="lverage_query",
            )

        accepted = []
        for prediction in predictions:
            if prediction["evalue"] > self.escore_threshold:
                continue
            record, species_ids = self._resolve_matrix(prediction)
            if request.ortholog_species_tax_id not in species_ids:
                continue
            accepted.append(record)
            if len(accepted) >= self.n_hits:
                break
        return accepted

    def _resolve_matrix(self, inference):
        matrix_id = inference["matrix_id"]
        match = MATRIX_ID_PATTERN.fullmatch(matrix_id)
        if match is None:
            raise ValueError(f"JASPAR matrix ID must be exact and versioned: {matrix_id!r}")
        try:
            with closing(self._readonly_connection(self.sqlite_db)) as connection:
                matrix = connection.execute(
                    "SELECT ID, NAME, COLLECTION FROM MATRIX "
                    "WHERE BASE_ID = ? AND VERSION = ?",
                    (match.group(1), int(match.group(2))),
                ).fetchone()
                if matrix is None or str(matrix[2]).upper() != "CORE":
                    raise ValueError(
                        f"exact matrix {matrix_id} is absent from local JASPAR 2026 CORE data"
                    )
                internal_id, motif_name, _ = matrix
                species = connection.execute(
                    "SELECT TAX_ID FROM MATRIX_SPECIES WHERE ID = ?", (internal_id,),
                ).fetchall()
                annotations = connection.execute(
                    "SELECT VAL FROM MATRIX_ANNOTATION WHERE ID = ? AND TAG = 'class'",
                    (internal_id,),
                ).fetchall()
                pfm = {}
                for base in "ACGT":
                    values = connection.execute(
                        "SELECT VAL FROM MATRIX_DATA WHERE ID = ? AND ROW = ? ORDER BY COL",
                        (internal_id, base),
                    ).fetchall()
                    if not values:
                        raise ValueError(f"local JASPAR matrix {matrix_id} lacks its {base} row")
                    pfm[base] = [float(value[0]) for value in values]
        except sqlite3.Error as error:
            raise RuntimeError(
                f"cannot resolve exact JASPAR matrix {matrix_id} from local data: {error}"
            ) from error

        lengths = {len(values) for values in pfm.values()}
        if len(lengths) != 1 or not lengths.pop():
            raise ValueError(f"local JASPAR matrix {matrix_id} has inconsistent PFM rows")
        if not annotations or not annotations[0][0]:
            raise ValueError(f"local JASPAR matrix {matrix_id} has no TF class annotation")
        try:
            species_ids = {
                int(tax_id)
                for (tax_id,) in species
                if tax_id is not None and str(tax_id).strip()
            }
        except (TypeError, ValueError) as error:
            raise ValueError(f"local JASPAR matrix {matrix_id} has invalid species annotations") from error

        return JasparRecord(
            matrix_id=matrix_id,
            motif_name=motif_name,
            pfm=pfm,
            motif_url=f"https://jaspar.elixir.no/matrix/{matrix_id}/",
            motif_class=annotations[0][0],
            inference_evalue=inference["evalue"],
            inference_tf_name=inference["tf_name"],
            inference_query_region=inference["query_region"],
            inference_tf_region=inference["tf_region"],
            inference_dbd_identity=inference["dbd_identity"],
        ), species_ids

    @staticmethod
    def _window_ortholog_sequence(sequence, ortholog_domain):
        """Return the same at-most-2,000-aa ortholog window used previously."""
        try:
            start = ortholog_domain.start
            end = ortholog_domain.end
        except AttributeError as error:
            raise ValueError("ortholog domains must provide start and end bounds") from error
        if not isinstance(sequence, str) or not sequence:
            raise ValueError("ortholog sequences must be nonempty strings")
        if isinstance(start, bool) or isinstance(end, bool):
            raise ValueError("ortholog domain bounds must be integers")
        if not isinstance(start, int) or not isinstance(end, int):
            raise ValueError("ortholog domain bounds must be integers")
        if start < 0 or end > len(sequence) or start >= end:
            raise ValueError("ortholog domain bounds must define a nonempty slice within the sequence")
        window_length = 2000
        if end - start > window_length:
            raise ValueError("ortholog domains must fit inside the JASPAR sequence window")
        if len(sequence) <= window_length:
            return sequence
        domain_center = (start + end) // 2
        window_start = domain_center - window_length // 2
        window_start = max(0, min(window_start, len(sequence) - window_length))
        return sequence[window_start:window_start + window_length]
