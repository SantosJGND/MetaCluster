# You need to install Biopython: pip install biopython
import gzip
import json
import logging
import os
import re
import subprocess
import time
from collections.abc import Generator
from dataclasses import dataclass
from functools import wraps

import dotenv
import numpy as np
import pandas as pd
from Bio import Entrez

dotenv.load_dotenv()

Entrez.email = os.getenv("NCBI_EMAIL", None)
if Entrez.email is None:
    raise ValueError("NCBI_EMAIL environment variable not set. Please set it to your email address.")

NCBI_TAXONOMY_LEVELS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
NCBI_TAXONOMY_LEVELS_EXTENDED = [
    "acellular root",
    "realm",
    "domain",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]


def is_valid_sequence_file(path: str) -> bool:
    """Return True when ``path`` is a non-empty FASTA file.

    Gzipped files are fully decompressed so truncated files and CRC failures are
    rejected before they can enter the reference cache.
    """
    try:
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            return False
        with open(path, "rb") as handle:
            is_gzip = handle.read(2) == b"\x1f\x8b"
        if is_gzip:
            with gzip.open(path, "rb") as handle:
                if not handle.readline().startswith(b">"):
                    return False
                while handle.read(1024 * 1024):
                    pass
            return True
        with open(path, "rb") as handle:
            return handle.readline().startswith(b">")
    except (OSError, EOFError, gzip.BadGzipFile):
        return False


def _temporary_path(output_path: str) -> str:
    """Return a process-unique temporary path next to ``output_path``."""
    return f"{output_path}.tmp.{os.getpid()}.{time.time_ns()}"


def retry_with_backoff(max_retries=3, initial_delay=1, backoff_factor=2):
    """
    Decorator for retrying functions with exponential backoff.
    Handles transient errors like rate limiting (429) and temporary network issues.
    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            delay = initial_delay
            last_exception = None
            for attempt in range(max_retries):
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    last_exception = e
                    error_msg = str(e).lower()
                    should_retry = any(
                        x in error_msg
                        for x in [
                            "429",
                            "rate limit",
                            "temporary failure",
                            "connection",
                            "timeout",
                            "service unavailable",
                            "500",
                            "502",
                            "503",
                            "504",
                        ]
                    )
                    if not should_retry and attempt > 0:
                        should_retry = attempt < 2

                    if should_retry and attempt < max_retries - 1:
                        logging.warning(
                            f"Attempt {attempt + 1}/{max_retries} failed for {func.__name__}: {e}. Retrying in {delay}s..."
                        )
                        time.sleep(delay)
                        delay *= backoff_factor
                    elif attempt == max_retries - 1:
                        logging.error(f"All {max_retries} attempts failed for {func.__name__}: {e}")
            raise last_exception

        return wrapper

    return decorator


def get_lineage(taxid: str) -> str | None:
    passport = Passport(taxid=taxid, accession=None)
    ncbi_tools = NCBITools()
    lineage = ncbi_tools.retrieve_passport_taxonomy(passport)
    return lineage


class NCBITaxonomistWrapper:
    TAX_SPECIES: str = "species"
    TAX_GENUS: str = "genus"
    TAX_FAMILY: str = "family"
    TAX_ORDER: str = "order"
    TAX_CLASS: str = "class"
    TAX_PHYLUM: str = "phylum"
    TAX_SUPERKINGDOM: str = "superkingdom"

    """Wrapper for ncbi-taxonomist command line tool."""

    def __init__(self, db: str | None = None, temp_dir: str | None = None):
        self.db_path = db
        self.temp_dir = temp_dir
        self.lineages: dict[int, dict[str, dict[str, str]]] = {}
        self.max_fetch = 30
        self.logger = logging.getLogger(__name__)
        logging.basicConfig(level=logging.INFO)
        self.logger.info("NCBI Taxonomist Wrapper initialized.")

    @staticmethod
    def _get_ncbi_bin() -> str:
        """Return path to ncbi-taxonomist binary, preferring venv over pyenv shims."""
        import sys
        from pathlib import Path
        venv_bin = Path(sys.executable).parent / "ncbi-taxonomist"
        if venv_bin.exists():
            return str(venv_bin)
        return "ncbi-taxonomist"

    def retrieve_lineages_cmd_local(self, taxids: list[int]) -> str:
        taxids_str = ",".join(str(t) for t in taxids)
        ncbi_bin = self._get_ncbi_bin()
        cmd = f"{ncbi_bin} resolve -t {taxids_str} -db {self.db_path} "
        self.logger.info("Retrieving taxids from local database")
        return cmd

    def retrieve_lineages_cmd_import(self, taxids: list[int]) -> str:
        taxids_str = ",".join(str(t) for t in taxids)
        ncbi_bin = self._get_ncbi_bin()
        cmd = f"{ncbi_bin} resolve -t {taxids_str} --remote | {ncbi_bin} import --database {self.db_path} "
        return cmd

    def split_taxids(self, taxids: list[int], chunk_size=50) -> Generator[list[int], None, None]:
        for i in range(0, len(taxids), chunk_size):
            yield taxids[i : i + chunk_size]

    def resolve_lineages(self, taxids: list[int]) -> None:
        if len(taxids) == 0:
            return
        taxids = [int(float(t)) for t in taxids]
        cmd = self.retrieve_lineages_cmd_local(taxids)
        stream = os.popen(cmd)
        output = stream.read()
        lineage_dict = self.parse_lineages_output(output)

        missing = set(taxids) - set(lineage_dict.keys())
        self.logger.info(
            f"Retrieved lineages for {len(lineage_dict)} taxids from local database. Missing {len(missing)} taxids."
        )
        if missing:
            for taxid_chunk in self.split_taxids(list(missing), chunk_size=self.max_fetch):
                cmd = self.retrieve_lineages_cmd_import(list(taxid_chunk))
                stream = os.popen(cmd)
                output = stream.read()
                lineage_dict_update = self.parse_lineages_output(output)
                lineage_dict.update(lineage_dict_update)
                time.sleep(3)

        self.lineages.update(lineage_dict)

        still_missing = set(taxids) - set(self.lineages.keys())
        if still_missing:
            ncbi_tools = NCBITools()
            self.logger.info(f"Still missing taxids: {len(still_missing)}")
            for taxid_chunk in self.split_taxids(list(still_missing), chunk_size=self.max_fetch):
                lineages = ncbi_tools.retrieve_taxonomies_batch(list(taxid_chunk))
                for taxid, lineage in lineages.items():
                    if lineage is not None:
                        parser = NCBIlineageParser(lineage)
                        self.lineages[taxid] = {
                            level: {"name": parser.get_level(level), "taxid": None, "level": i}
                            for i, level in enumerate(NCBI_TAXONOMY_LEVELS_EXTENDED)
                            if parser.get_level(level) is not None
                        }
                time.sleep(3)

        self.logger.info(f"Updated lineage dictionary with {len(self.lineages)} entries.")
        self.logger.info(f"Still missing taxids after NCBI fetch: {len(set(taxids) - set(self.lineages.keys()))}")

    def update_lineages(self, taxids: list[int]) -> None:
        missing = set(taxids) - set(self.lineages.keys())
        if missing:
            self.logger.info(f"Missing taxids: {len(missing)}. Trying to import them from NCBI...")
            for taxid_chunk in self.split_taxids(list(missing), chunk_size=self.max_fetch):
                cmd = self.retrieve_lineages_cmd_import(list(taxid_chunk))
                stream = os.popen(cmd)
                output = stream.read()
                lineage_dict_update = self.parse_lineages_output(output)
                self.lineages.update(lineage_dict_update)
                time.sleep(3)

        still_missing = set(taxids) - set(self.lineages.keys())
        if still_missing:
            self.logger.info(f"Still missing taxids: {len(still_missing)}.")
            for taxid in still_missing:
                lineage = get_lineage(str(taxid))
                if lineage is not None:
                    parser = NCBIlineageParser(lineage)
                    self.lineages[taxid] = {
                        level: {"name": parser.get_level(level), "taxid": None, "level": i}
                        for i, level in enumerate(NCBI_TAXONOMY_LEVELS_EXTENDED)
                        if parser.get_level(level) is not None
                    }

        self.logger.info(f"Updated lineage dictionary with {len(self.lineages)} entries.")
        self.logger.info(f"Still missing taxids after NCBI fetch: {len(set(taxids) - set(self.lineages.keys()))}")

    def parse_lineages_output(self, output: str) -> dict:
        lineage_dict = {}
        for line in output.splitlines():
            try:
                data = json.loads(line)
                if "query" in data and "lineage" in data:
                    len_lineage = len(data["lineage"])
                    lineage_dict[int(data["query"])] = {
                        x["rank"]: {"name": x["name"], "taxid": int(x["taxid"]), "level": len_lineage - i}
                        for i, x in enumerate(data["lineage"])
                    }
            except json.JSONDecodeError:
                continue
        return lineage_dict

    def get_name(self, taxid: int) -> str | None:
        if taxid not in self.lineages:
            return None
        lineage = self.lineages.get(taxid, {})
        lineage = sorted(lineage.items(), key=lambda x: x[1]["level"], reverse=True)
        if lineage:
            return lineage[0][1]["name"]
        return None

    def get_level(self, taxid: int, level: str) -> str | None:
        if taxid not in self.lineages:
            return None
        return self.lineages.get(taxid, {}).get(level, {}).get("name", None)

    def get_lineage(self, taxid: int) -> list:
        if taxid not in self.lineages:
            return []
        lineage = self.lineages.get(taxid, {})
        lineage = sorted(lineage.items(), key=lambda x: x[1]["level"])
        return [(rank, info["name"]) for rank, info in lineage]

    def compare_lineages_relative(self, taxid1: int, taxid2: int) -> tuple[float, str | None]:
        if taxid1 not in self.lineages:
            return 0.0, None
        if taxid2 not in self.lineages:
            return 0.0, None

        min_levels = min(
            min([x["level"] for x in self.lineages[taxid1].values()]),
            min([x["level"] for x in self.lineages[taxid2].values()]),
        )

        if min_levels == 0:
            lineage1 = self.get_lineage(taxid1)
            lineage2 = self.get_lineage(taxid2)
            tax, level = compare_lineages("; ".join([x[1] for x in lineage1]), "; ".join([x[1] for x in lineage2]))
            return tax, level

        lineage1 = sorted(self.lineages[taxid1].items(), key=lambda x: x[1]["level"])

        for tax, level in lineage1:
            if self.get_level(taxid2, tax) != level["name"]:
                return (level.get("level", 0) - 1) / len(lineage1), tax

        return 1.0, lineage1[-1][0] if lineage1 else None

    @staticmethod
    def level_is_atleast(level: str, threshold_level="species") -> bool:
        if threshold_level not in NCBI_TAXONOMY_LEVELS_EXTENDED:
            raise ValueError(f"Invalid threshold level: {threshold_level}")
        if level not in NCBI_TAXONOMY_LEVELS_EXTENDED:
            return True
        return NCBI_TAXONOMY_LEVELS_EXTENDED.index(level) >= NCBI_TAXONOMY_LEVELS_EXTENDED.index(threshold_level)

    @staticmethod
    def level_is_below(level: str, threshold_level="species") -> bool:
        if threshold_level not in NCBI_TAXONOMY_LEVELS_EXTENDED:
            raise ValueError(f"Invalid threshold level: {threshold_level}")
        if level not in NCBI_TAXONOMY_LEVELS_EXTENDED:
            return False
        return NCBI_TAXONOMY_LEVELS_EXTENDED.index(level) < NCBI_TAXONOMY_LEVELS_EXTENDED.index(threshold_level)


class NCBIlineageParser:
    def __init__(self, lineage: str | None):
        self.lineage = lineage
        self.levels = {}
        if lineage is not None:
            parts = [part.strip() for part in lineage.split(";")]
            for i, part in enumerate(parts):
                if i < len(NCBI_TAXONOMY_LEVELS_EXTENDED):
                    self.levels[NCBI_TAXONOMY_LEVELS_EXTENDED[i]] = part
                else:
                    self.levels[f"level_{i + 1}"] = part

    def get_level(self, level: str) -> str | None:
        return self.levels.get(level, None)

    def get_order(self) -> str | None:
        return self.get_level("order")

    def get_family(self) -> str | None:
        return self.get_level("family")

    def get_genus(self) -> str | None:
        return self.get_level("genus")

    def get_species(self) -> str | None:
        return self.get_level("species")


def compare_lineages(lineage1: str | None, lineage2: str | None) -> tuple[float, str | None]:
    if lineage1 is None or lineage2 is None:
        return 0.0, None

    if pd.isna(lineage1) or pd.isna(lineage2):
        return 0.0, None

    levels1 = lineage1.split("; ")
    levels2 = lineage2.split("; ")
    min_length = min(len(levels1), len(levels2))
    score = 0
    level = None
    for i in range(min_length):
        if levels1[i] == levels2[i]:
            score += 1
            if i < len(NCBI_TAXONOMY_LEVELS):
                level = NCBI_TAXONOMY_LEVELS[i]
            else:
                level = f"level_{i + 1}"
        else:
            break
    if min_length == 0:
        return 0.0, level
    return score / min_length, level


def normalize_accession(value: str | None) -> tuple[str | None, str]:
    """Return the bare accession and its Entrez database for an EsViritu identifier.

    EsViritu ``Assembly`` values may be prefixed with ``acc:`` / ``set:`` for
    custom database entries; RefSeq assembly accessions (``GCF_`` / ``GCA_``)
    live in the assembly database, everything else in nucleotide.
    """
    if value is None:
        return None, "nucleotide"
    accession = str(value).strip().strip('"').strip("'")
    for prefix in ("acc:", "set:"):
        if accession.lower().startswith(prefix):
            accession = accession[len(prefix) :].strip()
            break
    if not accession:
        return None, "nucleotide"
    db = "assembly" if re.match(r"^GC[AF]_\d+(?:\.\d+)?$", accession, re.IGNORECASE) else "nucleotide"
    return accession, db


@retry_with_backoff(max_retries=3, initial_delay=1)
def retrieve_taxid_from_accession(accession: str, db: str | None = None) -> int | None:
    """Resolve an NCBI nucleotide or assembly accession to its taxonomy id."""
    accession, inferred_db = normalize_accession(accession)
    if accession is None:
        return None
    db = "assembly" if inferred_db == "assembly" else (db or "nucleotide")
    if db not in {"nucleotide", "assembly"}:
        raise ValueError(f"Unsupported Entrez database: {db}")

    term = f"{accession}[Accession]" if db == "nucleotide" else f"{accession}[Assembly Accession]"
    handle = Entrez.esearch(db=db, term=term, retmax=1)
    record = Entrez.read(handle)
    handle.close()
    ids = record.get("IdList") or []
    if not ids:
        return None
    uid = ids[0]

    try:
        handle = Entrez.elink(dbfrom=db, db="taxonomy", id=uid)
        links = Entrez.read(handle)
        handle.close()
        for linkset in links:
            for linkdb in linkset.get("LinkSetDb", []):
                linked = linkdb.get("Link") or []
                if linked:
                    return int(linked[0]["Id"])
    except Exception:
        pass

    try:
        handle = Entrez.esummary(db=db, id=uid, report="full")
        summary = Entrez.read(handle)
        handle.close()
        if db == "assembly":
            docsum = summary["DocumentSummarySet"]["DocumentSummary"][0]
            taxid = docsum.get("Taxid") or docsum.get("SpeciesTaxid")
        else:
            taxid = summary[0].get("TaxId")
        if taxid:
            return int(taxid)
    except Exception:
        pass

    return None


ACCESSION_SPECIES_CACHE_PATH = os.path.join(
    os.path.expanduser("~"), ".cache", "metagenomics-evaluation-pipeline", "accession_species_taxid.json"
)


def _load_accession_cache(cache_path: str) -> dict[str, int | None]:
    if not os.path.exists(cache_path):
        return {}
    with open(cache_path) as handle:
        return {acc: (int(tax) if tax is not None else None) for acc, tax in json.load(handle).items()}


def _save_accession_cache(cache: dict[str, int | None], cache_path: str) -> None:
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    tmp_path = _temporary_path(cache_path)
    with open(tmp_path, "w") as handle:
        json.dump(cache, handle, indent=1, sort_keys=True)
    os.replace(tmp_path, cache_path)


@retry_with_backoff(max_retries=3, initial_delay=1)
def _fetch_species_taxid(accession: str) -> int | None:
    """Return the species-rank taxid of a nucleotide accession, or None if it cannot be determined."""
    handle = Entrez.esearch(db="nucleotide", term=accession, retmax=1)
    record = Entrez.read(handle)
    handle.close()
    if not record["IdList"]:
        return None

    handle = Entrez.esummary(db="nucleotide", id=record["IdList"][0])
    summary = Entrez.read(handle)
    handle.close()
    tax_id = summary[0].get("TaxId")
    if not tax_id:
        return None

    handle = Entrez.efetch(db="taxonomy", id=str(tax_id), retmode="xml")
    tax_records = Entrez.read(handle)
    handle.close()
    if not tax_records:
        return None
    node = tax_records[0]
    if node.get("Rank") == "species":
        return int(node["TaxId"])
    for ancestor in node.get("LineageEx", []):
        if ancestor.get("Rank") == "species":
            return int(ancestor["TaxId"])
    return None


def resolve_accession_species_taxids(
    accessions: list[str], cache_path: str = ACCESSION_SPECIES_CACHE_PATH, request_delay: float = 0.4
) -> dict[str, int | None]:
    """
    Map nucleotide/assembly accessions to species-rank taxids via Entrez, with a JSON cache.

    Successful lookups and definitive "no species" answers are cached. Network errors are not
    cached, so they are retried on the next run. Accessions that cannot be resolved map to None.
    """
    cache = _load_accession_cache(cache_path)
    result: dict[str, int | None] = {}
    new_entries = False
    for accession in accessions:
        if accession in cache:
            result[accession] = cache[accession]
            continue
        try:
            species_taxid = _fetch_species_taxid(accession)
        except Exception as e:
            logging.warning(f"Species taxid lookup failed for {accession}: {e}")
            result[accession] = None
            continue
        cache[accession] = species_taxid
        result[accession] = species_taxid
        new_entries = True
        time.sleep(request_delay)

    if new_entries:
        _save_accession_cache(cache, cache_path)
    return result


TAXID_LINEAGE_CACHE_PATH = os.path.join(os.path.dirname(ACCESSION_SPECIES_CACHE_PATH), "taxid_lineage.json")
ACCESSION_TAXID_CACHE_PATH = os.path.join(os.path.dirname(ACCESSION_SPECIES_CACHE_PATH), "accession_taxid.json")

# Ranks a retrieved taxid may have to be rewritten to a truth. Clades (e.g. 262897) are excluded.
SPECIES_OR_BELOW_RANKS = {"species", "subspecies", "strain", "no rank"}


@retry_with_backoff(max_retries=3, initial_delay=1)
def _fetch_accession_taxid(accession: str) -> int | None:
    """Return the taxid Entrez reports for a nucleotide accession (esummary TaxId), or None."""
    handle = Entrez.esearch(db="nucleotide", term=accession, retmax=1)
    record = Entrez.read(handle)
    handle.close()
    if not record["IdList"]:
        return None
    handle = Entrez.esummary(db="nucleotide", id=record["IdList"][0])
    summary = Entrez.read(handle)
    handle.close()
    tax_id = summary[0].get("TaxId")
    return int(tax_id) if tax_id else None


def resolve_accession_taxids(
    accessions: list[str], cache_path: str = ACCESSION_TAXID_CACHE_PATH, request_delay: float = 0.4
) -> dict[str, int | None]:
    """
    Map accessions to the taxid Entrez reports for them, with a JSON cache.

    Successful lookups and "no record" answers are cached. Network errors are not cached.
    """
    cache = _read_json_cache(cache_path)
    result: dict[str, int | None] = {}
    new_entries = False
    for accession in accessions:
        if accession in cache:
            result[accession] = cache[accession]
            continue
        try:
            taxid = _fetch_accession_taxid(accession)
        except Exception as e:
            logging.warning(f"Taxid lookup failed for {accession}: {e}")
            result[accession] = None
            continue
        cache[accession] = taxid
        result[accession] = taxid
        new_entries = True
        time.sleep(request_delay)

    if new_entries:
        _write_json_cache(cache, cache_path)
    return result


# Ranks above species. Only truth taxids at these ranks are rewritten to the retrieved taxid.
ABOVE_SPECIES_RANKS = {
    "genus",
    "subgenus",
    "tribe",
    "subfamily",
    "family",
    "superfamily",
    "suborder",
    "order",
    "subclass",
    "class",
    "subphylum",
    "phylum",
    "kingdom",
    "superkingdom",
    "domain",
    "realm",
}


def _read_json_cache(cache_path: str) -> dict:
    if not os.path.exists(cache_path):
        return {}
    with open(cache_path) as handle:
        return json.load(handle)


def _write_json_cache(cache: dict, cache_path: str) -> None:
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    tmp_path = _temporary_path(cache_path)
    with open(tmp_path, "w") as handle:
        json.dump(cache, handle, indent=1, sort_keys=True)
    os.replace(tmp_path, cache_path)


@retry_with_backoff(max_retries=3, initial_delay=1)
def _fetch_taxid_lineage(taxid: int) -> tuple[str | None, set[int], int | None]:
    """
    Return ``(rank, lineage_taxids, species_ancestor)`` for a taxid.

    The lineage includes the taxid itself. ``species_ancestor`` is the nearest species-rank node
    (the taxid itself when it is species rank), or None when there is none.
    """
    handle = Entrez.efetch(db="taxonomy", id=str(taxid), retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    if not records:
        return None, set(), None
    node = records[0]
    ancestors = node.get("LineageEx", [])
    lineage = {int(x["TaxId"]) for x in ancestors}
    lineage.add(int(node["TaxId"]))

    species_ancestor = None
    if node.get("Rank") == "species":
        species_ancestor = int(node["TaxId"])
    else:
        for x in reversed(ancestors):
            if x.get("Rank") == "species":
                species_ancestor = int(x["TaxId"])
                break
    return node.get("Rank"), lineage, species_ancestor


def resolve_taxid_lineages(
    taxids: list[int], cache_path: str = TAXID_LINEAGE_CACHE_PATH, request_delay: float = 0.4
) -> dict[int, tuple[str | None, set[int], int | None]]:
    """
    Map taxids to ``(rank, lineage_taxids, species_ancestor)`` via Entrez, with a JSON cache.

    Cache entries written before ``species_ancestor`` was stored are refetched. Failed lookups are
    not cached and are missing from the result.
    """
    raw = _read_json_cache(cache_path)
    result: dict[int, tuple[str | None, set[int], int | None]] = {}
    new_entries = False
    for taxid in taxids:
        key = str(taxid)
        if key in raw and "species_ancestor" in raw[key]:
            entry = raw[key]
            result[taxid] = (entry["rank"], set(entry["lineage"]), entry["species_ancestor"])
            continue
        try:
            rank, lineage, species_ancestor = _fetch_taxid_lineage(taxid)
        except Exception as e:
            logging.warning(f"Taxonomy lookup failed for taxid {taxid}: {e}")
            continue
        if rank is None:
            continue
        raw[key] = {"rank": rank, "lineage": sorted(lineage), "species_ancestor": species_ancestor}
        result[taxid] = (rank, lineage, species_ancestor)
        new_entries = True
        time.sleep(request_delay)

    if new_entries:
        _write_json_cache(raw, cache_path)
    return result


def _as_int(value) -> int | None:
    if value is None or pd.isna(value):
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def resolve_species_truth_taxids(
    input_df: pd.DataFrame,
    accession_cache_path: str = ACCESSION_TAXID_CACHE_PATH,
    lineage_cache_path: str = TAXID_LINEAGE_CACHE_PATH,
) -> pd.DataFrame:
    """
    Rewrite truth taxids of a simulation input table to the taxid Entrez reports for ``assembly_accession``.

    A row is rewritten only when all of these hold:
      - the input taxid is above species rank (see ``ABOVE_SPECIES_RANKS``);
      - the retrieved taxid is species rank or below (see ``SPECIES_OR_BELOW_RANKS``) and has a
        species ancestor, so clades such as 262897 are not rewritten;
      - the input taxid lies in the retrieved taxid's lineage.

    Every other row keeps its input taxid. The original value is kept in ``taxid_input_rank``.
    Outcome counts are logged per reason.
    """
    if "assembly_accession" not in input_df.columns:
        logging.warning("Input table has no assembly_accession column; truth taxids kept at input rank")
        return input_df
    if input_df.empty:
        return input_df

    accessions = input_df["assembly_accession"].dropna().astype(str).str.strip().unique().tolist()
    retrieved = resolve_accession_taxids(accessions, cache_path=accession_cache_path)
    retrieved_ids = [
        _as_int(retrieved.get(str(acc).strip())) if pd.notna(acc) else None for acc in input_df["assembly_accession"]
    ]
    input_ids = [_as_int(t) for t in input_df["taxid"]]

    needed = {t for t in input_ids if t is not None} | {r for r in retrieved_ids if r is not None}
    lineages = resolve_taxid_lineages(sorted(needed), cache_path=lineage_cache_path)

    reasons: list[str] = []
    new_ids: list[int | None] = []
    for t, r in zip(input_ids, retrieved_ids):
        if t is None:
            reason, new = "input_taxid_missing", None
        elif r is None:
            reason, new = "accession_unresolved", t
        elif t == r:
            reason, new = "already_retrieved", t
        elif t not in lineages:
            reason, new = "input_rank_unresolved", t
        elif lineages[t][0] not in ABOVE_SPECIES_RANKS:
            reason, new = "input_at_or_below_species", t
        elif r not in lineages:
            reason, new = "retrieved_lineage_unresolved", t
        elif lineages[r][0] not in SPECIES_OR_BELOW_RANKS or lineages[r][2] is None:
            reason, new = "retrieved_above_species", t
        elif t in lineages[r][1]:
            reason, new = "rewritten", r
        else:
            reason, new = "not_descendant", t
        reasons.append(reason)
        new_ids.append(new)

    out = input_df.copy()
    out["taxid_input_rank"] = input_df["taxid"]
    rewritten = pd.to_numeric(pd.Series(new_ids, index=input_df.index, dtype=object), errors="coerce")
    if rewritten.notna().all():
        rewritten = rewritten.astype(int)
    out["taxid"] = rewritten

    counts = pd.Series(reasons).value_counts().to_dict()
    logging.info(f"Truth taxid rewrite outcomes: {counts}")
    return out


@dataclass
class Passport:
    taxid: str | None
    accession: str | None = None
    lineage: str | None = None
    description: str | None = None

    def __init__(
        self,
        taxid: str | None,
        accession: str | None = None,
        lineage: str | None = None,
        description: str | None = None,
    ):
        if taxid and "." in str(taxid):
            taxid = str(taxid).split(".")[0]
        self.taxid = taxid
        self.accession = accession
        self.lineage = lineage
        self.description = description

    def __str__(self):
        return f"TaxID: {self.taxid}, Accession: {self.accession}"

    @property
    def prefix(self):
        if self.accession:
            return f"{self.taxid}_{self.accession}"
        else:
            return f"{self.taxid}"


@dataclass
class LocalAssembly(Passport):
    file_path: str | None = None

    def __str__(self):
        return f"TaxID: {self.taxid}, Accession: {self.accession}, File Path: {self.file_path}"


@dataclass
class ReferenceData(Passport):
    nucleotide_id: str | None = None
    assembly_id: str | None = None

    def __str__(self):
        return f"TaxID: {self.taxid}, Accession: {self.accession}, Description: {self.description}, Nucleotide ID: {self.nucleotide_id}, Assembly ID: {self.assembly_id}"


@retry_with_backoff(max_retries=3, initial_delay=1)
def retrieve_passport_taxonomy(taxid: str) -> str | None:
    try:
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        if not records:
            raise ValueError(f"No taxonomy found for taxid {taxid}")
        lineage = records[0]["Lineage"]
        return lineage
    except Exception:
        return None


@retry_with_backoff(max_retries=3, initial_delay=1)
def retrieve_reference_sequence_id(
    accID: str, include_term=None, exclude_term=None
) -> tuple[str | None, str | None, str | None]:
    try:
        term = accID
        if exclude_term is not None:
            term += f" NOT {exclude_term}"
        if include_term is not None:
            term += f" AND {include_term}"

        handle = Entrez.esearch(db="nucleotide", term=term, retmax=20)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            raise ValueError(f"No sequence found for accession {accID}")

        sequence_id = record["IdList"][0]

        handle = Entrez.esummary(db="nucleotide", id=sequence_id)
        summary = Entrez.read(handle)
        handle.close()
        if summary is None:
            raise ValueError(f"No summary found for sequence ID {sequence_id}")
        accession = summary[0]["AccessionVersion"]
        if accession != accID:
            print(f"Warning: Retrieved accession {accession} does not match requested {accID}")
        description = summary[0].get("Title", "No description available")
        return accession, description, sequence_id

    except ValueError as e:
        print(e)
        return None, None, None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None, None, None


@retry_with_backoff(max_retries=3, initial_delay=1)
def get_reference_sequence_url(
    taxid, include_term=None, exclude_term=None
) -> tuple[str | None, str | None, str | None]:
    try:
        term = f"txid{taxid}[Organism:exp] AND refseq"
        if include_term is not None:
            term += f" AND {include_term}"
        if exclude_term is not None:
            term += f" NOT {exclude_term}"

        handle = Entrez.esearch(db="nucleotide", term=term, retmax=1)
        record = Entrez.read(handle)
        handle.close()
        if not record["IdList"]:
            raise ValueError(f"No reference sequences found for taxid {taxid}")

        sequence_id = record["IdList"][0]

        handle = Entrez.esummary(db="nucleotide", id=sequence_id)
        summary = Entrez.read(handle)
        handle.close()

        if summary is None:
            raise ValueError(f"No summary found for sequence ID {sequence_id}")

        docsum = summary[0]
        accession = docsum["AccessionVersion"]
        description = docsum.get("Title", "No description available")

        return accession, description, sequence_id

    except ValueError as e:
        print(e)
        return None, None, None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None, None, None


@retry_with_backoff(max_retries=3, initial_delay=1)
def get_representative_assembly(
    taxid, include_term=None, exclude_term=None
) -> tuple[str | None, str | None, str | None]:
    try:
        term = f"txid{taxid}[Organism:exp]"
        if include_term is not None:
            term += f" AND {include_term}"
        if exclude_term is not None:
            term += f" NOT {exclude_term}"

        handle = Entrez.esearch(db="assembly", term=term, retmax=5)
        record = Entrez.read(handle)
        handle.close()
        if not record["IdList"]:
            raise ValueError(f"No representative genomes found for taxid {taxid}")

        assembly_id = record["IdList"][0]

        handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
        summary = Entrez.read(handle)
        handle.close()
        docsum = summary["DocumentSummarySet"]["DocumentSummary"][0]
        accession = docsum["AssemblyAccession"]
        description = docsum.get("SpeciesName", "No description available")

        return accession, description, assembly_id

    except ValueError as e:
        print(e)
        return None, None, None
    except Exception as e:
        print(f"An error occurred while fetching assembly for taxid {taxid}: {e}")
        return None, None, None


@retry_with_backoff(max_retries=3, initial_delay=1)
def retrieve_reference_sequence(nucleotide_id, output_path, gzipped=True) -> bool:
    tmp_path = _temporary_path(output_path)
    try:
        handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, rettype="fasta", retmode="text")
        fasta_data = handle.read()
        handle.close()
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        if gzipped:
            with gzip.open(tmp_path, "wt") as f:
                f.write(fasta_data)
        else:
            with open(tmp_path, "w") as f:
                f.write(fasta_data)
        if not is_valid_sequence_file(tmp_path):
            raise ValueError(f"Downloaded sequence is not a valid FASTA file: {output_path}")
        os.replace(tmp_path, output_path)
        return True
    except Exception as e:
        print(f"An error occurred while downloading sequence: {e}")
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        return False


@retry_with_backoff(max_retries=3, initial_delay=1)
def retrieve_assembly_sequence(assembly_id, output_path) -> bool:
    tmp_path = _temporary_path(output_path)
    try:
        handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
        summary = Entrez.read(handle)
        handle.close()
        docsum = summary["DocumentSummarySet"]["DocumentSummary"][0]
        ftp_path = docsum["FtpPath_RefSeq"] or docsum["FtpPath_GenBank"]
        if not ftp_path:
            print(f"No FTP path found for assembly ID {assembly_id}")
            return False
        asm_name = ftp_path.split("/")[-1]
        fasta_url = f"{ftp_path}/{asm_name}_genomic.fna.gz"
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        result = subprocess.run(["wget", "-O", tmp_path, fasta_url], capture_output=True, check=False)
        if result.returncode != 0:
            raise RuntimeError(f"Failed to download file: {result.stderr.decode()}")
        if not is_valid_sequence_file(tmp_path):
            raise ValueError(f"Downloaded assembly is not a valid gzipped FASTA file: {output_path}")
        os.replace(tmp_path, output_path)
        return True
    except Exception as e:
        print(f"An error occurred while downloading assembly: {e}")
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        return False


class NCBITools:
    def __init__(self):
        self.logger = logging.getLogger("NCBITools")
        self.logger.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)
        self.logger.propagate = False

    def retrieve_taxid_from_accession(self, accession: str, db: str | None = None) -> int | None:
        """Resolve an accession to a taxid, logging when it cannot be resolved."""
        taxid = retrieve_taxid_from_accession(accession, db)
        if taxid is None:
            self.logger.warning(f"Could not resolve accession {accession!r} to a taxid")
        return taxid

    @retry_with_backoff(max_retries=3, initial_delay=1)
    def retrieve_passport_taxonomy(self, passport: Passport) -> str | None:
        try:
            handle = Entrez.efetch(db="taxonomy", id=passport.taxid, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            if not records:
                raise ValueError(f"No taxonomy found for taxid {passport.taxid}")
            lineage = records[0]["Lineage"]
            return lineage
        except Exception as e:
            self.logger.error(f"An error occurred while fetching taxonomy for taxid {passport.taxid}: {e}")
            return None

    @retry_with_backoff(max_retries=3, initial_delay=1)
    def retrieve_taxonomies_batch(self, taxids: list[int]) -> dict:
        taxid_list = ",".join(map(str, taxids))
        fetch_cmd = [
            "efetch",
            "-db",
            "taxonomy",
            "-id",
            taxid_list,
            "-format",
            "xml",
            "|",
            "xtract",
            "-pattern",
            "Taxon",
            "-element",
            "TaxId",
            "Lineage",
        ]
        try:
            result = subprocess.run(" ".join(fetch_cmd), shell=True, capture_output=True, text=True)
            if result.returncode != 0:
                self.logger.error(f"Failed to retrieve taxonomies for taxids {taxids}: {result.stderr}")
                return {}
            taxonomies = {}
            for line in result.stdout.splitlines():
                line = line.split()
                taxids = []
                lineage = ""
                for part in line:
                    if part.isdigit():
                        taxids.append(int(part))
                    else:
                        lineage += part + " "

                if taxids:
                    taxonomies.update(dict.fromkeys(taxids, lineage.strip()))

            return taxonomies

        except Exception as e:
            import traceback

            self.logger.error(f"An error occurred while retrieving taxonomies for taxids {taxids}: {e}")
            self.logger.error(traceback.format_exc())
            return {}

    def query_sequence_databases(
        self, passport: Passport, include_term: str | None = None, exclude_term: str | None = None
    ) -> ReferenceData:
        if passport.taxid is None or pd.isna(passport.taxid) or np.isnan(float(passport.taxid)):
            self.logger.error("No taxid provided in passport.")
            return ReferenceData(taxid=None, accession=None, lineage=None, description=None)

        lineage = self.retrieve_passport_taxonomy(passport)
        self.logger.info(f"Lineage for taxid {passport.taxid}: {lineage}")
        passport.lineage = lineage

        if passport.accession is not None:
            accession, description, nucleotide_id = retrieve_reference_sequence_id(
                passport.accession, include_term, exclude_term
            )
            if nucleotide_id is not None:
                return ReferenceData(
                    taxid=passport.taxid,
                    accession=passport.accession,
                    lineage=passport.lineage,
                    description=description,
                    nucleotide_id=nucleotide_id,
                )
            else:
                self.logger.warning(
                    f"No nucleotide ID found for accession {passport.accession}, falling back to taxid search."
                )

        accession, description, nucleotide_id = get_reference_sequence_url(passport.taxid, include_term, exclude_term)

        if accession is not None and nucleotide_id is not None:
            return ReferenceData(
                taxid=passport.taxid,
                accession=accession,
                lineage=passport.lineage,
                description=description,
                nucleotide_id=nucleotide_id,
            )

        accession, description, assembly_id = get_representative_assembly(passport.taxid, include_term, exclude_term)

        return ReferenceData(
            taxid=passport.taxid,
            accession=accession,
            lineage=passport.lineage,
            description=description,
            assembly_id=assembly_id,
        )

    def retrieve_sequence_databases(
        self, reference_data: ReferenceData, output_path: str, gzipped: bool = True
    ) -> bool:
        if reference_data.nucleotide_id is not None:
            success = retrieve_reference_sequence(reference_data.nucleotide_id, output_path, gzipped)
            if not success:
                self.logger.warning(f"Failed to retrieve reference sequence for taxid {reference_data.taxid}")
            return success

        if reference_data.assembly_id is not None:
            success = retrieve_assembly_sequence(reference_data.assembly_id, output_path)
            if not success:
                self.logger.warning(f"Failed to retrieve assembly sequence for taxid {reference_data.taxid}")
            return success

        self.logger.warning(f"No sequence data found for taxid {reference_data.taxid}")
        return False
