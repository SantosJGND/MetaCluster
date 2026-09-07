import logging
import os
import re
import time

import pandas as pd

from metagenomics_utils.dataframe_utils import detect_id_columns, rename_columns_to_standard
from metagenomics_utils.ncbi_tools import LocalAssembly, NCBITools, Passport, ReferenceData


class RateLimiter:
    """
    Simple rate limiter to prevent hitting NCBI API rate limits.
    """

    def __init__(self, delay_between_calls: float = 0.5):
        self.delay = delay_between_calls
        self.last_call = 0

    def wait(self):
        elapsed = time.time() - self.last_call
        if elapsed < self.delay:
            time.sleep(self.delay - elapsed)
        self.last_call = time.time()


# Global rate limiter instance
_rate_limiter = RateLimiter(delay_between_calls=0.5)

_NCBI_ACCESSION_RE = re.compile(r"(GC[AF]_\d+\.\d+|[A-Z]{1,2}_?\d{5,}\.\d+)")


def _accession_from_filename(filename: str, taxid: str | None = None) -> str | None:
    """Best-effort accession extraction from a stored sequence filename."""
    m = _NCBI_ACCESSION_RE.search(filename)
    if m:
        return m.group(0)
    stem = filename
    for suffix in ("_sequence.fasta.gz", ".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    if taxid and stem.startswith(f"{taxid}_"):
        stem = stem[len(taxid) + 1 :]
    return stem or None


class AssemblyStore:
    """
    Class to manage assembly storage and retrieval.
    """

    def __init__(self, store_path: str):
        self.store_path = store_path
        os.makedirs(self.store_path, exist_ok=True)
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)
        handler = logging.StreamHandler()
        handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)
        self.logger.propagate = False

        self.ncbi = NCBITools()
        self.last_failed_taxids: pd.DataFrame = pd.DataFrame(columns=["taxid", "accession", "error"])

    def get_assembly_path(self, taxid: str) -> str:
        return os.path.join(self.store_path, taxid)

    _SEQUENCE_SUFFIXES = (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa")

    def _candidate_files(self, directory: str, accession: str | None) -> list[str]:
        """Sequence files in ``directory``, optionally restricted to those mentioning ``accession``."""
        if not os.path.isdir(directory):
            return []
        files = sorted(f for f in os.listdir(directory) if f.endswith(self._SEQUENCE_SUFFIXES))
        if accession:
            files = [f for f in files if accession in f]
        return [os.path.join(directory, f) for f in files]

    def retrieve_local_assembly(self, passport: Passport) -> LocalAssembly | None:
        """
        Locate an assembly for ``passport`` in the store.

        Lookup order:
        1. exact ``{store}/{taxid}/{taxid}_{accession}_sequence.fasta.gz``;
        2. any sequence file under ``{store}/{taxid}/`` containing the accession
           (or any sequence file there when no accession is known);
        3. any sequence file anywhere in the store containing the accession.
        """
        taxid_subdir = os.path.join(self.store_path, str(passport.taxid))
        exact = os.path.join(taxid_subdir, f"{passport.prefix}_sequence.fasta.gz")
        if os.path.exists(exact):
            return LocalAssembly(taxid=passport.taxid, accession=passport.accession, file_path=exact)

        candidates = self._candidate_files(taxid_subdir, passport.accession)
        if not candidates and passport.accession is None:
            candidates = self._candidate_files(taxid_subdir, None)
        if not candidates and passport.accession:
            for entry in sorted(os.listdir(self.store_path)) if os.path.isdir(self.store_path) else []:
                candidates.extend(self._candidate_files(os.path.join(self.store_path, entry), passport.accession))
                if candidates:
                    break

        if not candidates:
            self.logger.warning(f"No assembly file found for taxid {passport.taxid} and accession {passport.accession}")
            return None

        chosen = candidates[0]
        if len(candidates) > 1:
            self.logger.info(f"Multiple local assemblies for taxid {passport.taxid}; using {chosen}")
        accid = passport.accession or _accession_from_filename(os.path.basename(chosen), str(passport.taxid))
        return LocalAssembly(taxid=passport.taxid, accession=accid, file_path=chosen)

    def retrieve_assembly(
        self,
        passport: Passport,
        reference_data: ReferenceData | None = None,
        include_term: str | None = None,
        exclude_term: str | None = None,
    ) -> LocalAssembly | None:
        """
        Retrieve the assembly for the given taxid, either from local storage or NCBI.
        """
        # First, check if the assembly is available locally
        local_assembly = self.retrieve_local_assembly(passport)

        if local_assembly:
            self.logger.info(f"Using local assembly for taxid {passport.taxid}: {local_assembly.file_path}")
            return local_assembly

        # If not found locally, fetch from NCBI
        self.logger.info(f"Fetching assembly for taxid {passport.taxid} from NCBI...")
        #
        if reference_data is None:
            reference_data = self.ncbi.query_sequence_databases(
                passport, include_term=include_term, exclude_term=exclude_term
            )
        assembly_dir = os.path.join(self.store_path, str(passport.taxid))
        os.makedirs(assembly_dir, exist_ok=True)

        assembly_file_path = os.path.join(assembly_dir, f"{reference_data.prefix}_sequence.fasta.gz")
        success_dl = self.ncbi.retrieve_sequence_databases(reference_data, assembly_file_path, gzipped=True)

        if not success_dl:
            self.logger.error(f"Failed to download assembly for passport {passport.taxid}")
            return None

        return LocalAssembly(taxid=passport.taxid, accession=reference_data.accession, file_path=assembly_file_path)

    def match_taxid_to_assembly(
        self, classification_output_path: str, include_term: str | None = None, exclude_term: str | None = None
    ) -> pd.DataFrame:
        """
        Match taxids from the classification output to their respective assemblies.
        Tracks failed taxids for later retry or debugging.
        """
        if not os.path.exists(classification_output_path):
            raise FileNotFoundError(f"Classification output file not found: {classification_output_path}")

        df = pd.read_csv(classification_output_path, sep="\t", header=0)

        ids = detect_id_columns(df)
        if ids["taxid_col"]:
            df = rename_columns_to_standard(df, taxid_col=ids["taxid_col"], accid_col=ids["accid_col"])
        else:
            raise ValueError(
                "The classification output file must contain a taxonomic ID column [taxid, taxID or taxon]."
            )

        if not ids["taxid_col"] and not ids["accid_col"]:
            raise ValueError(
                "The classification output file must contain a taxonomic ID column "
                "[taxid, taxID or taxon] or an accession column "
                "[assembly_accession, accession, accID or accid]."
            )

        failed_taxids = []

        rate_limiter = RateLimiter(delay_between_calls=0.5)

        for index, row in df.iterrows():
            rate_limiter.wait()
            taxid = str(int(row["taxid"])) if ids["taxid_col"] and pd.notna(row.get("taxid")) else None

            accession = str(row["accid"]) if ids["accid_col"] and pd.notna(row.get("accid")) else None
            if taxid is None and accession is None:
                self.logger.warning(f"Skipping row {index} due to missing taxid and accession.")
                continue

            try:
                self.logger.info(f"Processing taxid {taxid}...")
                reference = None
                description = None
                if "description" in row and row["description"] is not None:
                    description = str(row["description"])
                if "nucleotide_id" in row and "assembly_id" in row:
                    if row["nucleotide_id"] is not None and not pd.isna(row["nucleotide_id"]):
                        reference = ReferenceData(
                            taxid=taxid,
                            accession=accession,
                            nucleotide_id=str(int(row["nucleotide_id"])),
                            assembly_id=None,
                        )
                    elif row["assembly_id"] is not None and not pd.isna(row["assembly_id"]):
                        reference = ReferenceData(
                            taxid=taxid,
                            accession=accession,
                            nucleotide_id=None,
                            assembly_id=str(int(row["assembly_id"])),
                        )
                passport = Passport(taxid=taxid, accession=accession)
                local_assembly = self.retrieve_assembly(
                    passport, reference_data=reference, include_term=include_term, exclude_term=exclude_term
                )

                if local_assembly:
                    df.at[index, "assembly_accession"] = local_assembly.accession
                    df.at[index, "assembly_file"] = local_assembly.file_path
                else:
                    self.logger.warning(f"No assembly found for taxid {taxid} and accession {accession}.")
                    failed_taxids.append({"taxid": taxid, "accession": accession, "error": "No assembly found"})
                    df.at[index, "assembly_accession"] = None
                    df.at[index, "assembly_file"] = None
            except Exception as e:
                self.logger.error(f"Error processing taxid {taxid}: {e}")
                failed_taxids.append({"taxid": taxid, "accession": accession, "error": str(e)})
                df.at[index, "assembly_accession"] = None
                df.at[index, "assembly_file"] = None

        # Save failed taxids for debugging and manual retry
        self.last_failed_taxids = pd.DataFrame(failed_taxids, columns=["taxid", "accession", "error"])
        n_rows = len(df)
        if failed_taxids:
            failed_file = os.path.join(self.store_path, "failed_taxids.tsv")
            self.last_failed_taxids.to_csv(failed_file, sep="\t", index=False)
            self.logger.warning(
                f"Assembly matching: {len(failed_taxids)}/{n_rows} taxids unmatched "
                f"({len(failed_taxids) / n_rows:.1%}); saved to {failed_file}"
            )
        else:
            self.logger.info(f"Assembly matching: all {n_rows} taxids matched")

        if df.empty:
            df = pd.DataFrame(columns=["taxid", "assembly_accession", "assembly_file"])

        return df

    def setup_mapping_references(
        self, classification_output_path: pd.DataFrame, mapping_references_dir: str = "references_to_map"
    ):

        if not os.path.exists(mapping_references_dir):
            os.makedirs(mapping_references_dir)

        if (
            "assembly_accession" not in classification_output_path.columns
            or "assembly_file" not in classification_output_path.columns
        ):
            print("The DataFrame must contain 'assembly_accession' and 'assembly_file' columns.")

        for _, row in classification_output_path.iterrows():
            accession = row["assembly_accession"]
            assembly_file = row["assembly_file"]

            if pd.isna(accession) or pd.isna(assembly_file):
                self.logger.warning(f"Skipping taxid {row['taxid']} due to missing assembly data.")
                continue

            dest_filename = f"{accession}.fna.gz"
            if "taxid" in row and not pd.isna(row["taxid"]):
                dest_filename = f"{row['taxid']}_{dest_filename}"
            dest_path = os.path.join(mapping_references_dir, dest_filename)

            if not os.path.exists(dest_path):
                self.logger.info(f"Copying {assembly_file} to {dest_path}")
                os.system(f"cp {assembly_file} {dest_path}")
            else:
                self.logger.warning(f"File {dest_path} already exists, skipping copy.")

        self.logger.info(f"Mapping references setup complete in {mapping_references_dir}.")
