import logging
import os
import sys
import pandas as pd
import time
from pathlib import Path
from typing import Optional

from metagenomics_utils.ncbi_tools import Passport, NCBITools, LocalAssembly, ReferenceData
from metagenomics_utils.dataframe_utils import detect_id_columns, rename_columns_to_standard


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

def dl_file(url: str, dest: str) -> None:
    """
    Download a file from a given URL to a specified destination.
    """
    import requests
    response = requests.get(url)
    if response.status_code == 200:
        with open(dest, 'wb') as f:
            f.write(response.content)
        print(f"Downloaded {url} to {dest}")
        return True
    else:
        print(f"Failed to download {url}: {response.status_code}")
        return False 

def dl_file_wget(url: str, dest: str) -> None:
    """
    Download a file using wget.
    """
    import os
    command = f"wget -O {dest} {url}"
    os.system(command)
    print(f"Downloaded {url} to {dest}")

    if not os.path.exists(dest):
        print(f"Failed to download {url} to {dest}")
        return False
    return True


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
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)
        self.logger.propagate = False

        self.ncbi = NCBITools()


    def get_assembly_path(self, taxid: str) -> str:
        return os.path.join(self.store_path, taxid)


    def retrieve_local_assembly(self, passport: Passport) -> Optional[LocalAssembly]:
        """
        Check if the assembly for the given taxid exists in the assembly store.
        """
        taxid_subdir = f"{self.store_path}/{passport.taxid}"

        # Assuming the first file is the assembly file
        assembly_file = os.path.join(taxid_subdir, f"{passport.prefix}_sequence.fasta.gz")

        if not os.path.exists(assembly_file):
            self.logger.warning(f"No assembly file found for taxid {passport.taxid} and accession {passport.accession}")
            return None
        accid = passport.accession

        return LocalAssembly(taxid=passport.taxid, accession=accid, file_path=assembly_file) if assembly_file else None


    def retrieve_assembly(self, passport: Passport, reference_data: Optional[ReferenceData] = None, include_term: Optional[str] = None, exclude_term: Optional[str] = None) -> Optional[LocalAssembly]:
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
            reference_data = self.ncbi.query_sequence_databases(passport, include_term=include_term, exclude_term=exclude_term)
        assembly_dir = os.path.join(self.store_path, str(passport.taxid))
        os.makedirs(assembly_dir, exist_ok=True)

        assembly_file_path = os.path.join(assembly_dir, f"{reference_data.prefix}_sequence.fasta.gz")
        success_dl = self.ncbi.retrieve_sequence_databases(reference_data, assembly_file_path, gzipped=True)

        if not success_dl:
            self.logger.error(f"Failed to download assembly for passport {passport.taxid}")
            return None
        
        return LocalAssembly(taxid=passport.taxid, accession=reference_data.accession, file_path=assembly_file_path)

    def match_taxid_to_assembly(
        self,
        classification_output_path: str,
        include_term: Optional[str] = None,
        exclude_term: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Match taxids from the classification output to their respective assemblies.
        Tracks failed taxids for later retry or debugging.
        """
        if not os.path.exists(classification_output_path):
            raise FileNotFoundError(f"Classification output file not found: {classification_output_path}")

        df = pd.read_csv(classification_output_path, sep='\t', header=0)
        
        ids = detect_id_columns(df)
        if ids['taxid_col']:
            df = rename_columns_to_standard(df, taxid_col=ids['taxid_col'], accid_col=ids['accid_col'])
        else:
            raise ValueError(
                "The classification output file must contain a taxonomic ID column "
                "[taxid, taxID or taxon]."
            )

        if not ids['taxid_col'] and not ids['accid_col']:
            raise ValueError(
                "The classification output file must contain a taxonomic ID column "
                "[taxid, taxID or taxon] or an accession column "
                "[assembly_accession, accession, accID or accid]."
            )

        failed_taxids = []
        
        rate_limiter = RateLimiter(delay_between_calls=0.5)
        
        for index, row in df.iterrows():
            rate_limiter.wait()
            taxid = str(int(row['taxid'])) if ids['taxid_col'] and pd.notna(row.get('taxid')) else None
            accession = str(row['accid']) if ids['accid_col'] and pd.notna(row.get('accid')) else None
            if taxid is None and accession is None:
                self.logger.warning(f"Skipping row {index} due to missing taxid and accession.")
                continue
            
            try:
                self.logger.info(f"Processing taxid {taxid}...")
                reference = None
                description = None
                if "description" in row and row['description'] is not None:
                    description = str(row['description'])
                if 'nucleotide_id' in row and 'assembly_id' in row:
                    if row['nucleotide_id'] is not None and not pd.isna(row['nucleotide_id']):
                        reference = ReferenceData(
                            taxid=taxid,
                            accession=accession,
                            nucleotide_id=str(int(row['nucleotide_id'])),
                            assembly_id=None
                        )
                    elif row['assembly_id'] is not None and not pd.isna(row['assembly_id']):
                        reference = ReferenceData(
                            taxid=taxid,
                            accession=accession,
                            nucleotide_id=None,
                            assembly_id=str(int(row['assembly_id']))
                        )
                passport = Passport(taxid=taxid, accession=accession)
                local_assembly = self.retrieve_assembly(passport, reference_data=reference, include_term=include_term, exclude_term=exclude_term)

                if local_assembly:
                    df.at[index, 'assembly_accession'] = local_assembly.accession
                    df.at[index, 'assembly_file'] = local_assembly.file_path
                else:
                    self.logger.warning(f"No assembly found for taxid {taxid} and accession {accession}.")
                    failed_taxids.append({'taxid': taxid, 'accession': accession, 'error': 'No assembly found'})
                    df.at[index, 'assembly_accession'] = None
                    df.at[index, 'assembly_file'] = None
            except Exception as e:
                self.logger.error(f"Error processing taxid {taxid}: {e}")
                failed_taxids.append({'taxid': taxid, 'accession': accession, 'error': str(e)})
                df.at[index, 'assembly_accession'] = None
                df.at[index, 'assembly_file'] = None
        
        # Save failed taxids for debugging and manual retry
        if failed_taxids:
            failed_df = pd.DataFrame(failed_taxids)
            failed_file = os.path.join(self.store_path, 'failed_taxids.tsv')
            failed_df.to_csv(failed_file, sep='\t', index=False)
            self.logger.warning(f"Saved {len(failed_taxids)} failed taxids to {failed_file}")
        
        if df.empty:
            df = pd.DataFrame(columns=['taxid', 'assembly_accession', 'assembly_file'])

        return df

    
    def setup_mapping_references(
        self,
        classification_output_path: pd.DataFrame,
        mapping_references_dir: str = "references_to_map"):

        if not os.path.exists(mapping_references_dir):
            os.makedirs(mapping_references_dir)
        
        if "assembly_accession" not in classification_output_path.columns or \
        "assembly_file" not in classification_output_path.columns:
            print("The DataFrame must contain 'assembly_accession' and 'assembly_file' columns.")
        
        for _, row in classification_output_path.iterrows():
            accession = row['assembly_accession']
            assembly_file = row['assembly_file']

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


