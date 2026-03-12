from metagenomics_utils.ncbi_tools import Passport, NCBITools
from metagenomics_utils.reference_utils import AssemblyStore
import os
import pandas as pd
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "metagenomics_utils"))
from metagenomics_utils.dataframe_utils import detect_id_columns, rename_columns_to_standard, validate_required_columns

import argparse

def get_args():
    """
    Define the argument parser with subcommands.
    """
    parser = argparse.ArgumentParser(description="Manage taxid-to-assembly retrieval and reference setup.")
    subparsers = parser.add_subparsers(dest="command", required=True, help="Subcommands: retrieve or check")

    # Subcommand: retrieve
    retrieve_parser = subparsers.add_parser("retrieve", help="Retrieve assemblies based on the input table.")
    retrieve_parser.add_argument(
        "--input_table",
        type=str,
        required=True,
        help="Path to the classification output file."
    )
    retrieve_parser.add_argument(
        "--assembly_store",
        type=str,
        default="assemblies",
        help="Directory to store downloaded assemblies."
    )
    retrieve_parser.add_argument(
        "--mapping_references_dir",
        type=str,
        default="references_to_map",
        help="Directory to store mapping references."
    )

    retrieve_parser.add_argument(
        "--include_term",
        type=str,
        default=None,
        help="Term to include in NCBI search."
    )
    retrieve_parser.add_argument(
        "--exclude_term",
        type=str,
        default=None,
        help="Term to exclude from NCBI search."
    )


    # Subcommand: check
    check_parser = subparsers.add_parser("check", help="Check if mapping ids can be retrieved.")
    check_parser.add_argument(
        "--input_table",
        type=str,
        required=True,
        help="Path to the classification output file."
    )

    check_parser.add_argument(
        "--assessment",
        type=str,
        default="assembly_assessment.tsv",
        help="Path to the assessment file to check assemblies."
    )
    check_parser.add_argument(
        "--include_term",
        type=str,
        default=None,
        help="Term to include in NCBI search."
    )
    check_parser.add_argument(
        "--exclude_term",
        type=str,
        default=None,
        help="Term to exclude from NCBI search."
    )

    return parser.parse_args()


def retrieve_assemblies(args):
    """
    Retrieve assemblies based on the input table and store them in the specified directory.
    """
    classification_output_path = args.input_table
    assembly_store = args.assembly_store
    mapping_references_dir = args.mapping_references_dir

    assembly_store = AssemblyStore(assembly_store)
    df = assembly_store.match_taxid_to_assembly(classification_output_path)
    
    assembly_store.setup_mapping_references(df, mapping_references_dir=mapping_references_dir)
    if "assembly_accession" not in df.columns or "assembly_file" not in df.columns:
        df['assembly_accession'] = None
        df['assembly_file'] = None
    df = df.dropna(subset=['assembly_accession', 'assembly_file'])
    df.to_csv(
        os.path.join(mapping_references_dir, "matched_assemblies.tsv"),
        index=False,
        sep='\t'
    )

def check_assemblies_exist(args):
    """
    Check if the assemblies can be retrieved from the input table.
    """
    input_table = args.input_table
    include_term = args.include_term if hasattr(args, 'include_term') else None
    exclude_term = args.exclude_term if hasattr(args, 'exclude_term') else None
    df = pd.read_csv(input_table, sep='\t')

    ids = detect_id_columns(df)
    if ids['taxid_col'] or ids['accid_col']:
        df = rename_columns_to_standard(df, taxid_col=ids['taxid_col'], accid_col=ids['accid_col'])
    else:
        raise ValueError(
            "The classification output file must contain a taxonomic ID column "
            "[taxid, taxID or taxon] or an accession column "
            "[assembly_accession, accession, accID or accid]."
        )

    def check_assembly_exists(row, include_term=None, exclude_term=None):
        """
        Check if the assembly for the given taxid exists.
        """
        taxid = str(int(row["taxid"])) if "taxid" in row and pd.notna(row["taxid"]) else None
        accid = str(row["accid"]) if "accid" in row and pd.notna(row["accid"]) else None

        passport = Passport(taxid=taxid, accession=accid)
        ncbi_tools = NCBITools()

        reference_data = ncbi_tools.query_sequence_databases(passport, include_term=include_term, exclude_term=exclude_term)

        row['assembly_accession'] = reference_data.accession
        row['description'] = reference_data.description
        row['nucleotide_id'] = reference_data.nucleotide_id
        row['assembly_id'] = reference_data.assembly_id
        row['lineage'] = reference_data.lineage

        return row

    df = df.apply(lambda row: check_assembly_exists(row, include_term=include_term, exclude_term=exclude_term), axis=1)
    df.to_csv(args.assessment, index=False, sep='\t')

def main():

    args = get_args()

    if args.command == "retrieve":
        retrieve_assemblies(args)
    elif args.command == "check":
        check_assemblies_exist(args)



if __name__ == "__main__":
    main()