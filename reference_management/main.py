import os
import sys
from pathlib import Path

import pandas as pd

from metagenomics_utils.ncbi_tools import NCBITools, Passport
from metagenomics_utils.reference_utils import AssemblyStore

sys.path.insert(0, str(Path(__file__).parent.parent / "metagenomics_utils"))
import argparse

from metagenomics_utils.dataframe_utils import detect_id_columns, rename_columns_to_standard


def get_args():
    """
    Define the argument parser with subcommands.
    """
    parser = argparse.ArgumentParser(description="Manage taxid-to-assembly retrieval and reference setup.")
    subparsers = parser.add_subparsers(dest="command", required=True, help="Subcommands: retrieve or check")

    # Subcommand: retrieve
    retrieve_parser = subparsers.add_parser("retrieve", help="Retrieve assemblies based on the input table.")
    retrieve_parser.add_argument(
        "--input_table", type=str, required=True, help="Path to the classification output file."
    )
    retrieve_parser.add_argument(
        "--assembly_store", type=str, default="assemblies", help="Directory to store downloaded assemblies."
    )
    retrieve_parser.add_argument(
        "--mapping_references_dir", type=str, default="references_to_map", help="Directory to store mapping references."
    )

    retrieve_parser.add_argument("--include_term", type=str, default=None, help="Term to include in NCBI search.")
    retrieve_parser.add_argument("--exclude_term", type=str, default=None, help="Term to exclude from NCBI search.")

    retrieve_parser.add_argument(
        "--min_uniq_reads",
        type=int,
        default=1,
        help="Minimum uniq_reads for a classified taxid to require a matched assembly (default: 1).",
    )
    retrieve_parser.add_argument(
        "--no_fail_on_missing",
        action="store_true",
        help="Do not exit non-zero when classified taxids (uniq_reads >= --min_uniq_reads) lack a matched assembly.",
    )
    retrieve_parser.add_argument(
        "--max_missing_pct",
        type=float,
        default=5.0,
        help=(
            "Skip the dataset (exit code 3) when the fraction of classified references "
            "(uniq_reads >= --min_uniq_reads) lacking a matched assembly exceeds this "
            "percentage (default: 5.0). Datasets at or below the threshold proceed normally."
        ),
    )

    # Subcommand: check
    check_parser = subparsers.add_parser("check", help="Check if mapping ids can be retrieved.")
    check_parser.add_argument("--input_table", type=str, required=True, help="Path to the classification output file.")

    check_parser.add_argument(
        "--assessment",
        type=str,
        default="assembly_assessment.tsv",
        help="Path to the assessment file to check assemblies.",
    )
    check_parser.add_argument("--include_term", type=str, default=None, help="Term to include in NCBI search.")
    check_parser.add_argument("--exclude_term", type=str, default=None, help="Term to exclude from NCBI search.")

    return parser.parse_args()


def missing_pct_exceeds(n_missing: int, n_qualified: int, max_pct: float = 5.0) -> bool:
    """
    True when more than ``max_pct`` percent of qualified references lack a match.

    A dataset with no qualified references has 0% missing and never exceeds the
    threshold. The threshold is strict: exactly ``max_pct`` proceeds.
    """
    if not n_qualified:
        return False
    return (100.0 * n_missing / n_qualified) > max_pct


def retrieve_assemblies(args):
    """
    Retrieve assemblies based on the input table and store them in the specified directory.
    """
    classification_output_path = args.input_table
    assembly_store = args.assembly_store
    mapping_references_dir = args.mapping_references_dir
    min_uniq_reads = args.min_uniq_reads
    max_missing_pct = args.max_missing_pct
    fail_on_missing = not args.no_fail_on_missing

    assembly_store = AssemblyStore(assembly_store)
    df = assembly_store.match_taxid_to_assembly(classification_output_path)

    assembly_store.setup_mapping_references(df, mapping_references_dir=mapping_references_dir)
    if "assembly_accession" not in df.columns or "assembly_file" not in df.columns:
        df["assembly_accession"] = None
        df["assembly_file"] = None

    unmatched = df[df["assembly_accession"].isna() | df["assembly_file"].isna()]
    if not unmatched.empty:
        unmatched_path = os.path.join(mapping_references_dir, "unmatched_taxids.tsv")
        unmatched.to_csv(unmatched_path, index=False, sep="\t")
        print(
            f"WARNING: {len(unmatched)}/{len(df)} classification taxids have no matched assembly "
            f"(these cannot be recalled via read mapping). Saved to {unmatched_path}"
        )

    # Assembly-completeness gate: a classified taxid with uniq_reads >= min_uniq_reads
    # requires a matched assembly (from the local store or the NCBI fallback). When the
    # share of qualified taxids lacking a match exceeds max_missing_pct%, the dataset is
    # SKIPPED (exit code 3, DATASET_SKIPPED.txt marker); at or below the threshold the
    # extracted subset proceeds normally.
    qualified = df
    if "uniq_reads" in df.columns:
        qualified = qualified[qualified["uniq_reads"] >= min_uniq_reads]

    missing_classified = qualified[qualified["assembly_accession"].isna() | qualified["assembly_file"].isna()]
    n_qualified = len(qualified)
    n_missing = len(missing_classified)
    missing_pct = (100.0 * n_missing / n_qualified) if n_qualified else 0.0
    skip_dataset = fail_on_missing and missing_pct_exceeds(n_missing, n_qualified, max_missing_pct)
    if not missing_classified.empty:
        missing_path = os.path.join(mapping_references_dir, "unmatched_classified_taxids.tsv")
        missing_classified.to_csv(missing_path, index=False, sep="\t")
        action = "SKIP" if skip_dataset else "WARNING"
        print(
            f"{action}: {n_missing}/{n_qualified} classified taxids "
            f"(uniq_reads >= {min_uniq_reads}) have no matched assembly ({missing_pct:.1f}% "
            f"> max {max_missing_pct}%). Saved to {missing_path}"
        )
        if skip_dataset:
            skip_path = os.path.join(mapping_references_dir, "DATASET_SKIPPED.txt")
            with open(skip_path, "w") as fh:
                fh.write(
                    f"Skipped: {n_missing}/{n_qualified} classified references "
                    f"({missing_pct:.1f}%) lack a matched assembly, exceeding "
                    f"max_missing_pct={max_missing_pct}%. Dataset excluded from analysis.\n"
                )
            print(
                f"Assembly validation failed: {missing_pct:.1f}% missing > "
                f"max_missing_pct={max_missing_pct}% -> dataset skipped (excluded from analysis)."
            )
    else:
        print(f"Assembly validation: all {n_qualified} classified taxids matched")

    df.dropna(subset=["assembly_accession", "assembly_file"]).to_csv(
        os.path.join(mapping_references_dir, "matched_assemblies.tsv"), index=False, sep="\t"
    )

    if skip_dataset:
        sys.exit(3)


def check_assemblies_exist(args):
    """
    Check if the assemblies can be retrieved from the input table.
    """
    input_table = args.input_table
    include_term = args.include_term if hasattr(args, "include_term") else None
    exclude_term = args.exclude_term if hasattr(args, "exclude_term") else None
    df = pd.read_csv(input_table, sep="\t")

    ids = detect_id_columns(df)
    if ids["taxid_col"] or ids["accid_col"]:
        df = rename_columns_to_standard(df, taxid_col=ids["taxid_col"], accid_col=ids["accid_col"])
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

        reference_data = ncbi_tools.query_sequence_databases(
            passport, include_term=include_term, exclude_term=exclude_term
        )

        row["assembly_accession"] = reference_data.accession
        row["description"] = reference_data.description
        row["nucleotide_id"] = reference_data.nucleotide_id
        row["assembly_id"] = reference_data.assembly_id
        row["lineage"] = reference_data.lineage

        return row

    df = df.apply(lambda row: check_assembly_exists(row, include_term=include_term, exclude_term=exclude_term), axis=1)
    df.to_csv(args.assessment, index=False, sep="\t")


def main():

    args = get_args()

    if args.command == "retrieve":
        retrieve_assemblies(args)
    elif args.command == "check":
        check_assemblies_exist(args)


if __name__ == "__main__":
    main()
