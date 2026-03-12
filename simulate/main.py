import argparse
import os
import pandas as pd


def get_args():
    """
    Define the argument parser with subcommands.
    """
    parser = argparse.ArgumentParser(description="Simulate reads files using wgsim from reference genomes provided in a table.")
    parser.add_argument(
        "--input_table",
        type=str,
        required=True,
        help="Path to the classification output file."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Path to the output directory."
    )

    parser.add_argument(
        "--prefix",
        type=str,
        default="fastq",
        help="Prefix for the output files."
    )

    parser.add_argument(
        "--wgsim_args",
        type=str,
        default="-e 0.02 -1 150 -2 150 -r 0.02 -R 0.05 -S 100",
        help="Arguments to pass to wgsim."
    ) 

    return parser.parse_args()


def simulate_wgsim_illumina(reference_genome, output_prefix, wgsim_args):
    """
    Simulate reads using wgsim.
    """
    log_run_file = f"{output_prefix}_wgsim.log"
    log_plan_file = f"{output_prefix}_wgsim_plan.log"
    command = f"wgsim {wgsim_args} {reference_genome} {output_prefix}_1.fq {output_prefix}_2.fq > {log_plan_file} 2> {log_run_file}"

    os.system(command)
    
    if os.path.exists(f"{output_prefix}_1.fq") and os.path.exists(f"{output_prefix}_2.fq"):
        return True
    else:
        return False

def concatenate_fastq_files(file_list, output_file):
    """
    Concatenate multiple FASTQ files into a single file.
    """
    with open(output_file, 'wb') as wfd:
        for f in file_list:
            with open(f, 'rb') as fd:
                wfd.write(fd.read())
    return os.path.exists(output_file)

def main():
    args = get_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Read the input table
    df = pd.read_csv(args.input_table, sep = "\t")

    all_read1_files = []
    all_read2_files = []

    original_args = args.wgsim_args
    
    for idx, row in df.iterrows():
        genome_id = idx
        reference_genome = row['path']
        nreads = row['reads']
        
        output_prefix = os.path.join(args.output_dir, "fastq" + f"_{genome_id}")
        
        mut_rate= None
        wgsim_args = original_args
        
        if "mutation_rate" in row and not pd.isna(row["mutation_rate"]):
            mut_rate = row["mutation_rate"]
            if "-r" in wgsim_args:
                wgsim_args = wgsim_args.split("-r")
                wgsim_args_part_1 = wgsim_args[1].split(' ')
                wgsim_args_part_1 = [arg for arg in wgsim_args_part_1 if arg != '']
                wgsim_args_part_1 = ' '.join(wgsim_args_part_1[1:])

                wgsim_args[1] = f"{mut_rate} {wgsim_args_part_1}"
                wgsim_args = "-r ".join(wgsim_args)

            else:
                wgsim_args = f"-r {mut_rate} {wgsim_args}"

        wgsim_args = wgsim_args.replace("  ", " ").strip()
        wgsim_args = wgsim_args + " -N " + str(nreads)
        print(f"Simulating {nreads} reads for genome {genome_id} with mutation rate {mut_rate if mut_rate else 'default'}...")

        success = simulate_wgsim_illumina(reference_genome, output_prefix, wgsim_args)

        if success:
            all_read1_files.append(f"{output_prefix}_1.fq")
            all_read2_files.append(f"{output_prefix}_2.fq")
            print(f"Successfully simulated reads for {genome_id}.")
        else:
            print(f"Failed to simulate reads for {genome_id}.")
    
    
    # Concatenate all read files
    final_read1_file = os.path.join(args.output_dir, f"{args.prefix}_R1.fq")
    final_read2_file = os.path.join(args.output_dir, f"{args.prefix}_R2.fq")

    if concatenate_fastq_files(all_read1_files, final_read1_file):
        print(f"All read 1 files concatenated into {final_read1_file}.")
    else:
        print("Failed to concatenate read 1 files.")
    
    if concatenate_fastq_files(all_read2_files, final_read2_file):
        print(f"All read 2 files concatenated into {final_read2_file}.")
    else:
        print("Failed to concatenate read 2 files.")



if __name__ == "__main__":
    main()