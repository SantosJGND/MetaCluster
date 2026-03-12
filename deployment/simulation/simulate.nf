params.input_table = params.input_table ?: ""

workflow {

    if (params.input_table == "") {
        error "Input table path is not provided. Please set the 'input_table' parameter."
    }
    input_table_ch = Channel
        .fromPath(params.input_table)
        .ifEmpty { error("Cannot find the input table: ${params.input_table}") }

    // Get reference sequences for each taxonomic ID in the input table
    matched_table_ch = ExtractFastaSequences(input_table_ch)
    input_table_ch = FormatToMess(input_table_ch, matched_table_ch.input_table_with_sequences)

    // Simulate reads based on the input table
    reads_ch = WgsimSimulateReads(params.technology, input_table_ch)
}


/*
* Match FormatSequences file to mess input table
*/
process FormatToMess {
    tag "FormatToMess ${input_table.baseName}"
    publishDir "${params.output_dir}/${input_table.baseName}/input", mode: 'copy'

    input:
    path input_table
    path matched_table

    output:
    path "${input_table.baseName}.tsv", emit: formatted_input_table

    script:
    """
    #!/usr/bin/env python3
    import os
    import pandas as pd
    matched_table = pd.read_csv("${matched_table}", sep="\\t")
    if 'path' not in matched_table.columns:
        if 'assembly_file' not in matched_table.columns:
            raise ValueError("The matched table must contain a 'path' or 'assembly_file' column.")
        matched_table.rename(columns = {'assembly_file': 'path'}, inplace=True)
    matched_table['fasta'] = matched_table['path'].apply(lambda x: os.path.basename(x))
    matched_table.to_csv("${input_table.baseName}.tsv", sep="\\t", index=False)
    """
}


/*
* extract fasta sequences corresponding to the taxonomic IDs in the input_table
*/
process ExtractFastaSequences {
    tag "ExtractFastaSequences ${input_table.baseName}"

    input:
    path input_table

    output:
    path "reference_sequences/matched_assemblies.tsv", emit: input_table_with_sequences

    script:
    """
    ${params.python_bin} ${params.references_extract_script} retrieve \
    --input_table ${input_table} \
    --assembly_store "${params.assembly_store}" \
    --mapping_references_dir "reference_sequences" 
    """
}


/*
* Simulate reads using wgsim
*/
process WgsimSimulateReads {
    tag "WgsimSimulateReads ${input_table.baseName}"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    val technology
    path input_table

    output:
    tuple val("${input_table.baseName}"), path("${input_table.baseName}/fastq/*.fq.gz")

    script:
    """
    ${params.python_bin} ${params.wgsim_python_path} \
    --input_table ${input_table} \
    --prefix "${input_table.baseName}" \
    --output_dir "${input_table.baseName}/fastq" \
    --wgsim_args "${params.wgsim_args}"

    bgzip ${input_table.baseName}/fastq/${input_table.baseName}_R1.fq
    bgzip ${input_table.baseName}/fastq/${input_table.baseName}_R2.fq

    """
}


/*
* Simulate reads using the mess package and conda environment
*/
process SimulateReadsMess {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    val technology
    path input_table

    output:
    tuple val("${input_table.baseName}"), path("${input_table.baseName}/fastq/*.fq.gz")

    script:
    """
    mess simulate --input ${input_table} \
    --output "${input_table.baseName}" \
    --threads 3 \
    --tech ${technology} \
    --bam

    """
}
