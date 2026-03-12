/*
 * Classify reads using Diamond in paired-end mode
 */
process DiamondClassificationPaired {
    tag "DiamondClassificationPaired ${sample_id}"
    publishDir "${params.output_dir}/${sample_id}/classification/diamond", mode: 'copy'

    input:
    val sample_id
    tuple path(fastq1), path(fastq2)
    val diamond_bin
    val diamond_index
    val diamond_params
    val classifier_process_script
    val python_bin
    val minimum_uniq_reads

    output:
    path "${sample_id}_diamond_classification.tsv"
    path "${sample_id}_diamond_processed_classifier_output.tsv"
    val "${sample_id}"

    script:
    def fastq2_arg = fastq2 != 'NO_FILE' && fastq2.exists() ? "--query ${fastq2}" : ""
    """
    ${diamond_bin} blastx \
        --db ${diamond_index} \
        --query ${fastq1} \
        ${fastq2_arg} \
        --out ${sample_id}_diamond_classification.tsv \
        --outfmt 6 \
        ${diamond_params}

    ${python_bin} ${classifier_process_script} \
        --input "${sample_id}_diamond_classification.tsv" \
        --output "${sample_id}_diamond_processed_classifier_output.tsv" \
        --type "diamond" \
        --nuniq_threshold ${minimum_uniq_reads}
    """
}
