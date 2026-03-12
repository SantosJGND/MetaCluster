/*
 * Run Kraken2 classification on paired-end reads.
 */
process Kraken2ClassificationPaired {
    tag "Kraken2Classification ${sample_id}"
    publishDir "${params.output_dir}/${sample_id}/classification/kraken2", mode: 'copy'

    input:
    val sample_id
    tuple path(fastq1), path(fastq2)
    val kraken2_bin
    val kraken2_index
    val kraken2_params
    val classifier_process_script
    val python_bin
    val minimum_uniq_reads

    output:
    path "${sample_id}_kraken2_report.txt"
    path "${sample_id}_krk2_processed_classifier_output.tsv"
    val sample_id

    script:
    def fastq2_arg = fastq2 != 'NO_FILE' && fastq2.exists() ? "${fastq1} ${fastq2}" : "${fastq1}"
    """
    ${kraken2_bin} --db ${kraken2_index} \
        --report ${sample_id}_kraken2_report.txt \
        --output ${sample_id}_kraken2_classification.txt \
        ${fastq2_arg} ${kraken2_params}

    ${python_bin} ${classifier_process_script} \
        --input "${sample_id}_kraken2_report.txt" \
        --output "${sample_id}_krk2_processed_classifier_output.tsv" \
        --type "kraken2" \
        --nuniq_threshold ${minimum_uniq_reads}
    """
}

/*
 * Run Kraken2 classification on single-end reads.
 */
process Kraken2ClassificationSingle {
    tag "Kraken2ClassificationSingle ${sample_id}"
    publishDir "${params.output_dir}/${sample_id}/classification/kraken2", mode: 'copy'

    input:
    val sample_id
    path fastq1
    val kraken2_bin
    val kraken2_index
    val kraken2_params
    val classifier_process_script
    val python_bin
    val minimum_uniq_reads

    output:
    path "${sample_id}_kraken2_report.txt"
    path "${sample_id}_krk2_processed_classifier_output.tsv"
    val sample_id

    script:
    """
    ${kraken2_bin} --db ${kraken2_index} \
        --report ${sample_id}_kraken2_report.txt \
        --output ${sample_id}_kraken2_classification.txt \
        ${fastq1} ${kraken2_params}

    ${python_bin} ${classifier_process_script} \
        --input "${sample_id}_kraken2_report.txt" \
        --output "${sample_id}_krk2_processed_classifier_output.tsv" \
        --type "kraken2" \
        --nuniq_threshold ${minimum_uniq_reads}
    """
}
