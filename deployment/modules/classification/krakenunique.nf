/*
 * Classify reads using KrakenUnique in paired-end mode
 */
process KrakenUniqueClassificationPaired {
    tag "KrakenUniqueClassificationPaired ${sample_id}"
    publishDir "${params.output_dir}/${sample_id}/classification/krakenunique", mode: 'copy'

    input:
    val sample_id
    tuple path(fastq1), path(fastq2)
    val krakenunique_bin
    val krakenunique_index
    val krakenunique_params
    val classifier_process_script
    val python_bin
    val minimum_uniq_reads

    output:
    path "${sample_id}_krakenunique_classification.txt"
    path "${sample_id}_krakenunique_processed_classifier_output.tsv"
    val "${sample_id}"

    script:
    def fastq2_arg = fastq2 != 'NO_FILE' && fastq2.exists() ? "${fastq2}" : ""
    """
    if [ -f ${fastq2} ]; then
        ${krakenunique_bin} --db ${krakenunique_index} \
            --report ${sample_id}_krakenunique_report.txt \
            --output ${sample_id}_krakenunique_classification.txt \
            ${fastq1} ${fastq2} ${krakenunique_params}
    else
        ${krakenunique_bin} --db ${krakenunique_index} \
            --report ${sample_id}_krakenunique_report.txt \
            --output ${sample_id}_krakenunique_classification.txt \
            ${fastq1} ${krakenunique_params}
    fi

    ${python_bin} ${classifier_process_script} \
        --input "${sample_id}_krakenunique_classification.txt" \
        --output "${sample_id}_krakenunique_processed_classifier_output.tsv" \
        --type "kuniq" \
        --nuniq_threshold ${minimum_uniq_reads}
    """
}
