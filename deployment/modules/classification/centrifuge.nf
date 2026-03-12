/*
 * Classify reads using Centrifuge in paired-end mode
 */
process CentrifugeClassificationPaired {
    tag "CentrifugeClassificationPaired ${sample_id}"
    publishDir "${params.output_dir}/${sample_id}/classification/centrifuge", mode: 'copy'

    input:
    val sample_id
    tuple path(fastq1), path(fastq2)
    val centrifuge_index
    val centrifuge_params
    val classifier_process_script
    val python_bin
    val minimum_uniq_reads

    output:
    path "${sample_id}_centrifuge_report.txt"
    path "${sample_id}_centrifuge_processed_classifier_output.tsv"
    val "${sample_id}"

    script:
    def fastq2_arg = fastq2 != 'NO_FILE' && fastq2.exists() ? "-2 ${fastq2}" : ""
    """
    if [ -f ${fastq1} ] && [ -f ${fastq2} ]; then
        centrifuge -x ${centrifuge_index} -1 ${fastq1} -2 ${fastq2} \
            -S ${sample_id}_centrifuge_classification.tsv \
            --output ${sample_id}_centrifuge_classification.txt \
            --report-file ${sample_id}_centrifuge_report.txt \
            ${centrifuge_params}
    elif [ -f ${fastq1} ]; then
        centrifuge -x ${centrifuge_index} -U ${fastq1} \
            -S ${sample_id}_centrifuge_classification.tsv \
            --output ${sample_id}_centrifuge_classification.txt \
            --report-file ${sample_id}_centrifuge_report.txt \
            ${centrifuge_params}
    fi

    ${python_bin} ${classifier_process_script} \
        --input "${sample_id}_centrifuge_report.txt" \
        --output "${sample_id}_centrifuge_processed_classifier_output.tsv" \
        --type "centrifuge" \
        --nuniq_threshold ${minimum_uniq_reads}
    """
}
