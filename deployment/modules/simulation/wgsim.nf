/*
 * Simulate reads using wgsim
 */
process WgsimSimulateReads {
    tag "WgsimSimulateReads ${sample_id}"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    val sample_id
    path input_table
    val wgsim_python_path
    val wgsim_args
    val output_prefix

    output:
    tuple val("${sample_id}"), path("${sample_id}/fastq/*.fq.gz")

    script:
    """
    ${wgsim_python_path} \
        --input_table ${input_table} \
        --prefix "${sample_id}" \
        --output_dir "${sample_id}/fastq" \
        --wgsim_args "${wgsim_args}"

    bgzip ${sample_id}/fastq/${sample_id}_R1.fq
    bgzip ${sample_id}/fastq/${sample_id}_R2.fq
    """
}
