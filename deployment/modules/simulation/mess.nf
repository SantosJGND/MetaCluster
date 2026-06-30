/*
 * Simulate reads using the mess package
 */
process SimulateReadsMess {
    tag "SimulateReadsMess ${sample_id}"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    val sample_id
    path input_table
    val technology

    output:
    tuple val("${sample_id}"), path("${sample_id}/fastq/*.fq.gz")

    script:
    """
    mess simulate --input ${input_table} \
        --output "${sample_id}" \
        --threads 3 \
        --tech ${technology} \
        --bam
    """
}
