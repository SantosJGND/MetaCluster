/*
 * Map paired-end reads to reference using minimap2
 * Input: tuple(sample_id, fastq1, fastq2, reference)
 */
process MapMinimap2Paired {
    tag "MapMinimap2Paired ${sample_id} ${reference.baseName}"

    input:
    tuple val(sample_id), path(fastq1), path(fastq2), path(reference)
    val minimap2_params

    output:
    tuple path("${sample_id}_${reference.baseName}.bam"), val(sample_id), val(reference.baseName)

    script:
    def fastq2_arg = fastq2 != 'NO_FILE' && fastq2.exists() ? "${fastq2}" : ""
    def map_mode = fastq2_arg ? "-ax sr" : "-ax sr"
    """
    mkdir -p ${sample_id}
    minimap2 ${minimap2_params} ${map_mode} ${reference} ${fastq1} ${fastq2_arg} | samtools view -bS -F 4 - > ${sample_id}_${reference.baseName}.bam
    """
}
