/*
 * Sort and index BAM file
 * Input: tuple(bamfile, sample_id, reference_id)
 * Output: tuple(sorted_bam, bai, sample_id, reference_id)
 */
process SortBam {
    tag "SortBam ${sample_id}_${reference_id}"

    input:
    tuple path(bamfile), val(sample_id), val(reference_id)

    output:
    tuple path("${sample_id}_${reference_id}.sorted.bam"), path("${sample_id}_${reference_id}.sorted.bam.bai"), val(sample_id), val(reference_id)

    script:
    """
    samtools sort ${bamfile} > ${sample_id}_${reference_id}.sorted.bam
    samtools index ${sample_id}_${reference_id}.sorted.bam
    samtools addreplacerg -r "ID:${sample_id}" -r "SM:${sample_id}" -o named.bam ${sample_id}_${reference_id}.sorted.bam
    mv named.bam ${sample_id}_${reference_id}.sorted.bam
    """
}
