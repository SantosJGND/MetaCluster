/*
 * Filter bam file using msamtools
 */
process FilterBamMsamtools {
    tag "FilterBamMsamtools ${sample_id}_${reference_id}"

    input:
    tuple path(bamfile), val(sample_id), val(reference_id)
    val msamtools_params

    output:
    tuple path("${sample_id}_${reference_id}.filtered.bam"), val(sample_id), val(reference_id)

    script:
    """
    msamtools filter -b ${msamtools_params} ${bamfile} > ${sample_id}_${reference_id}.filtered.bam
    """
}
