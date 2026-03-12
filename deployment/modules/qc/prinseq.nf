/*
 * Quality control of paired-end reads using prinseq++
 */
process QCReadsPrinseqPaired {
    debug true

    input:
    val sample_id
    tuple path(fastq1), path(fastq2)
    val prin output:
    tupleseq_params

    path("${sample_id}_good_R1.fastq.gz"), path("${sample_id}_good_R2.fastq.gz"), emit: reads

    script:
    """
    prinseq++ ${prinseq_params} -fastq ${fastq1} -fastq2 ${fastq2} \
        -out_good ${sample_id}_good_R1.fastq -out_bad ${sample_id}_bad_R1.fastq \
        -out_good2 ${sample_id}_good_R2.fastq -out_bad2 ${sample_id}_bad_R2.fastq
    bgzip ${sample_id}_good_R1.fastq && bgzip ${sample_id}_good_R2.fastq
    bgzip ${sample_id}_bad_R1.fastq && bgzip ${sample_id}_bad_R2.fastq
    """
}

/*
 * Quality control of single-end reads using prinseq++
 */
process QCReadsPrinseqSingle {
    debug true

    input:
    val sample_id
    path fastq1
    val prinseq_params

    output:
    tuple path("${sample_id}_good_R1.fastq.gz")

    script:
    """
    prinseq++ ${prinseq_params} -fastq ${fastq1} -out_good ${sample_id}_good_R1.fastq -out_bad ${sample_id}_bad_R1.fastq
    bgzip ${sample_id}_good_R1.fastq && bgzip ${sample_id}_bad_R1.fastq
    """
}
