/*
 * Reference-truth mapping workflow (stripped-down, no classification).
 *
 * Maps a dataset's simulated reads against the references that were used to
 * simulate it. The reference set is reconstructed directly from the study
 * input table with the same retrieval used by simulate.nf (no include/exclude
 * term) and the reads are mapped with the same mapping / coverage parameters
 * as the production classification pipeline (classify.nf). Clustering and
 * clade matching are intentionally omitted.
 *
 * Outputs are written inside each dataset directory:
 *   $output_dir/$analysis_id/mapping/
 *     - matched_assemblies.tsv
 *     - merged_coverage_statistics.tsv
 *     - mapping_statistics.tsv                 (one row per dataset)
 *     - mapping_statistics_per_reference.tsv   (one row per reference)
 */

params.analysis_id = params.analysis_id ?: ""
params.reads = params.reads ?: ""
params.input_table = params.input_table ?: ""
params.output_dir = params.output_dir ?: "output"

workflow {

    if (params.input_table == "") {
        error("Input table path is not provided. Please set the 'input_table' parameter.")
    }
    if (params.reads == "") {
        error("Reads directory path is not provided. Please set the 'reads' parameter.")
    }

    analysis_id = params.analysis_id

    input_table_ch = Channel
        .fromPath(params.input_table)
        .ifEmpty { error("Cannot find the input table: ${params.input_table}") }

    // Read paired-end simulated reads
    Channel
        .fromPath(["${params.reads}/*1.fq.gz", "${params.reads}/*1.fastq.gz"])
        .map { file -> tuple(file) }
        .ifEmpty { error("Cannot find any paired-end fastq files in ${params.reads}") }
        .set { reads_ch1 }

    Channel
        .fromPath(["${params.reads}/*2.fq.gz", "${params.reads}/*2.fastq.gz"])
        .map { file -> tuple(file) }
        .ifEmpty { error("Cannot find R2 fastq files in ${params.reads}") }
        .set { reads_ch2 }

    reads_ch = reads_ch1.combine(reads_ch2)

    // Reconstruct the references used to simulate this dataset.
    // Same retrieval as simulate.nf: no include_term / exclude_term.
    reference_sequences_ch = ExtractReferenceSequences(analysis_id, input_table_ch)

    flattened_reference_sequences_ch = reference_sequences_ch.reference_sequences.flatMap { ref_list -> ref_list }

    combined_ch = reads_ch.combine(flattened_reference_sequences_ch)

    mapped_reads_ch = MapMinimap2Paired(analysis_id, combined_ch)
    filtered_alignments_ch = FilterBamMsamtools(mapped_reads_ch)
    sorted_reads_ch = sortBam(filtered_alignments_ch)
    coverage_ch = SamtoolsCoverage(sorted_reads_ch)

    coverage_ch = coverage_ch.map { file1, file2, _filename, _refname, group -> tuple(group, file1, file2) }
    coverage_ch = coverage_ch.groupTuple()

    merged_coverage_ch = MergeCoverageStatistics(coverage_ch)

    MappingStatistics(
        analysis_id,
        input_table_ch,
        reference_sequences_ch.matched_assemblies,
        merged_coverage_ch.merged_coverage_statistics,
    )
}


/*
 * Extract the references used to simulate a dataset.
 * Invocation mirrors simulate.nf's ExtractFastaSequences (no include/exclude).
 */
process ExtractReferenceSequences {
    tag "ExtractReferenceSequences ${analysis_id}"

    input:
    val analysis_id
    path classifier_output

    output:
    path "reference_sequences/*gz", emit: reference_sequences, optional: true
    path "reference_sequences/matched_assemblies.tsv", emit: matched_assemblies

    script:
    """
    ${params.python_bin} ${params.references_extract_script} retrieve \
    --input_table ${classifier_output} \
    --assembly_store "${params.assembly_store}" \
    --mapping_references_dir "reference_sequences"
    """
}


/*
 * Map reads to references using minimap2
 */
process MapMinimap2Paired {
    tag "MapMinimap2Paired ${query_id}"

    input:
    val query_id
    tuple path(fastq1), path(fastq2), path(reference)

    output:
    tuple path("${query_id}_${reference.baseName}.bam"), val(query_id), val(reference.baseName)

    script:
    """
    minimap2 -t ${task.cpus} ${params.minimap2_illumina_params} -a ${reference} ${fastq1} ${fastq2} | \
        samtools view -bS -F 4 - > ${query_id}_${reference.baseName}.bam
    """
}


/*
 * Filter BAM files
 */
process FilterBamMsamtools {
    tag "FilterBamMsamtools ${bamfile.baseName}"

    input:
    tuple path(bamfile), val(sample_id), val(reference_id)

    output:
    tuple path("${bamfile.baseName}.filtered.bam"), val(sample_id), val(reference_id)

    script:
    """
    msamtools filter -b ${params.msamtools_params} ${bamfile} > ${bamfile.baseName}.filtered.bam
    """
}


/*
 * sort and index bam file, maintain tuple file, query_id, reference_id in channel
 */
process sortBam {
    tag "sortMapping"

    input:
    tuple path(bamfile), val(query_id), val(reference_id)

    output:
    tuple path("${query_id}_${reference_id}.sorted.bam"), path("${query_id}_${reference_id}.sorted.bam.bai"), val(query_id), val(reference_id)

    script:
    """
    samtools sort ${bamfile} > ${query_id}_${reference_id}.sorted.bam
    samtools index ${query_id}_${reference_id}.sorted.bam
    samtools addreplacerg -r "ID:${query_id}" -r "SM:${query_id}" -o named.bam ${query_id}_${reference_id}.sorted.bam
    mv named.bam ${query_id}_${reference_id}.sorted.bam
    """
}


/*
 * Use samtools to extract BAM file coverage statistics
 */
process SamtoolsCoverage {

    tag "SamtoolsCoverage ${bamfile.baseName}"

    input:
    tuple path(bamfile), path(bamindex), val(query_id), val(reference_id)

    output:
    tuple path("${bamfile.baseName}.stats.txt"), path("${bamfile.baseName}.coverage.txt"), val(bamfile.baseName), val(reference_id), val(query_id)

    script:
    """
    samtools coverage -o ${bamfile.baseName}.coverage.txt ${bamfile}
    samtools stats ${bamfile} | grep ^SN | cut -f 2- > ${bamfile.baseName}.stats.txt
    """
}


/*
 * merge coverage statistics from different BAM files
 */
process MergeCoverageStatistics {
    tag "MergeCoverageStatistics ${query_id}"

    input:
    tuple val(query_id), path(stats_files), path(coverage_files)

    output:
    path "merged_coverage_statistics.tsv", emit: merged_coverage_statistics

    script:
    def coverage_files_string = coverage_files.collect { it[0] }.join(',')
    def stats_files_string = stats_files.collect { it[0] }.join(',')
    """
    #!/usr/bin/env python3
    import pandas as pd
    import glob


    cov_files = "${coverage_files_string}".split(',')
    stats_files = "${stats_files_string}".split(',')
    coverage_data = []
    for ix, file in enumerate(cov_files):
        stats_filename = stats_files[ix]
        stats_df = pd.read_csv(stats_filename, sep="\\t", header=None, names=["stat", "value", "comment"])

        coverage_df = pd.read_csv(file, sep="\\t")
        coverage_df['file'] = file.split('/')[-1]
        coverage_df['error_rate'] = stats_df.loc[stats_df['stat'] == 'error rate:', 'value'].values[0]
        coverage_df['sample'] = stats_filename.replace("_coverage.txt", "").split("/")[-1]

        coverage_data.append(coverage_df)

    merged_df = pd.concat(coverage_data)
    merged_df.to_csv("merged_coverage_statistics.tsv", sep="\\t", index=False)
    """
}


/*
 * Compute mapping statistics for the dataset.
 * The two input tables are copied to the publish directory alongside the
 * derived statistics so everything needed for recall analysis is colocated.
 */
process MappingStatistics {
    tag "MappingStatistics ${analysis_id}"
    publishDir { "${params.output_dir}/${analysis_id}/mapping" }, mode: 'copy'

    input:
    val analysis_id
    path input_table
    path matched_assemblies
    path coverage_statistics

    output:
    path "matched_assemblies.tsv"
    path "merged_coverage_statistics.tsv"
    path "mapping_statistics.tsv"
    path "mapping_statistics_per_reference.tsv"

    script:
    def stats_script = params.mapping_stats_script ?: "${projectDir}/mapping_statistics.py"
    """
    ${params.python_bin} ${stats_script} \
        --input-table ${input_table} \
        --matched-assemblies matched_assemblies.tsv \
        --coverage-statistics merged_coverage_statistics.tsv \
        --dataset ${analysis_id} \
        --output-dir .
    """
}
