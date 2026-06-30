include { MergeClassificationResults } from '../modules/classification/merge'
include { Kraken2ClassificationPaired } from '../modules/classification/kraken2'
include { DiamondClassificationPaired } from '../modules/classification/diamond'
include { KrakenUniqueClassificationPaired } from '../modules/classification/krakenunique'
include { CentrifugeClassificationPaired } from '../modules/classification/centrifuge'
include { QCReadsPrinseqPaired } from '../modules/qc/prinseq'
include { MatchCladeReportWithReferenceSequences } from '../modules/evaluation/match'

params.analysis_id = params.analysis_id ?: UUID.randomUUID().toString()
params.reads = params.reads ?: ""
params.output_dir = params.output_dir ?: "output"
params.qc = params.qc ?: true
params.centrifuge = params.centrifuge ?: true
params.kraken2 = params.kraken2 ?: true
params.diamond = params.diamond ?: false
params.krakenunique = params.krakenunique ?: false

workflow {

    analysis_id = params.analysis_id

    if (params.reads == "") {
        error("Reads directory path is not provided. Please set the 'reads' parameter.")
    }

    // Read input reads
    Channel
        .fromPath(["${params.reads}/*1.fq.gz", "${params.reads}/*1.fastq.gz"])
        .map { file -> tuple(file) }
        .ifEmpty { error('Cannot find any paired-end fastq files') }
        .set { reads_ch1 }

    Channel
        .fromPath(["${params.reads}/*2.fq.gz", "${params.reads}/*2.fastq.gz"])
        .map { file -> tuple(file) }
        .ifEmpty { file('NO_FILE') }
        .set { reads_ch2 }

    reads_ch = reads_ch1.combine(reads_ch2)

    // Conditional QC
    if (params.qc) {
        reads_ch = QCReadsPrinseqPaired(analysis_id, reads_ch, params.prinseq_params).reads
    }

    // Conditional classifiers - invoke module processes
    // Each classifier returns: (path(report), path(processed), val(sample_id))
    kraken2_ch = params.kraken2
        ? Kraken2ClassificationPaired(
            analysis_id,
            reads_ch,
            params.kraken2_index,
            params.kraken2_params,
            params.classifier_process_script,
            params.python_bin,
            params.minimum_uniq_reads,
        )
        : null

    centrifuge_ch = params.centrifuge
        ? CentrifugeClassificationPaired(
            analysis_id,
            reads_ch,
            params.centrifuge_index,
            params.centrifuge_params,
            params.classifier_process_script,
            params.python_bin,
            params.minimum_uniq_reads,
        )
        : null

    diamond_ch = params.diamond
        ? DiamondClassificationPaired(
            analysis_id,
            reads_ch,
            params.diamond_index,
            params.diamond_params,
            params.classifier_process_script,
            params.python_bin,
            params.minimum_uniq_reads,
        )
        : null

    krakenunique_ch = params.krakenunique
        ? KrakenUniqueClassificationPaired(
            analysis_id,
            reads_ch,
            params.krakenunique_index,
            params.krakenunique_params,
            params.classifier_process_script,
            params.python_bin,
            params.minimum_uniq_reads,
        )
        : null



    // Merge classification results - pass report and processed files to module
    // Module expects: sample_id, kraken2_report, kraken2_processed, centrifuge_report, centrifuge_processed, krakenunique_report, krakenunique_processed, diamond_report, diamond_processed
    merge_classification_results_ch = MergeClassificationResults(
        analysis_id,
        params.kraken2 ? kraken2_ch.report : file('NO_FILE_kraken2_report.tsv'),
        params.kraken2 ? kraken2_ch._processed : file('NO_FILE_kraken2_processed.tsv'),
        params.centrifuge ? centrifuge_ch.report : file('NO_FILE_centrifuge_report.tsv'),
        params.centrifuge ? centrifuge_ch._processed : file('NO_FILE_centrifuge_processed.tsv'),
        params.krakenunique ? krakenunique_ch.report : file('NO_FILE_krakenunique_report.tsv'),
        params.krakenunique ? krakenunique_ch._processed : file('NO_FILE_krakenunique_processed.tsv'),
        params.diamond ? diamond_ch.report : file('NO_FILE_diamond_report.tsv'),
        params.diamond ? diamond_ch._processed : file('NO_FILE_diamond_processed.tsv'),
    )

    // Extract reference sequences
    reference_sequences_ch = ExtractReferenceSequences(analysis_id, merge_classification_results_ch)

    flattened_reference_sequences_ch = reference_sequences_ch.reference_sequences.flatMap { ref_list -> ref_list }

    combined_ch = reads_ch.combine(flattened_reference_sequences_ch)

    mapped_reads_ch = MapMinimap2Paired(analysis_id, combined_ch)

    filtered_alignments_ch = FilterBamMsamtools(mapped_reads_ch)
    sorted_reads_ch = sortBam(filtered_alignments_ch)
    coverage_ch = SamtoolsCoverage(sorted_reads_ch)

    coverage_ch = coverage_ch.map { file1, file2, _filename, _refname, group -> tuple(group, file1, file2) }
    coverage_ch = coverage_ch.groupTuple()

    mapping_files_info = mapped_reads_ch.map { file, group, _refname -> tuple(group, file) }
    mapping_files_info = mapping_files_info.groupTuple()

    merged_coverage_ch = MergeCoverageStatistics(coverage_ch)

    clustering_ch = ClusterMappedReads(analysis_id, mapping_files_info)

    MatchCladeReportWithReferenceSequences(
        analysis_id,
        clustering_ch.clade_report.ifEmpty([]),
        reference_sequences_ch.matched_assemblies,
        merged_coverage_ch.merged_coverage_statistics,
        merge_classification_results_ch,
    )
}


/*
* Add software column to classifier output to track source
*/
process AddSoftwareColumn {
    tag "AddSoftwareColumn ${analysis_id}_${software}"

    input:
    val analysis_id
    tuple val(query_id), path(report_file), path(processed_file)
    val software

    output:
    tuple val(query_id), path(report_file), path("${query_id}_${software}_with_software.tsv")

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    
    df = pd.read_csv("${processed_file}", sep="\t")
    df['software'] = '${software}'
    df.to_csv("${query_id}_${software}_with_software.tsv", sep="\t", index=False)
    """
}


/*
* Check reads_ch, if tuple contains only one file, append a placeholder NO_FILE for the second file
*/
process CheckReads {
    input:
    tuple path(fastq1), path(fastq2)

    output:
    tuple path(fastq1), path(fastq2)

    script:
    """
    if [ ! -f ${fastq2} ]; then
        echo "WARNING: Only one read file found, using NO_FILE placeholder"
        touch NO_FILE
    fi
    """
}


/*
* Quality control of the reads using prinseq++ (single-end)
*/
process QCReadsPrinseqSingle {
    tag "QCReadsPrinseqSingle ${analysis_id}"

    input:
    val analysis_id
    tuple path(fastq1)

    output:
    path "${analysis_id}_good_R1.fastq.gz", emit: _good_R1
    path "${analysis_id}_bad_R1.fastq.gz", emit: _bad_R1

    script:
    """
    prinseq++ ${params.prinseq_params} -fastq ${fastq1} -out_good ${analysis_id}_good_R1.fastq -out_bad ${analysis_id}_bad_R1.fastq
    bgzip ${analysis_id}_good_R1.fastq && bgzip ${analysis_id}_bad_R1.fastq
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
        #coverage_df.columns = ["chrom", "start", "end", "coverage"]
        coverage_df['file'] = file.split('/')[-1]
        coverage_df['error_rate'] = stats_df.loc[stats_df['stat'] == 'error rate:', 'value'].values[0]
        coverage_df['sample'] = stats_filename.replace("_coverage.txt", "").split("/")[-1]
        
        coverage_data.append(coverage_df)
    
    merged_df = pd.concat(coverage_data)
    merged_df.to_csv("merged_coverage_statistics.tsv", sep="\\t", index=False)
    """
}


/*
* Get coverage statistics using samtools
*/
process SamtoolsStats {
    tag "SamtoolsStats ${bamfile.baseName}"

    input:
    path bamfile

    output:
    path "${bamfile.baseName}_stats.txt"

    script:
    """
    samtools stats ${bamfile} > ${bamfile.baseName}_stats.txt
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
* Merge mapping statistics
*/
process MergeMappingStatistics {
    tag "MergeMappingStatistics ${query_id}"

    input:
    tuple val(query_id), path(stats_files)

    output:
    path "merged_mapping_statistics.tsv", emit: merged_mapping_statistics

    script:
    def stats_files_string = stats_files.collect { it[0] }.join(',')
    """
    #!/usr/bin/env python3
    import pandas as pd
    import glob
    
    stats_files_list = "${stats_files_string}".split(',')
    
    dfs = []
    for sf in stats_files_list:
        df = pd.read_csv(sf, sep="\\t")
        dfs.append(df)
    
    merged = pd.concat(dfs)
    merged.to_csv("merged_mapping_statistics.tsv", sep="\\t", index=False)
    """
}


/*
* Extract mapping statistics
*/
process ExtractMappingStatistics {
    tag "ExtractMappingStatistics ${bamfile.baseName}"

    input:
    path bamfile

    output:
    path "${bamfile.baseName}_mapping_stats.tsv"

    script:
    """
    #!/usr/bin/env python3
    import pysam
    import pandas as pd
    
    bam = pysam.AlignmentFile("${bamfile}", "rb")
    
    stats = []
    for read in bam:
        stats.append({
            'read_name': read.query_name,
            'reference': read.reference_name,
            'mapping_quality': read.mapping_quality,
            'flag': read.flag,
            'is_proper_pair': read.is_proper_pair,
        })
    
    df = pd.DataFrame(stats)
    df.to_csv("${bamfile.baseName}_mapping_stats.tsv", sep="\\t", index=False)
    bam.close()
    """
}


/*
* Cluster mapped reads
*/
process ClusterMappedReads {
    tag "ClusterMappedReads ${analysis_id}"
    publishDir { "${params.output_dir}/${analysis_id}" }, mode: 'copy'

    input:
    val analysis_id
    tuple val(query_id), path(mapped_reads)

    output:
    path "clustering/clade_report.tsv", emit: clade_report, optional: true
    path "clustering/sample_report.tsv", emit: sample_report, optional: true
    path "clustering/distance_matrix.tsv", emit: distance_matrix, optional: true
    path "clustering/all_node_statistics.tsv", emit: all_node_statistics, optional: true
    path "clustering/nj_tree_edges.txt", emit: nj_tree_edges, optional: true

    script:
    def mapped_reads_string = mapped_reads.collect { it[0] }.join(',')
    """
    ${params.map_to_matrix_bin} \
    -i . \
    -o clustering \
    ${params.map_to_matrix_params} 
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
* Extract reference sequences from classification results
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
    --mapping_references_dir "reference_sequences" \
    --include_term "complete" \
    --exclude_term "plasmid"
    """
}
