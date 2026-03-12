/*
 * Unified Classification-Only Workflow
 * 
 * This workflow classifies existing FASTQ reads, maps to references,
 * and clusters the results. No simulation step.
 * 
 * Input: FASTQ reads directory (reads)
 * Output: Clustered classification results with reference mapping
 * 
 * Usage:
 *   nextflow run workflows/classify_only.nf -params-file params.json
 * 
 * Parameters required:
 *   - reads: Path to directory containing FASTQ files
 *   - output_dir: Output directory
 *   - assembly_store: Assembly database directory
 *   - use_plus: Include Diamond and KrakenUnique classifiers
 *   - analysis_id: Optional analysis ID (auto-generated if not provided)
 */

include { QCReadsPrinseqPaired } from '../modules/qc/prinseq'
include { Kraken2ClassificationPaired } from '../modules/classification/kraken2'
include { CentrifugeClassificationPaired } from '../modules/classification/centrifuge'
include { DiamondClassificationPaired } from '../modules/classification/diamond'
include { KrakenUniqueClassificationPaired } from '../modules/classification/krakenunique'
include { MergeClassificationResults } from '../modules/classification/merge'
include { MapMinimap2Paired } from '../modules/mapping/minimap2'
include { FilterBamMsamtools } from '../modules/mapping/filter'
include { SortBam } from '../modules/mapping/sort'
include { SamtoolsCoverage; MergeCoverageStatistics } from '../modules/mapping/coverage'
include { ClusterMappedReads } from '../modules/clustering/cluster'
include { ExtractReferenceSequences } from '../modules/reference/extract'
include { MatchCladeReportWithReferenceSequences } from '../modules/evaluation/match'

params.reads = ""
params.output_dir = "output"
params.analysis_id = ""
params.use_plus = false

workflow CLASSIFY_ONLY {

    analysis_id = params.analysis_id ?: UUID.randomUUID().toString()

    if (params.reads == "") {
        error("Reads directory path is not provided. Please set the 'reads' parameter.")
    }

    // Read FASTQ files
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

    // Quality control
    qc_ch = QCReadsPrinseqPaired(analysis_id, reads_ch, params.prinseq_params)
    
    r1_ch = qc_ch.out

    // Classify reads
    kraken2_ch = Kraken2ClassificationPaired(
        analysis_id,
        r1_ch,
        params.kraken2_bin,
        params.kraken2_index,
        params.kraken2_params,
        params.classifier_process_script,
        params.python_bin,
        params.minimum_uniq_reads
    )

    centrifuge_ch = CentrifugeClassificationPaired(
        analysis_id,
        r1_ch,
        params.centrifuge_index,
        params.centrifuge_params,
        params.classifier_process_script,
        params.python_bin,
        params.minimum_uniq_reads
    )

    // Merge classification results
    if (params.use_plus) {
        diamond_ch = DiamondClassificationPaired(
            analysis_id,
            r1_ch,
            params.diamond_bin,
            params.diamond_index,
            params.diamond_params,
            params.classifier_process_script,
            params.python_bin,
            params.minimum_uniq_reads
        )

        krakenunique_ch = KrakenUniqueClassificationPaired(
            analysis_id,
            r1_ch,
            params.krakenunique_bin,
            params.krakenunique_index,
            params.krakenunique_params,
            params.classifier_process_script,
            params.python_bin,
            params.minimum_uniq_reads
        )

        merge_ch = MergeClassificationResults(
            analysis_id,
            kraken2_ch[1],
            kraken2_ch[0],
            centrifuge_ch[1],
            centrifuge_ch[0],
            krakenunique_ch[1],
            krakenunique_ch[0],
            diamond_ch[1],
            diamond_ch[0]
        )
    } else {
        merge_ch = MergeClassificationResults(
            analysis_id,
            kraken2_ch[1],
            kraken2_ch[0],
            centrifuge_ch[1],
            centrifuge_ch[0],
            file('NO_FILE'),
            file('NO_FILE'),
            file('NO_FILE'),
            file('NO_FILE')
        )
    }

    // Extract reference sequences from classification
    classified_refs_ch = ExtractReferenceSequences(
        analysis_id,
        merge_ch,
        params.python_bin,
        params.references_extract_script,
        params.assembly_store,
        "complete",
        "plasmid"
    )

    flattened_refs_ch = classified_refs_ch.reference_sequences
        .flatMap { ref_list -> ref_list }

    // Map reads to references
    combined_ch = r1_ch.combine(flattened_refs_ch)
    mapped_ch = MapMinimap2Paired(combined_ch, params.minimap2_illumina_params)

    // Process BAM files
    filtered_ch = FilterBamMsamtools(mapped_ch, params.msamtools_params)
    sorted_ch = SortBam(filtered_ch)
    coverage_ch = SamtoolsCoverage(sorted_ch)

    // Group coverage files by sample
    coverage_grouped_ch = coverage_ch
        .map { stats, cov, filename, refname, sample -> tuple(sample, stats, cov) }
        .groupTuple()

    merged_coverage_ch = MergeCoverageStatistics(coverage_grouped_ch)

    // Group mapping files for clustering
    mapping_grouped_ch = mapped_ch
        .map { bam, sample, refname -> tuple(sample, bam) }
        .groupTuple()

    // Cluster mapped reads
    cluster_ch = ClusterMappedReads(
        analysis_id,
        mapping_grouped_ch,
        params.map_to_matrix_bin,
        params.map_to_matrix_params
    )

    cluster_ch.clade_report.ifEmpty {
        log.info("No clustering results generated. Ending workflow.")
        System.exit(0)
    }

    // Match results with references
    MatchCladeReportWithReferenceSequences(
        analysis_id,
        cluster_ch.clade_report,
        classified_refs_ch.matched_assemblies,
        merged_coverage_ch.merged_coverage_statistics,
        merge_ch
    )
}
