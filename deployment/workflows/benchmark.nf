/*
 * Unified Benchmark Workflow
 * 
 * This workflow simulates metagenomic reads from taxonomic IDs,
 * classifies them, maps to references, and clusters the results.
 * 
 * Input: Taxonomic ID table (input_table)
 * Output: Clustered classification results with reference mapping
 * 
 * Usage:
 *   nextflow run workflows/benchmark.nf -params-file params.json
 * 
 * Parameters required:
 *   - input_table: Path to taxonomic ID table
 *   - output_dir: Output directory
 *   - assembly_store: Assembly database directory
 *   - technology: Sequencing technology (illumina, nanopore, etc.)
 *   - use_wgsim: Use wgsim instead of mess for simulation
 *   - use_plus: Include Diamond and KrakenUnique classifiers
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
include { WgsimSimulateReads } from '../modules/simulation/wgsim'
include { SimulateReadsMess } from '../modules/simulation/mess'
include { MatchCladeReportWithReferenceSequences } from '../modules/evaluation/match'

params.input_table = ""
params.output_dir = "output"
params.technology = "illumina"
params.use_wgsim = false
params.use_plus = false

workflow BENCHMARK {

    if (params.input_table == "") {
        error("Input table path is not provided. Please set the 'input_table' parameter.")
    }

    input_table_ch = Channel
        .fromPath(params.input_table)
        .ifEmpty { error("Cannot find the input table: ${params.input_table}") }

    // Extract reference sequences for the input table
    reference_sequences_ch = ExtractReferenceSequences(
        input_table_ch.baseName,
        input_table_ch,
        params.python_bin,
        params.references_extract_script,
        params.assembly_store,
        "complete",
        "plasmid"
    )

    // Simulate reads
    reads_ch = params.use_wgsim ? 
        WgsimSimulateReads(input_table_ch.baseName, input_table_ch, params.wgsim_python_path, params.wgsim_args, input_table_ch.baseName) :
        SimulateReadsMess(input_table_ch.baseName, input_table_ch, params.technology)

    reads_ch = reads_ch.map { sample_id, fastq_files ->
        def (fastq1, fastq2) = fastq_files
        tuple(sample_id, fastq1, fastq2)
    }

    // Quality control
    qc_ch = QCReadsPrinseqPaired(input_table_ch.baseName, reads_ch, params.prinseq_params)

    // Classify reads
    kraken2_ch = Kraken2ClassificationPaired(
        input_table_ch.baseName,
        qc_ch,
        params.kraken2_bin,
        params.kraken2_index,
        params.kraken2_params,
        params.classifier_process_script,
        params.python_bin,
        params.minimum_uniq_reads
    )

    centrifuge_ch = CentrifugeClassificationPaired(
        input_table_ch.baseName,
        qc_ch,
        params.centrifuge_index,
        params.centrifuge_params,
        params.classifier_process_script,
        params.python_bin,
        params.minimum_uniq_reads
    )

    // Merge classification results
    if (params.use_plus) {
        diamond_ch = DiamondClassificationPaired(
            input_table_ch.baseName,
            qc_ch,
            params.diamond_bin,
            params.diamond_index,
            params.diamond_params,
            params.classifier_process_script,
            params.python_bin,
            params.minimum_uniq_reads
        )

        krakenunique_ch = KrakenUniqueClassificationPaired(
            input_table_ch.baseName,
            qc_ch,
            params.krakenunique_bin,
            params.krakenunique_index,
            params.krakenunique_params,
            params.classifier_process_script,
            params.python_bin,
            params.minimum_uniq_reads
        )

        merge_ch = MergeClassificationResults(
            input_table_ch.baseName,
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
            input_table_ch.baseName,
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
        input_table_ch.baseName,
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
    combined_ch = qc_ch.combine(flattened_refs_ch)
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
        input_table_ch.baseName,
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
        input_table_ch.baseName,
        cluster_ch.clade_report,
        classified_refs_ch.matched_assemblies,
        merged_coverage_ch.merged_coverage_statistics,
        merge_ch
    )
}
