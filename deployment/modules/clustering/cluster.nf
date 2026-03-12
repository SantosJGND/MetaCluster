/*
 * Cluster mapped reads across alignment files
 */
process ClusterMappedReads {
    tag "ClusterMappedReads ${sample_id}"
    publishDir "${params.output_dir}/${sample_id}", mode: 'copy'

    input:
    val sample_id
    tuple val(query_id), path(mapped_reads)
    val map_to_matrix_bin
    val map_to_matrix_params

    output:
    path "clustering/clade_report.tsv", emit: clade_report, optional: true
    path "clustering/sample_report.tsv", emit: sample_report, optional: true
    path "clustering/distance_matrix.tsv", emit: distance_matrix, optional: true
    path "clustering/all_node_statistics.tsv", emit: all_node_statistics, optional: true
    path "clustering/nj_tree_edges.txt", emit: nj_tree_edges, optional: true

    script:
    def mapped_reads_string = mapped_reads.collect { it[0] }.join(',')
    """
    ${map_to_matrix_bin} \
        --files ${mapped_reads_string} \
        -o clustering \
        ${map_to_matrix_params}
    """
}
