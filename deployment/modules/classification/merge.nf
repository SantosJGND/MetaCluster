/*
 * Merge classification results from different classifiers
 * Supports 2-classifier (basic) and 4-classifier (plus) modes
 */
process MergeClassificationResults {
    tag "MergeClassificationResults ${sample_id}"
    publishDir { "${params.output_dir}/${sample_id}/classification" }, mode: 'copy'

    input:
    val sample_id
    path kraken2_report
    path kraken2_processed
    path centrifuge_report
    path centrifuge_processed
    path krakenunique_report
    path krakenunique_processed
    path diamond_report
    path diamond_processed

    output:
    path "${sample_id}_merged_classification.tsv"

    script:
    def has_krakenunique = krakenunique_report.toString() != 'NO_FILE' && !krakenunique_report.toString().startsWith('NO_FILE')
    def has_diamond = diamond_report.toString() != 'NO_FILE' && !diamond_report.toString().startsWith('NO_FILE')
    def kuniq_processed = has_krakenunique ? "krakenunique_df," : ""
    def diamond_processed_line = has_diamond ? "diamond_df," : ""
    """
    #!/usr/bin/env python3
    import pandas as pd

    kraken2_df = pd.read_csv("${kraken2_processed}", sep="\t").rename(columns={"taxID": "taxid", "name": "description"})

    centrifuge_df = pd.read_csv("${centrifuge_processed}", sep="\t").rename(columns={"taxID": "taxid", "name": "description"})

    ${has_krakenunique ? "krakenunique_df = pd.read_csv(\"${krakenunique_processed}\", sep=\"\\t\").rename(columns={\"taxID\": \"taxid\", \"name\": \"description\"})" : ""}
    ${has_diamond ? "diamond_df = pd.read_csv(\"${diamond_processed}\", sep=\"\\t\").rename(columns={\"taxID\": \"taxid\", \"name\": \"description\"})" : ""}

    all_dfs = [kraken2_df, centrifuge_df${has_krakenunique ? ", " + kuniq_processed : ""}${has_diamond ? " " + diamond_processed_line : ""}]
    combined_df = pd.concat(all_dfs, ignore_index=True)

    merged_df = combined_df.groupby('taxid', as_index=False).agg({
        'description': 'first',
        'uniq_reads': 'sum',
        'software_name': lambda x: '/'.join(sorted(set(x)))
    })

    merged_df = merged_df.rename(columns={'software_name': 'classifiers'})
    merged_df = merged_df.sort_values(by='uniq_reads', ascending=False)

    merged_df.to_csv(f"${sample_id}_merged_classification.tsv", sep="\t", index=False)
    """
}
