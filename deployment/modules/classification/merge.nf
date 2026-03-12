/*
 * Merge classification results from different classifiers
 * Supports 2-classifier (basic) and 4-classifier (plus) modes
 */
process MergeClassificationResults {
    tag "MergeClassificationResults ${sample_id}"
    publishDir "${params.output_dir}/${sample_id}/classification", mode: 'copy'

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
    def has_krakenunique = krakenunique_report.toString() != 'NO_FILE' && !krakenunique_report.toString().endsWith('NO_FILE')
    def has_diamond = diamond_report.toString() != 'NO_FILE' && !diamond_report.toString().endsWith('NO_FILE')
    def kuniq_processed = has_krakenunique ? ", krakenunique_df" : ""
    def kuniq_cols = has_krakenunique ? "'kuniq_reads'" : ""
    def kuniq_merge = has_krakenunique ? "krakenunique_df = pd.read_csv(\"${krakenunique_processed}\", sep=\"\\t\").rename(columns={\"taxID\": \"taxid\", \"name\": \"description\"})\n    krakenunique_df['classifier'] = 'krakenunique'" : ""
    def diamond_processed_line = has_diamond ? ", diamond_df" : ""
    def diamond_cols = has_diamond ? "'diamond_reads'" : ""
    def diamond_merge = has_diamond ? "diamond_df = pd.read_csv(\"${diamond_processed}\", sep=\"\\t\").rename(columns={\"taxID\": \"taxid\", \"name\": \"description\"})\n    diamond_df['classifier'] = 'diamond'" : ""
    def suffixes = has_krakenunique && has_diamond ? '("_kraken2", "_centrifuge", "_krakenunique", "_diamond")' : has_krakenunique ? '("_kraken2", "_centrifuge", "_krakenunique")' : has_diamond ? '("_kraken2", "_centrifuge", "_diamond")' : '("_kraken2", "_centrifuge")'
    """
    #!/usr/bin/env python3
    import pandas as pd

    kraken2_df = pd.read_csv("${kraken2_processed}", sep="\\t").rename(columns={"taxID": "taxid", "name": "description"})
    kraken2_df['classifier'] = 'kraken2'

    centrifuge_df = pd.read_csv("${centrifuge_processed}", sep="\\t").rename(columns={"taxID": "taxid", "name": "description"})
    centrifuge_df['classifier'] = 'centrifuge'

    ${kuniq_merge}
    ${diamond_merge}

    # Build merge list
    dfs_to_merge = [kraken2_df, centrifuge_df${kuniq_processed}${diamond_processed_line}]

    # Start with kraken2 and centrifuge
    merged_df = pd.merge(kraken2_df, centrifuge_df, on=["description", "taxid"], how="outer", suffixes=("_kraken2", "_centrifuge"))

    # Add krakenunique if present
    ${has_krakenunique ? "merged_df = pd.merge(merged_df, krakenunique_df, on=[\"description\", \"taxid\"], how=\"outer\", suffixes=(\"_kraken2\", \"_krakenunique\"))" : ""}

    # Add diamond if present
    ${has_diamond ? "merged_df = pd.merge(merged_df, diamond_df, on=[\"description\", \"taxid\"], how=\"outer\", suffixes=(\"_kraken2\", \"_diamond\"))" : ""}

    merged_df = merged_df.drop_duplicates(subset=["taxid"])

    def classify(row):
        classifiers = []
        if pd.notna(row.get('classifier_kraken2')):
            classifiers.append('kraken2')
        if pd.notna(row.get('classifier_centrifuge')):
            classifiers.append('centrifuge')
        ${has_krakenunique ? "if pd.notna(row.get('classifier_krakenunique')):\n            classifiers.append('krakenunique')" : ""}
        ${has_diamond ? "if pd.notna(row.get('classifier_diamond')):\n            classifiers.append('diamond')" : ""}
        
        if not classifiers:
            return 'unclassified'
        return '/'.join(classifiers)

    # Calculate total unique reads
    kraken2_reads = merged_df.get('uniq_reads_kraken2', pd.Series([0]*len(merged_df))).fillna(0)
    cent_reads = merged_df.get('uniq_reads_centrifuge', pd.Series([0]*len(merged_df))).fillna(0)
    merged_df['total_uniq_reads'] = kraken2_reads + cent_reads
    ${has_krakenunique ? "kuniq_reads = merged_df.get('uniq_reads_krakenunique', pd.Series([0]*len(merged_df))).fillna(0)\n    merged_df['total_uniq_reads'] = merged_df['total_uniq_reads'] + kuniq_reads" : ""}
    ${has_diamond ? "diamond_reads = merged_df.get('uniq_reads_diamond', pd.Series([0]*len(merged_df))).fillna(0)\n    merged_df['total_uniq_reads'] = merged_df['total_uniq_reads'] + diamond_reads" : ""}

    merged_df = merged_df.sort_values(by='total_uniq_reads', ascending=False)
    merged_df['classification'] = merged_df.apply(classify, axis=1)

    merged_df.to_csv(f"${sample_id}_merged_classification.tsv", sep="\\t", index=False)
    """
}
