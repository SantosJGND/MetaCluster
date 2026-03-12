/*
 * Match clustering clade_report.txt and reference_sequences/matched_assemblies.tsv
 */
process MatchCladeReportWithReferenceSequences {
    tag "MatchCladeReportWithReferenceSequences ${sample_id}"
    publishDir "${params.output_dir}/${sample_id}/output", mode: 'copy'

    input:
    val sample_id
    path clade_report
    path matched_assemblies
    path coverage_report
    path merge_classification_results

    output:
    path "clade_report_with_references.tsv", emit: clade_report_with_references
    path clade_report
    path matched_assemblies
    path coverage_report

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    clade_report = pd.read_csv("${clade_report}", sep="\\t", header=None,  names=["clade", "nuniq", "freq", "min_pair_dist", "nfiles", "files"])
    clade_report['clade']
    clade_report['files'] = clade_report['files'].str.split(',')
    clade_report = clade_report.explode('files')

    matched_assemblies = pd.read_csv("${matched_assemblies}", sep="\\t")
    matched_assemblies['filename'] = matched_assemblies['assembly_file'].str.split('/').str[-1]

    coverage_report = pd.read_csv("${coverage_report}", sep="\\t")
    merged_classification_results = pd.read_csv("${merge_classification_results}", sep="\\t")

    def find_assembly_mapping(row):
        accession = row['assembly_accession']
        if accession is None or pd.isna(accession):
            row['clade'] = 'unmapped'
            row['nuniq'] = 0
            row['freq'] = 0
            return row
        match = clade_report[clade_report['files'].str.contains(accession, na=False)]
        if match.empty:
            row['clade'] = 'unmapped'
            row['nuniq'] = 0
            row['freq'] = 0
        else:
            row['clade'] = match['clade'].values[0]
            row['nuniq'] = match['nuniq'].values[0]
            row['freq'] = match['freq'].values[0]
        
        return row
    
    def find_assembly_coverage(row):
        accession = row['assembly_accession']
        if accession is None or pd.isna(accession):
            row['coverage'] = 0
            return row
        match = coverage_report[coverage_report['file'].str.contains(accession, na=False)]
        if match.empty:
            row['coverage'] = 0
        else:
            row['coverage'] = match['coverage'].values[0]
        
        return row
    
    def find_assembly_classification(row):
        taxid = row['taxid']
        match = merged_classification_results[merged_classification_results['taxid'] == taxid]
        if match.empty:
            row['classifier'] = 'unclassified'
        else:
            row['classifier'] = match['classification'].values[0]
        
        return row
    
    clade_report_with_references = matched_assemblies.apply(find_assembly_mapping, axis=1)
    clade_report_with_references = clade_report_with_references.apply(find_assembly_classification, axis=1)
    clade_report_with_references = clade_report_with_references.apply(find_assembly_coverage, axis=1)
    clade_report_with_references = clade_report_with_references[['description', 'taxid', 'assembly_accession', \
            'coverage', 'clade', 'nuniq', 'freq', 'classifier']]

    clade_report_with_references.to_csv("clade_report_with_references.tsv", sep="\\t", index=False)
    """
}
