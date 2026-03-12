/*
 * Use samtools to extract BAM file coverage statistics
 */
process SamtoolsCoverage {
    tag "SamtoolsCoverage ${sample_id}_${reference_id}"

    input:
    tuple path(bamfile), path(bamindex), val(sample_id), val(reference_id)

    output:
    tuple path("${bamfile.baseName}.stats.txt"), path("${bamfile.baseName}.coverage.txt"), val(bamfile.baseName), val(reference_id), val(sample_id)

    script:
    """
    samtools coverage -o ${bamfile.baseName}.coverage.txt ${bamfile}
    samtools stats ${bamfile} | grep ^SN | cut -f 2- > ${bamfile.baseName}.stats.txt
    """
}

/*
 * Merge coverage statistics from different BAM files
 */
process MergeCoverageStatistics {
    tag "MergeCoverageStatistics ${sample_id}"

    input:
    tuple val(sample_id), path(stats_files), path(coverage_files)

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
        df = pd.read_csv(file, sep="\\t")
        df['file'] = file.split('/')[-1]
        df['error_rate'] = stats_df.loc[stats_df['stat'] == 'error rate:', 'value'].values[0]
        coverage_data.append(df)
    merged_df = pd.concat(coverage_data, ignore_index=True)
    

    merged_df.reset_index(inplace=True)
    merged_df.to_csv("merged_coverage_statistics.tsv", sep="\\t", index=False)
    """
}
