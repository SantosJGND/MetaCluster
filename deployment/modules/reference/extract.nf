/*
 * Extract reference sequences from classifier output results to global reference database
 */
process ExtractReferenceSequences {
    tag "ExtractReferenceSequences ${sample_id}"

    input:
    val sample_id
    path classifier_output
    val python_bin
    val references_extract_script
    val assembly_store
    val include_term
    val exclude_term

    output:
    path "reference_sequences/*gz", emit: reference_sequences, optional: true
    path "reference_sequences/matched_assemblies.tsv", emit: matched_assemblies

    script:
    def include_arg = include_term ? "--include_term \"${include_term}\"" : ""
    def exclude_arg = exclude_term ? "--exclude_term \"${exclude_term}\"" : ""
    """
    ${python_bin} ${references_extract_script} retrieve \
        --input_table ${classifier_output} \
        --assembly_store "${assembly_store}" \
        --mapping_references_dir "reference_sequences" \
        ${include_arg} \
        ${exclude_arg}
    """
}
