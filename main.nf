/*
 * Mapping Cluster: Metagenomics Cross-Hit Resolution Benchmark
 *
 * Entry point for the metagenomics evaluation pipeline.
 * Supports the full benchmark workflow and a classification-only workflow.
 *
 * Usage:
 *   nextflow run main.nf -profile conda -params-file deployment/params.json
 *   nextflow run main.nf -profile conda -params-file deployment/params.json --workflow classify_only
 */

include { BENCHMARK } from './deployment/workflows/benchmark'
include { CLASSIFY_ONLY } from './deployment/workflows/classify_only'

params.workflow = "benchmark"

workflow {
    switch (params.workflow) {
        case "classify_only":
            CLASSIFY_ONLY()
            break
        case "benchmark":
        default:
            BENCHMARK()
            break
    }
}
