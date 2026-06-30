# sam_to_matrix

Cluster references based on read mapping similarity.

## Approach

1. Parse BAM/SAM files to get read-to-reference mappings
2. Build distance matrix from presence/absence or match counts
3. (Optional) Build NJ tree and identify private clades

## Usage

```bash
sam_to_matrix --input-directory /path/to/bam/files [options]
```

**Key options:**
- `--input-directory` - Directory containing BAM files
- `--max-reads N` - Max reads to process (default: 500000)
- `--frequency-threshold F` - Min frequency across files (default: 0.1)
- `--threshold T` - Private reads threshold for clades (default: 0.6)
- `--cluster-analysis/--no-cluster-analysis` - Enable/disable tree & clade analysis (default: true)
- `-o, --output-directory` - Output folder (default: "output")

## Output

**Always generated:**
- `distance_matrix.tsv`
- `presence_absence_matrix.tsv`

**With `--cluster-analysis` (default):**
- `nj_tree.newick`
- `nj_tree_edges.txt`
- `all_node_statistics.tsv`
- `clade_report.tsv`
- `sample_report.tsv`

## Acknowledgements

- [rust-phylogeny](https://github.com/RagnarGrootKoerkamp/rust-phylogeny) - Original NJ implementation