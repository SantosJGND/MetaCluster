# Report: truth-taxid vs classifier-taxid mismatch (genus-rank truths)

Status: observation report only. No code, data, or configuration was changed.

## 1. Sample inspected

Path: `/home/bioinf/Desktop/INSA/Manuscripts/Clustering_Clinical_Metagenomics/data/data_samples/output_sample`

- 10 datasets: `dataset_0001_plan` ... `dataset_0010_plan` (170 files, about 1.8 GB).
- Per dataset: `fastq/` (R1/R2), `input/<ds>.tsv` (ground truth), `classification/` (per-classifier raw output and `<ds>_merged_classification.tsv`), `clustering/` (distance matrices), `output/` (`merged_coverage_statistics.tsv`, `matched_assemblies.tsv`, `clade_report_with_references.tsv`).
- Input rows record the simulated reference in `accid`, `assembly_accession`, and `path`.

## 2. Observed mismatches

### 2.1 Norovirus, truth taxid 142786

| Dataset | Simulated reads | Simulated accession | Centrifuge taxid reported | Centrifuge reads | In merged classification as 142786 |
|---|---|---|---|---|---|
| dataset_0004_plan | 14,561 | NC_045762.1 | 262897 (`norovirus_giv`) | 10,104 | no |
| dataset_0008_plan | 10 | NC_045762.1 | 262897 | 8 | no |
| dataset_0010_plan | 66,059 | NC_045762.1 | 262897 | 13,914 | no |

- NCBI resolves NC_045762.1 to taxid 262897 (*Norovirus GIV*). Taxid 142786 is the genus *Norovirus*.
- Classifiers detect the simulated organism at species taxid 262897. The genus taxid 142786 never appears in classifier output, so the truth taxon is scored `absent_classification`.
- Dataset_0010 also reports other norovirus taxa (552592 `norovirus_gii.17`, 1529918 `norovirus_gii/hu/jp/2007/...`), which are not in the truth table.

### 2.2 Rotavirus, truth taxid 10912

| Dataset | Simulated reads | Simulated accession |
|---|---|---|
| dataset_0004_plan | 92,639 | NC_135476.1 |
| dataset_0009_plan | 116,475 | NC_135476.1 |
| dataset_0010_plan | 71,640 | NC_135476.1 |

- NCBI resolves NC_135476.1 to taxid 2682575 (*Rotavirus K*). Taxid 10912 is the genus *Rotavirus*.
- Taxid 2682575 does not occur in any file in the sample.
- No rotavirus or sedoreovirus taxon appears in any classifier output for these datasets. The only related node is `Reovirales` (order), with 0 reads in dataset_0009.
- This is a classification miss, not only a taxid mismatch, at least in this sample. The classifier databases are not in the sample, so the cause (for example, reference absent from the index) is not established.

### 2.3 Wider pattern in the in-house reference tables

These come from the full test-set outputs at `data/results/viral/full_analysis/filter_direct_xgb_order/` in the manuscript package, not from this sample.

- `test_datasets_reference_leaves.tsv`: 15,304 of 18,561 leaf rows have `assembly_taxid != taxid`, across 41 input taxids.
- `test_datasets_input_df.tsv`: 857 of 4,603 input rows (18.6%) have `loss_stage = absent_classification`. All of them belong to the 41 mismatched taxids. Norovirus (111 rows) and rotavirus (116 rows) are two of the largest.
- Ancestry check using `taxa.db` (read-only), for leaf rows where the two taxids differ:
  - 193 rows: input taxid is an ancestor of the assembly taxid (e.g. 142786 is an ancestor of 262897).
  - 1,056 rows: assembly is not under the input taxid (e.g. HSV-2 10310, HHV-5 10359, P2 phage 10679, Echovirus 12078).
  - 12,076 rows: assembly taxid not in the partial `taxa.db`, so undetermined.
  - 1,979 rows: blank taxid.
  - Caveat: `taxa.db` is partial, so the 1,056 figure is a lower bound on non-ancestor cases only for the rows it could resolve.
- Internal inconsistency to check: 15 norovirus rows have `survived_post_cleanup = True` while `in_classification = False`.

### 2.4 Related evidence from the EsViritu evaluation

- The EsViritu identifier cache maps `nucleotide:NC_045762.1` to taxid 262897. EsViritu does detect this reference. The norovirus false negatives in the simulated EsViritu evaluation (396) are therefore a key mismatch, not absence.
- No cache entry exists for NC_135476.1, so EsViritu never reported that accession.

## 3. Likely cause (hypothesis, not verified)

Ground-truth taxids in the simulation input (`input/<ds>.tsv`) are genus-rank for these taxa (142786, 10912, and others), while classifiers, the merged classification, and the matched assemblies report species- or strain-rank taxids (262897, and assembly taxids for the simulated genomes). Matching is by exact taxid, so the genus truth is never equal to the reported taxid and is counted as absent.

## 4. Limitations

- The 1,998-dataset simulation root named for the EsViritu evaluation (`CODE/.../metagenomics-evaluation-pipeline`) does not contain them: `work/` is empty, and only two test fixtures have `input/*.tsv`.
- Classifier index contents (Centrifuge, Kraken2, Diamond, KrakenUnique) are not in the sample, so the rotavirus absence cannot be traced to the database.
- NCBI taxonomy results were taken from live queries (esummary, efetch) at the time of writing, not from a pinned taxonomy dump.

## 5. Not done (for the next step)

- No change to the evaluator, truth-matching rule, or manuscript.
- No re-scoring of the 10 sample datasets.
- No check of the classifier database contents.

## 6. Suggested next steps

1. Decide the truth-matching rule: lineage-ancestry (a predicted taxon matches a genus truth if its lineage contains the genus), or accession-derived species.
2. Locate the classifier database paths to test whether rotavirus is indexed.
3. Re-score the 10 sample datasets under the chosen rule and compare with the EsViritu truth for the same dataset names.
