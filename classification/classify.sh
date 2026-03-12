#!/bin/bash

CENTRIGUE_BIN="/path/to/centrigue/bin"
CENTRIFUGE_XPATH="/path/to/centrifuge/db"
CENTRIGUE_ARGS="--report-file centrifuge_report.txt --threads 4"

KRAKEN2_BIN="/path/to/kraken2/bin"
KRAKEN2_XPATH="/path/to/kraken2/db"

INPUT_SAMPLE_NAME="sample1"
INPUT_R1="sample1_R1.fastq.gz"
INPUT_R2="sample1_R2.fastq.gz"

OUTPUT_DIR="output/classification_results/"
CENTRIFUGE_OUTPUT_DIR="output/centrifuge_results/"
KRAKEN2_OUTPUT_DIR="output/kraken2_results/"

mkdir -p $CENTRIFUGE_OUTPUT_DIR
mkdir -p $KRAKEN2_OUTPUT_DIR

# Run Centrifuge classification
$CENTRIGUE_BIN"centrifuge" \
    -x $CENTRIFUGE_XPATH \
    -1 $INPUT_R1 \
    -2 $INPUT_R2 \
    -S $CENTRIFUGE_OUTPUT_DIR$INPUT_SAMPLE_NAME".tsv" \
    --output $CENTRIFUGE_OUTPUT_DIR$INPUT_SAMPLE_NAME".txt" \
    --report-file $CENTRIFUGE_OUTPUT_DIR$INPUT_SAMPLE_NAME"_report.txt" \
    $CENTRIGUE_ARGS

# Run Kraken2 classification
$KRAKEN2_BIN"kraken2" \
    --db $KRAKEN2_XPATH \
    --threads 4 \
    --report $KRAKEN2_OUTPUT_DIR$INPUT_SAMPLE_NAME"_report.txt" \
    "--gzip-compressed" \
    --output $KRAKEN2_OUTPUT_DIR$INPUT_SAMPLE_NAME".txt" \
    $INPUT_R1 $INPUT_R2 