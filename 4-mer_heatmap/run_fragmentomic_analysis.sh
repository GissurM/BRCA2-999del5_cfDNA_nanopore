#!/bin/bash

# Comprehensive cfDNA Fragmentomics Analysis Pipeline
# Processing all 24 T2T-aligned BAM files for CCCA motif analysis

set -e

SCRIPT_DIR="/location/of/this/script"
BAM_DIR="/path/to/your/bam/dir"
REFERENCE="/path/to/your/reference/fasta"
OUTPUT_DIR="$SCRIPT_DIR/fragmentomics_full"
RESULTS_DIR="$SCRIPT_DIR/results"
CHR_LIST="$SCRIPT_DIR/chromosomes.txt"

cd $SCRIPT_DIR

echo "=== cfDNA Fragmentomics Analysis Pipeline ==="
echo "Processing 24 hg38-aligned cfDNA samples"
echo "Focus: CCCA motif analysis for cancer detection"
echo "Reference: hg38"
echo ""

# Create output directories
mkdir -p $OUTPUT_DIR
mkdir -p $RESULTS_DIR

# Process each sample
# Process all BAM files matching the pattern
for BAM_FILE in "$BAM_DIR"/*_aligned.bam; do
    # Check if file exists (in case no files match the pattern)
    if [ ! -f "$BAM_FILE" ]; then
        continue
    fi

    # Extract sample name from filename
    SAMPLE=$(basename "$BAM_FILE" .bam)

    echo "=== Processing $SAMPLE ==="
    echo "BAM: $BAM_FILE"

    # Define all output file paths
    STATS_FILE="$OUTPUT_DIR/${SAMPLE}.stats"
    MOTIF_FILE="$OUTPUT_DIR/${SAMPLE}.motif"
    MOTIF_R_FILE="$RESULTS_DIR/${SAMPLE}.motif.R"
    MOTIF_CSV_FILE="$RESULTS_DIR/${SAMPLE}.motif.csv"

    # Check if final output already exists (fully processed)
    if [ -f "$MOTIF_CSV_FILE" ] && [ -f "$MOTIF_R_FILE" ]; then
        echo "✓ SKIPPING: $SAMPLE already fully processed"
        echo "  Found: ${SAMPLE}.motif.csv and ${SAMPLE}.motif.R"
        echo ""
        continue
    fi

    # Step 1: Extract fragments and motifs (skip if already done)
    if [ -f "$STATS_FILE" ] && [ -f "$MOTIF_FILE" ]; then
        echo "✓ Step 1: Using existing stats and motif files"
        echo "  Stats: ${SAMPLE}.stats ($(wc -l < "$STATS_FILE") lines)"
        echo "  Motif: ${SAMPLE}.motif ($(wc -l < "$MOTIF_FILE") lines)"
    else
        echo "Step 1: Extracting fragments and 5' motifs..."
        python nanopore_fragmentomics.py extract \
            --bam "$BAM_FILE" \
            --reference "$REFERENCE" \
            --output-dir "$OUTPUT_DIR" \
            --sample-name "$SAMPLE"
        
        if [ ! -f "$STATS_FILE" ] || [ ! -f "$MOTIF_FILE" ]; then
            echo "ERROR: Extraction failed for $SAMPLE"
            echo ""
            continue
        fi
        echo "✓ Extraction complete"
    fi

    # Step 2: Count motifs (only if not already done)
    if [ -f "$MOTIF_CSV_FILE" ]; then
        echo "✓ Step 2: Using existing motif counts"
        echo "  Found: ${SAMPLE}.motif.csv"
    else
        echo "Step 2: Counting motifs..."
        Rscript $SCRIPT_DIR/count_motif.R \
            "$CHR_LIST" \
            "$OUTPUT_DIR/" \
            "$OUTPUT_DIR/" \
            "$RESULTS_DIR/" \
            "$SAMPLE"
        

        if [ ! -f "$MOTIF_CSV_FILE" ]; then
            echo "ERROR: Motif counting failed for $SAMPLE"
            echo ""
            continue
        fi
        echo "✓ Motif counting complete"
    fi

    echo "✓ Completed: $SAMPLE"
    echo ""
done

echo "=== Analysis Complete ==="
echo "Results directory: $RESULTS_DIR"
echo "Summary files: *.motif.csv"

# Generate summary report
echo "Generating summary report..."
python -c "
import pandas as pd
import glob
import os

results_dir = '$RESULTS_DIR'
csv_files = glob.glob(os.path.join(results_dir, '*.motif.csv'))

if not csv_files:
    print('No motif CSV files found in results directory')
    exit(1)

summary_data = []
for file in sorted(csv_files):
    df = pd.read_csv(file)
    sample = os.path.basename(file).replace('.motif.csv', '')

    # Get CCCA frequency
    ccca_row = df[df['motif'] == 'CCCA']
    if not ccca_row.empty:
        ccca_freq = ccca_row['frequency'].iloc[0]
        ccca_count = ccca_row['count'].iloc[0]
    else:
        ccca_freq = 0
        ccca_count = 0

    total_frags = df['count'].sum()
    summary_data.append({
        'Sample': sample,
        'CCCA_Count': ccca_count,
        'CCCA_Frequency': ccca_freq,
        'Total_Fragments': total_frags
    })

summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(os.path.join(results_dir, 'CCCA_summary.csv'), index=False)
print('Summary saved to:', os.path.join(results_dir, 'CCCA_summary.csv'))
print(summary_df)
"

echo "Pipeline complete!"
