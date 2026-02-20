#!/bin/bash

# Simple bash script that replicates 06_blast.nf functionality
# This script downloads UniProt data, creates a BLAST database, and runs BLASTP

module load blast-plus

set -e  # Exit on error
set -u  # Exit on undefined variable

# Parameters (equivalent to params in Nextflow)
QUERY_FASTA="${1:-data/query.fasta}"
DB_QUERY="reviewed:true%20AND%20organism_id:559292"
OUTDIR="${2:-results}"
WORK_DIR="work_bash"

# Create working and output directories
mkdir -p "$WORK_DIR"
mkdir -p "$OUTDIR"

echo "=== Step 1: Downloading UniProt yeast proteins ==="
curl -L \
  "https://rest.uniprot.org/uniprotkb/stream?query=${DB_QUERY}&format=fasta&compressed=false" \
  -o "${WORK_DIR}/uniprot_yeast.fasta"

echo "=== Step 2: Creating BLAST database ==="
makeblastdb \
  -in "${WORK_DIR}/uniprot_yeast.fasta" \
  -dbtype prot \
  -out "${WORK_DIR}/uniprot_yeast"

echo "=== Step 3: Running BLASTP search ==="
blastp \
  -query "$QUERY_FASTA" \
  -db "${WORK_DIR}/uniprot_yeast" \
  -outfmt "6 qseqid sseqid pident length evalue bitscore" \
  > "${WORK_DIR}/blast.tsv"

# Copy results to output directory
cp "${WORK_DIR}/blast.tsv" "${OUTDIR}/blast.tsv"

echo "=== Pipeline complete! ==="
echo "Results saved to: ${OUTDIR}/blast.tsv"
echo "Working directory: ${WORK_DIR}"
