#!/bin/bash

BAM_FILE="/media/raafat/Dr.Manal/S_littoralis/mapped_reads/combined_mapped.bam"  # Path to your combined BAM file
OUTPUT_DIR="trinity_output"
THREADS=16

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Run genome-guided Trinity assembly
Trinity --genome_guided_bam ${BAM_FILE} \
    --max_memory 50G \
    --CPU ${THREADS} \
    --output ${OUTPUT_DIR} \
    --genome_guided_max_intron 10000
