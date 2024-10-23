#!/bin/bash

# Define directories and parameters
GENOME_DIR="/media/raafat/Dr.Manal/S_littoralis/Reference"
GTF_FILE="/media/raafat/Dr.Manal/S_littoralis/GCA_902850265.1_PGI_Spodlit_v1_genomic.gtf"
FASTQ_DIR="/media/raafat/Dr.Manal/S_littoralis"
OUTPUT_DIR="mapped_reads"
THREADS=16
MIN_OVERLAP=0.3
MISMATCH_MAX=2

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Loop through all paired-end fastq.gz files
for sample in $(ls ${FASTQ_DIR}/*120327_preprocessed.fastq.gz); do
    base=$(basename ${sample} 120327_preprocessed.fastq.gz)  # Get the base name of the sample
    R1="${FASTQ_DIR}/${base}120327_preprocessed.fastq.gz"    # Define the read 1 file


    # Run STAR for each pair of fastq files
    STAR --runThreadN ${THREADS} \
        --sjdbGTFfile ${GTF_FILE} \
        --readFilesCommand gunzip -c \
        --readFilesIn ${R1} \
        --genomeSAindexNbases 15 \
        --outSAMtype BAM SortedByCoordinate \
        --genomeDir ${GENOME_DIR} \
        --outFilterScoreMinOverLread ${MIN_OVERLAP} \
        --outFilterMatchNminOverLread ${MIN_OVERLAP} \
        --outFilterMismatchNmax ${MISMATCH_MAX} \
        --outFileNamePrefix ${OUTPUT_DIR}/mapped_${base}_

done

# Combine all BAM files into one using samtools
samtools merge -@ ${THREADS} ${OUTPUT_DIR}/combined_mapped.bam ${OUTPUT_DIR}/mapped_*_Aligned.sortedByCoord.out.bam

# Index the combined BAM file (optional)
samtools index ${OUTPUT_DIR}/combined_mapped.bam
