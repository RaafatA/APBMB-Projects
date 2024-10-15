#!/bin/bash
#SBATCH --job-name=star_mapping               # Job name
#SBATCH --output=star_mapping_%j.log          # Output log
#SBATCH --error=star_mapping_%j.err           # Error log
#SBATCH --ntasks=1                            # Number of tasks (generally 1 for this kind of job)
#SBATCH --cpus-per-task=16                    # Number of CPU cores per task
#SBATCH --mem=64G                             # Total memory
#SBATCH --time=48:00:00                       # Time limit (hrs:min:sec)
#SBATCH --partition=cpu                       # Partition (use 'cpu')

# Load STAR and samtools modules or activate conda environment

module load Anaconda3  

# Alternatively, activate your conda environment if STAR and samtools are installed there
source /eb/software/Anaconda3/2021.05/bin/activate
conda init 
conda activate align-seq
conda install -c bioconda samtools

# Define directories and parameters
GENOME_DIR="/home/ag003u1/data/Reference"
GTF_FILE="/home/ag003u1/data/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.gtf"
FASTQ_DIR="/home/ag003u1/data/Transcriptome/"
OUTPUT_DIR="mapped_reads"
THREADS=16
MIN_OVERLAP=0.3
MISMATCH_MAX=2

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Loop through all paired-end fastq.gz files
for sample in $(ls ${FASTQ_DIR}/*_1.fq.gz); do
    base=$(basename ${sample} _1.fq.gz)  # Get the base name of the sample
    R1="${FASTQ_DIR}/${base}_1.fq.gz"    # Define the read 1 file
    R2="${FASTQ_DIR}/${base}_2.fq.gz"    # Define the read 2 file

    # Run STAR for each pair of fastq files
    /home/ag003u1/.conda/envs/align-seq/bin/STAR --runThreadN ${THREADS} \
        --sjdbGTFfile ${GTF_FILE} \
        --readFilesCommand gunzip -c \
        --readFilesIn ${R1} ${R2} \
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
