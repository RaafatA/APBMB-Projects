#!/bin/bash
#SBATCH --job-name=trinity_Denovo_assembly_guided          
#SBATCH --output=trinity_assembly_%j.log     
#SBATCH --error=trinity_assembly_%j.err      
#SBATCH --ntasks=1                           
#SBATCH --cpus-per-task=16                   
#SBATCH --mem=128G                           
#SBATCH --time=72:00:00                      
#SBATCH --partition=cpu                     

# Load Anaconda or necessary modules
module load Anaconda3  

# Activate conda environment (if Trinity is installed there)
source /eb/software/Anaconda3/2021.05/bin/activate
conda init
conda activate align-seq

# Install Trinity if not already installed in the environment
conda install -c bioconda trinity

# Define directories and parameters
BAM_FILE="mapped_reads/combined_mapped.bam"  # Path to your combined BAM file
OUTPUT_DIR="trinity_output"
THREADS=16

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Run genome-guided Trinity assembly
/home/ag003u1/.conda/envs/align-seq/bin/Trinity --genome_guided_bam ${BAM_FILE} \
    --max_memory 128G \
    --CPU ${THREADS} \
    --output ${OUTPUT_DIR} \
    --genome_guided_max_intron 10000
