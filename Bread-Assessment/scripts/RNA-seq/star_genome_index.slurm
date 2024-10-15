#!/bin/bash
#SBATCH --job-name=genome_index          # Job name
#SBATCH --output=genome_index_%j.log     # Standard output log file
#SBATCH --error=genome_index_%j.err      # Standard error log file
#SBATCH --ntasks=1                       # Number of tasks (generally 1 for this kind of job)
#SBATCH --cpus-per-task=12               # Number of CPU cores per task
#SBATCH --mem=64G                        # Total memory per node (64 GB)
#SBATCH --time=72:00:00                  # Time limit hrs:min:sec
#SBATCH --partition=cpu               # Partition (or queue) name, can be different for your system

# Load STAR module or conda environment if necessary
module load Anaconda3 


source activate /home/ag003u1/.conda/envs/align-seq

# Running STAR to generate genome index
/home/ag003u1/.conda/envs/align-seq/bin/STAR --runThreadN 12 \
    --runMode genomeGenerate \
    --genomeDir /home/ag003u1/data/Reference \
    --genomeFastaFiles /home/ag003u1/data/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna \
    --sjdbGTFfile /home/ag003u1/data/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.gtf \
    --genomeSAsparseD 3 \
    --genomeChrBinNbits 14 \
    --limitGenomeGenerateRAM 64000000000 \
    --limitBAMsortRAM 64000000000
