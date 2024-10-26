#!/bin/bash

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n shotgun fastp fastqc multiqc trimmomatic samtools bowtie megahit prokka humann prodigal 

source activate shotgun
