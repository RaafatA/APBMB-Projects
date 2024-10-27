#!/bin/bash

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n meta-shotgun snakemake fastqc trimmomatic multiqc bowtie2 samtools megahit quast prokka prodigal humann eggnog-mapper
conda activate meta-shotgun
