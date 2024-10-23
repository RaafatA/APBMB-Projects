#!/bin/bash

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n trinity-rna_seq star trinity samtools fastqc trimmomatic fastp 
