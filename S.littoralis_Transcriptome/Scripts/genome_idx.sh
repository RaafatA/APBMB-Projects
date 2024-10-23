#!/bin/bash

STAR --runThreadN 12 \
    --runMode genomeGenerate \
    --genomeDir /media/raafat/Dr.Manal/S_littoralis/Reference \
    --genomeFastaFiles /media/raafat/Dr.Manal/S_littoralis/GCA_902850265.1_PGI_Spodlit_v1_genomic.fna \
    --sjdbGTFfile /media/raafat/Dr.Manal/S_littoralis/GCA_902850265.1_PGI_Spodlit_v1_genomic.gtf
