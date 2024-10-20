## **2. Quality Control**
#### **2.1. Initial Quality Check (FastQC)**:
```bash
fastqc lake1_sample1_R1.fastq lake1_sample1_R2.fastq -o ./qc_reports/
```

#### **2.2. Adapter Trimming and Quality Filtering (Trimmomatic)**:
```bash 
trimmomatic PE -threads 8 \ lake1_sample1_R1.fastq lake1_sample1_R2.fastq \ lake1_sample1_R1_trimmed.fastq lake1_sample1_R1_unpaired.fastq \ lake1_sample1_R2_trimmed.fastq lake1_sample1_R2_unpaired.fastq \ ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \ LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
```

#### **2.3. Post-trimming Quality Check (FastQC & MultiQC)**:
```bash 
fastqc lake1_sample1_R1_trimmed.fastq lake1_sample1_R2_trimmed.fastq -o ./qc_reports/ multiqc ./qc_reports/ -o ./multiqc_report/
```

## **3. Host DNA Removal (Optional)**

```bash 
# 3.1. Build Index for Host Genome (Bowtie2)
bowtie2-build host_genome.fa host_index

# 3.2. Align Reads to Host Genome (Bowtie2)
bowtie2 -x host_index -1 lake1_sample1_R1_trimmed.fastq -2 lake1_sample1_R2_trimmed.fastq -S lake1_sample1_host_aligned.sam --threads 8

# 3.3. Remove Host Reads (Extract Unaligned Reads with SAMtools):
samtools view -b -f 12 -F 256 lake1_sample1_host_aligned.sam | samtools bam2fq - > lake1_sample1_host_removed.fastq
```

## 4. Read Assembly
``` bash 
# 4.1. Assemble Reads into Contigs (MEGAHIT):
megahit -1 lake1_sample1_host_removed.fastq -2 lake1_sample2_host_removed.fastq -o ./megahit_output/ --min-contig-len 500

# Su
megahit -1 lake1_sample1_R1_trimmed.fastq -2 lake1_sample1_R2_trimmed.fastq -o ./megahit_output/ --k-min 27 --k-step 10 --min-contig-len 500


# Assembly Quality Evaluation 
'''Use QUAST to evaluate the quality of the assembly, including metrics like N50, total length, and number of contigs.'''
quast megahit_output/final.contigs.fa -o quast_report/

```

## 5. Gene Prediction (on Assembled Contigs)
``` bash 
# 5.1. Predict Genes (PROKKA):
prokka ./megahit_output/final.contigs.fa --outdir ./prokka_output/ --metagenome --cpus 8

# Prodigal v2.6.3
'''Use **Prodigal v2.6.3** to predict protein-coding sequences (CDS) from the assembled contigs. Run in metagenomic mode (`-p meta`) to account for fragmented contigs.'''
prodigal -i megahit_output/final.contigs.fa -a predicted_proteins.faa -d predicted_genes.fna -o gene_predictions.gff -p meta

```

## **6. Functional Annotation**
```bash 
# 6.1. Functional Annotation of Predicted Genes (HUMAnN3 on assembled contigs):
'''If you’re using assembled contigs (from `PROKKA`), map the predicted gene products (protein FASTA) to functional databases:'''
humann --input prokka_output/protein.faa --output humann_output/ --threads 8
```
#### Functional Annotation with eggNOG-mapper
Assign functions to the predicted protein-coding genes using **eggNOG-mapper v2.1.8** with the **DIAMOND** search mode. This step maps the predicted genes to orthologous groups in the **eggNOG v5.0.2** database.
```bash 
emapper.py -i predicted_proteins.faa -o eggNOG_annotations --cpu 8 --data_dir /path/to/eggnog_db --usemem --itype proteins

```
