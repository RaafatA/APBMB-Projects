# Load config file
configfile: "config.yaml"

# Define sample names dynamically
SAMPLES = config["samples"].keys()

rule all:
    input:
        expand("qc_reports/{sample}_fastqc.html", sample=SAMPLES),
        "qc_reports/multiqc_report.html",
        "megahit_output/final.contigs.fa",
        "quast_report/report.txt",
        "prokka_output/prokka.gff",
        "humann_output/pathabundance.tsv",
        "eggNOG_annotations.emapper.annotations"

## 2. Quality Control
rule fastqc_initial:
    input:
        r1=lambda wildcards: config["samples"][wildcards.sample]["r1"],
        r2=lambda wildcards: config["samples"][wildcards.sample]["r2"]
    output:
        "qc_reports/{sample}_R1_fastqc.html",
        "qc_reports/{sample}_R2_fastqc.html"
    shell:
        "fastqc {input.r1} {input.r2} -o qc_reports/"

rule trimmomatic:
    input:
        r1=lambda wildcards: config["samples"][wildcards.sample]["r1"],
        r2=lambda wildcards: config["samples"][wildcards.sample]["r2"]
    output:
        trimmed_r1="{sample}_R1_trimmed.fastq",
        unpaired_r1="{sample}_R1_unpaired.fastq",
        trimmed_r2="{sample}_R2_trimmed.fastq",
        unpaired_r2="{sample}_R2_unpaired.fastq"
    params:
        adapters=config["adapters"]
    shell:
        "trimmomatic PE -threads {config[threads]} {input.r1} {input.r2} "
        "{output.trimmed_r1} {output.unpaired_r1} {output.trimmed_r2} {output.unpaired_r2} "
        "ILLUMINACLIP:{params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50"

rule fastqc_post_trimming:
    input:
        trimmed_r1="{sample}_R1_trimmed.fastq",
        trimmed_r2="{sample}_R2_trimmed.fastq"
    output:
        "qc_reports/{sample}_R1_trimmed_fastqc.html",
        "qc_reports/{sample}_R2_trimmed_fastqc.html"
    shell:
        "fastqc {input.trimmed_r1} {input.trimmed_r2} -o qc_reports/"

rule multiqc:
    input:
        expand("qc_reports/{sample}_fastqc.html", sample=SAMPLES)
    output:
        "qc_reports/multiqc_report.html"
    shell:
        "multiqc qc_reports/ -o qc_reports/"

## 4. Read Assembly (no host removal)
rule megahit_assembly:
    input:
        trimmed_r1="{sample}_R1_trimmed.fastq",
        trimmed_r2="{sample}_R2_trimmed.fastq"
    output:
        "megahit_output/final.contigs.fa"
    shell:
        "megahit -1 {input.trimmed_r1} -2 {input.trimmed_r2} -o megahit_output/ --min-contig-len 500"

rule quast_evaluation:
    input:
        contigs="megahit_output/final.contigs.fa"
    output:
        report="quast_report/report.txt"
    shell:
        "quast {input.contigs} -o quast_report/"

## 5. Gene Prediction
rule prokka_annotation:
    input:
        contigs="megahit_output/final.contigs.fa"
    output:
        "prokka_output/prokka.gff"
    shell:
        "prokka {input.contigs} --outdir prokka_output/ --metagenome --cpus {config[threads]}"

rule prodigal_prediction:
    input:
        contigs="megahit_output/final.contigs.fa"
    output:
        proteins="predicted_proteins.faa",
        genes="predicted_genes.fna",
        gff="gene_predictions.gff"
    shell:
        "prodigal -i {input.contigs} -a {output.proteins} -d {output.genes} -o {output.gff} -p meta"

## 6. Functional Annotation
rule humann_annotation:
    input:
        proteins="prokka_output/protein.faa"
    output:
        "humann_output/pathabundance.tsv"
    shell:
        "humann --input {input.proteins} --output humann_output/ --threads {config[threads]}"

rule eggnog_mapper:
    input:
        proteins="predicted_proteins.faa"
    output:
        "eggNOG_annotations.emapper.annotations"
    params:
        eggnog_db="/path/to/eggnog_db"
    shell:
        "emapper.py -i {input.proteins} -o eggNOG_annotations --cpu {config[threads]} --data_dir {params.eggnog_db} --usemem --itype proteins"
