# getting read counts using kallisto
import os
import glob
import pandas as pd
import numpy as np

data_dir = "../../data/20TF_RNA_fastq/01.RawData"
output_dir = "../../output/03-rna/00/"
genome = "../../data/rna_resources/homo_sapiens_ensembl_v96/Homo_sapiens.GRCh38.cdna.all.synKLF1XBP1FOXO1.fa"
idx = "../../data/rna_resources/homo_sapiens_ensembl_v96/transcriptome_synKLF1XBP1FOXO1.idx" # kallisto index to use, automatically generate if does not exist
sample_sheet = "samples.csv"
config = pd.read_csv(sample_sheet, header=0)
sample_labels = list(config.sampleName)
os.chdir(output_dir)

localrule: all
rule all:
    input:
        expand("output/fastqs/qc/{sample_label}_R1_untrimmed_fastqc.html", sample_label = sample_labels),
        expand("output/fastqs/qc/{sample_label}_R1_trimmed_fastqc.html", sample_label = sample_labels),
        expand("output/kallisto/{sample_label}/abundance.tsv", sample_label = sample_labels),
    output:
        "rna_success.txt"
    shell:
        "echo $(date)  > {output};" 
        "echo created by Betty Liu and the Greenleaf lab >> {output};"

rule fastqc_unmapped_untrimmed:
    input:
        left = lambda w: config.loc[config.sampleName == w.sample_label]["R1Path"]
    output:
        lh = "output/fastqs/qc/{sample_label}_R1_untrimmed_fastqc.html",
        # stuff we don't really care about but want to eliminate when run is botched
        lz = "output/fastqs/qc/{sample_label}_R1_untrimmed_fastqc.zip",
    params:
        error_out_file = "error_files/{sample_label}_untrimmmed_fastqc",
        cores = "1",
        memory = "8000",
        job_name = "fastqc"
    benchmark:
        "benchmarks/fastqc/{sample_label}_untrimmed.txt"
    run:
        # files must be deleted if a previous run failed -- they do not have sample_label in them and cannot be listed as temp files
        shell("rm -f output/fastqs/qc/" + input[0].split('/')[-1].replace(".gz","").replace(".fastq","").replace(".fq","")+ "_fastqc.html"),
        shell("rm -f output/fastqs/qc/" + input[0].split('/')[-1].replace(".gz","").replace(".fastq","").replace(".fq","")+ "_fastqc.zip"),
        # run fastqc, which annoyingly doesnt allow you to format output files
        shell("fastqc {input.left} --outdir=output/fastqs/qc/"),
        # so rename output files manually
        shell("mv output/fastqs/qc/" + input[0].split('/')[-1].replace(".gz","").replace(".fastq","").replace(".fq","")  + "_fastqc.html {output.lh};"),
        shell("mv output/fastqs/qc/" + input[0].split('/')[-1].replace(".gz","").replace(".fastq","").replace(".fq","")  + "_fastqc.zip {output.lz};")

rule fastqc_unmapped_trimmed:
    input:
        left = "output/fastqs/trimmed/{sample_label}_R1_trimmed.fastq.gz",
    output:
        lh = "output/fastqs/qc/{sample_label}_R1_trimmed_fastqc.html",
        # stuff we don't really care about but want to eliminate when run is botched
        lz = "output/fastqs/qc/{sample_label}_R1_trimmed_fastqc.zip",
    params:
        error_out_file = "error_files/{sample_label}_trimmmed_fastqc",
        cores = "1",
        memory = "8000",
        job_name = "fastqc"
    benchmark:
        "benchmarks/fastqc/{sample_label}_trimmed.txt"
    run:
        # files must be deleted if a previous run failed -- they do not have sample_label in them and cannot be listed as temp files
        shell("rm -f output/fastqs/qc/" + input[0].split('/')[-1].replace(".gz","").replace(".fastq","").replace(".fq","")+ "_fastqc.html"),
        shell("rm -f output/fastqs/qc/" + input[0].split('/')[-1].replace(".gz","").replace(".fastq","").replace(".fq","")+ "_fastqc.zip"),
        # run fastqc
        shell("fastqc {input.left} --outdir=output/fastqs/qc/")

rule trim_adapters:
    input:
        left = lambda w: config.loc[config.sampleName == w.sample_label]["R1Path"],
    output:
        left = "output/fastqs/trimmed/{sample_label}_R1_trimmed.fastq.gz",
    params:
        error_out_file = "error_files/{sample_label}_trim",
        run_time = "2:30:00",
        cores = "1",
        memory = "8000",
        job_name = "trimming"
    benchmark:
        "benchmarks/trimming/{sample_label}.txt"
    shell:
        "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --minimum-length 20 --overlap=5 -o {output.left} {input.left} "


rule build_kallisto_index:
    input:
        genome
    output:
        idx
    benchmark:
        "benchmarks/build_kallisto_index/build_kallisto_index.txt"
    shell:
        "kallisto index -i {output} {input}"

rule kallisto:
    input:
        left = "output/fastqs/trimmed/{sample_label}_R1_trimmed.fastq.gz",
        idx = rules.build_kallisto_index.output
    output:
        "output/kallisto/{sample_label}/abundance.tsv"
    params:
        error_out_file = "error_files/{sample_label}_kallisto",
        cores = "2",
        memory = "8000",
        job_name = "kallisto",
        avg_read_length = 150,
        sd_read_length = 15
    benchmark:
        "benchmarks/kallisto/{sample_label}.txt"
    shell:
         "kallisto quant -i {input.idx} -o output/kallisto/{wildcards.sample_label} -b 20 -l {params.avg_read_length} -s {params.sd_read_length} --single <(zcat {input.left})"
