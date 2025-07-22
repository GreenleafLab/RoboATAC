#!/bin/bash

# set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate nucleoatac
ml biology samtools

ncores=8
genome_fasta=../../data/snakeatac_resources/genomes/hg38-no-haps/hg38-no-haps.fa

base_dir=../../
sample_meta=$base_dir/code/01-preprocessing/meta.txt
peak_dir=$base_dir/output/01-preprocessing/output/peaks/
bam_dir=$base_dir/output/01-preprocessing/output/bams/deduped/
out_dir=$base_dir/output/02-atac/12
mkdir -p $out_dir

#sample_ls=$(awk 'NR>1 {print $1}' $sample_meta)
sample_ls=$(awk 'NR>1 {print $1}' $sample_meta | grep '^HEK293T' | grep 'NR4A1')
#sample_ls=$(awk 'NR>1 {print $1}' $sample_meta | grep '^HEK293T_P3G7_SPI1_d100')
#sample_ls=$(awk 'NR>1 {print $1}' $sample_meta | grep 'HEK293T_P3H9')

for sample in $sample_ls; do
    echo "----------------------"
    echo $sample
    broadpeak="$peak_dir/${sample}_peaks.broadPeak"
    bam="$bam_dir/$sample.noMT.filtered.deduped.bam"
    prefix="$out_dir/$sample/$sample"
    split_dir="$out_dir/$sample/splits"
    mkdir -p $out_dir/$sample
    mkdir -p $split_dir

    # success=$prefix.nfrpos.bed.gz.tbi
    success=$prefix.occpeaks.bed.gz.tbi
    merge_flag=1
    if [ ! -f "$success" ]; then
        for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX; do
        #for chr in chr1; do 
        #for chr in chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX; do
            chr_dir="$split_dir/$chr"
            mkdir -p $chr_dir

            bam_chr=$chr_dir/$chr.bam
            peak_chr=$chr_dir/$chr.broadPeak
            prefix_chr=$chr_dir/$chr
            
            # rerun any incomplete chromosomes
            if [ ! -f $prefix_chr.occpeaks.bed.gz.tbi ]; then
                merge_flag=0
                cmd="samtools view -bo $bam_chr $bam $chr; samtools index $bam_chr; cat $broadpeak | grep $chr > $peak_chr; nucleoatac run --bed $peak_chr --bam $bam_chr --fasta $genome_fasta --out $prefix_chr --cores $ncores"
                echo $cmd
                if [ $chr = "chr1" ]; then # more time for chr1 which takes most time
                    sbatch -p wjg,sfgf,biochem --mem-per-cpu=4g --time=168:00:00 --job-name=nuc-$sample-$chr --out=slurm-logs/nucleoatac/slurm-%j-nucleoatac-$sample-$chr --cpus-per-task=$ncores --wrap="$cmd"
                else
                    sbatch -p wjg,sfgf,biochem --mem-per-cpu=4g --time=72:00:00 --job-name=nuc-$sample-$chr --out=slurm-logs/nucleoatac/slurm-%j-nucleoatac-$sample-$chr --cpus-per-task=$ncores --wrap="$cmd"
                fi
                
            else
                echo "skipping $chr"
            fi
        done

        # merge into a single file if all chromosomes finished
        if [ $merge_flag -eq 1 ]; then
            zcat $split_dir/*/*occpeaks.bed.gz > $prefix.occpeaks.bed
            bgzip $prefix.occpeaks.bed
            tabix -p bed $prefix.occpeaks.bed.gz
        fi
    else
        echo "skipping $sample"    
    fi    
done

