# args ---------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

dsa_path <- args[1]
region <- args[2]
motifDb <- args[3]
motif_score_file <- args[4]
dds_file <- args[5]
group1 <- args[6] # reference level
group2 <- args[7] # sample level
contrast.col <- args[8]
p.cutoff <- as.numeric(args[9])
log2fc.cutoff <- as.numeric(args[10])
deseq_dir <- args[11]
plot_dir <- args[12]
motif_dir <- args[13]

# main ---------------------------------------------------------------------------------

# script modified from Sidney V.
# Import Libraries
suppressWarnings({
  library(ChrAccR)
  library(muLogR)
  library(muRtools)
  library(ggplot2)
  library(chromVAR)
  library(motifmatchr)
  library(JASPAR2018)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(dplyr)
  library(DESeq2)
  library(GenomicFeatures)
  library(ChrAccRAnnotationHg38)
  library(tidyverse)
  library(ggpubr)
  library(ggrastr)
  library(ggrepel)
  library(viridis)
  library(Gviz)
  library(TFBSTools)
  library(reshape2)
})

source("../utils/tf_funcs.R")

## read motif matching to JASPAR ------------------------------------------------
dsa_norm <- loadDsAcc(dsa_path)

motifs <- readRDS(motifDb)
Peaks <- as.data.frame(dsa_norm@coord[[region]],row.names = NULL, stringsAsFactors=FALSE) %>% 
  dplyr::select(c("seqnames","start","end"))
motif.scores <- readRDS(motif_score_file)

# TF Motif Enrichment on differential peaks by {tf} -------------------------------
dds <- readRDS(dds_file) 

message(paste0(contrast.col, ": ", group2, " vs ", group1))
tf_motif_de_analysis(dds, dsa_norm, motifs, motif.scores, Peaks, 
                      group1, group2, contrast.col, region, p.cutoff=p.cutoff, log2fc.cutoff=log2fc.cutoff,
                      deseq_dir=deseq_dir, plot_dir=plot_dir, motif_dir=motif_dir) # function in tf_funcs.R
