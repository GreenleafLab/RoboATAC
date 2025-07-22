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
  #library(biomaRt)
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

outdir <- "../../output/02-atac/03/"
plotdir <- "../../plots/02-atac/03/"
dir.create(plotdir, recursive=T, showWarnings=F)
dir.create(outdir, recursive=T, showWarnings=F)


dsa_norm <- loadDsAcc('../../output/02-atac/01/dsa_norm_rpkmlog2quantile/')

region <- "filteredConsensus"
  
## motif matching to JASPAR ------------------------------------------------
motifs <- getMotifs(database=JASPAR2020::JASPAR2020) 
saveRDS(motifs, paste0(outdir, "/JASPAR2020_Motifs.rds")

genome <- "BSgenome.Hsapiens.UCSC.hg38"
Peaks <- as.data.frame(dsa_norm@coord[[region]],row.names = NULL, stringsAsFactors=FALSE) %>% 
  dplyr::select(c("seqnames","start","end"))
Peaks.Ranges <- makeGRangesFromDataFrame(Peaks, ignore.strand = T)
motif.scores <- matchMotifs(motifs, Peaks.Ranges, genome = genome, out="scores")
saveRDS(motif.scores, paste0(outdir, "/motif.scores.JASPAR2020.rds"))

# matching with position information
motif.scores <- matchMotifs(motifs, Peaks.Ranges, genome = genome, out="positions")
saveRDS(motif.scores, paste0(outdir, "/motif.scores.JASPAR2020.pos.rds"))