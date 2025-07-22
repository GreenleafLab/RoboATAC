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
  library(patchwork)
  library(phastCons100way.UCSC.hg38)
})
fig <- function(width, heigth){
 options(repr.plot.width = width, repr.plot.height = heigth)
 }

source("../utils/tf_funcs.R")

deseq_dir <- "../../output/02-atac/04/deseq_data_Plate_TF_Dose/"
plotdir <- "../../plots/02-atac/11/"
outdir <- "../../output/02-atac/11/"

dir.create(plotdir, showWarnings = F, recursive=T)
dir.create(outdir, showWarnings = F, recursive=T)

padj_cutoff <- 0.05
log2fc_cutoff <- 0.58

# read motif matching results --------------------------------------------
motif.scores <- readRDS("../../output/02-atac/03/motif.scores.JASPAR2020.rds") # this is JASPAR motifs
#rownames(motif.scores)

# read in the motif peak names
dsa_norm <- loadDsAcc('../../output/02-atac/01/dsa_norm_rpkmlog2quantile/')
region <- "filteredConsensus"

Peaks <- as.data.frame(dsa_norm@coord[[region]],row.names = NULL, stringsAsFactors=FALSE) %>%
 dplyr::select(c("seqnames","start","end","name"))
rownames(motif.scores) <- Peaks$name
Peaks <- GRanges(Peaks)

# conservation score --------------------------------------------
phast <- phastCons100way.UCSC.hg38
allPeaksCons <- gscores(phast, Peaks, summaryFun=mean)

Peaks <- as.data.frame(allPeaksCons,row.names = NULL, stringsAsFactors=FALSE) %>% 
  dplyr::select(c("seqnames","start","end","name","default")) %>% dplyr::rename(conservation_score="default")
rownames(Peaks) <- Peaks$name
head(Peaks)

# correlate motif score and sensitivity params --------------------------------------------
motif_code_ls <- c(
                  #  "ALX4"="MA0681.2_PHOX2B",
                  #  "ELF1"="MA0473.3_ELF1",
                  #  "IRF4"="MA1419.1_IRF4",
                  #  "KLF1"="MA0039.4_KLF4",
                  #  "KLF4"="MA0039.4_KLF4",
                  #  "LEF1"="MA0768.1_LEF1",
                  #  "NR4A1"="MA1112.2_NR4A1",
                  #  "OCT4"="MA1115.1_POU5F1",
                  #  "PRDM1"="MA0508.3_PRDM1",
                  #  "SOX2"="MA0143.4_SOX2",
                  #  "SP4"="MA0685.1_SP4",
                  #  "SPI1"="MA0080.5_SPI1")
                   "TCF3"="MA0522.3_TCF3")


group1 = "GFP_1" # reference level
region = "filteredConsensus"

for (i in seq_along(motif_code_ls)){
  TF = names(motif_code_ls)[i]
  group2 = paste0(TF, "_1")
  motif_code = motif_code_ls[i]
  message(TF)
  
  
  # read sensitivity fitting data --------------------------------------------
  t <- read.csv(sprintf("../../output/02-atac/10/sensitivity_%s_log2fc0.58_fit.csv",TF), row.names = 1)
  # subset motif matching matrix to fitted peaks
  subset <- motif.scores[(rownames(motif.scores) %in% rownames(t)), motif_code]
  subset
  
  h_pval_cutoff_ls <- c(1, 0.2, 0.1)
  ka_pval_cutoff_ls <- c(1, 0.2, 0.1)
  for (k in seq_along(h_pval_cutoff_ls)){
    h_pval_cutoff <- h_pval_cutoff_ls[k]
    ka_pval_cutoff <- ka_pval_cutoff_ls[k]
    
    tsub <- t %>% dplyr::filter((res_linear - res_hill) > 0) # saturating sensitive
    #tsub <- t %>% dplyr::filter((res_linear - res_hill) < 0) # nonsaturating sensitive
    tsub <- tsub %>% dplyr::filter(h_pval<h_pval_cutoff & ka_pval<ka_pval_cutoff)
    scores <- motifScores(subset)[rownames(tsub),] %>% as.matrix %>% as.data.frame
    counts <- motifCounts(subset)[rownames(tsub),] %>% as.matrix %>% as.data.frame
    colnames(scores) <- paste0(colnames(subset), "_motif_score") 
    colnames(counts) <- paste0(colnames(subset), "_motif_counts")
    
    merged.tsub <- cbind(tsub,scores,counts)
    merged.tsub["conservation_score"] <- Peaks[rownames(merged.tsub), "conservation_score"]
    colnames(merged.tsub) <- gsub("\\/|\\|", "_",colnames(merged.tsub)) # clean up column names
    print(dim(merged.tsub))
    head(merged.tsub)
    
    currmotif_score <- paste0(motif_code, "_motif_score")
    currmotif_count <- paste0(motif_code, "_motif_counts")
    
    merged.tsub.subset <- merged.tsub %>% dplyr::filter(merged.tsub[[currmotif_count]]>0)
    merged.tsub.subset <- merged.tsub.subset %>%
      mutate(motif_score_bin = cut(merged.tsub.subset[[currmotif_score]], breaks = 8), 
             motif_count_bin = !!sym(currmotif_count))
    merged.tsub.subset[merged.tsub.subset[[currmotif_count]]>3, "motif_count_bin"] <- ">3"
    merged.tsub.subset$motif_score_bin <- factor(merged.tsub.subset$motif_score_bin, levels=merged.tsub.subset$motif_score_bin %>% unique %>% sort)
    merged.tsub.subset$motif_count_bin <- factor(merged.tsub.subset$motif_count_bin, levels=c("1", "2", "3", ">3"))
    
    dir.create(sprintf("%s/%s_subset", outdir, TF), recursive = T, showWarnings = F)
    saveRDS(merged.tsub.subset, sprintf("%s/%s_subset/%s_vs_ka_h_hpval%.2f_kapval%.2f.rds", outdir, TF, currmotif_count, h_pval_cutoff, ka_pval_cutoff))
    
    # conservation score vs motif score --------------------------------------------
    fig(8,4)
    p1 <- ggplot(merged.tsub.subset, aes(x=!!sym(currmotif_score), y=conservation_score)) + geom_point() +
      theme(text = element_text(size=20)) + ylab("conservation score") + xlab(currmotif_score) +
      geom_smooth(method='lm', formula= y~x) + stat_cor(aes(label =  ..rr.label..)) 

    p2 <- ggplot(merged.tsub.subset, aes(x=!!sym(currmotif_count), y=conservation_score)) + geom_point() +
      theme(text = element_text(size=20)) + ylab("conservation score") + xlab(currmotif_count) +
      geom_smooth(method='lm', formula= y~x) + stat_cor(aes(label =  ..rr.label..)) 

    dir.create(sprintf("%s/%s_subset", plotdir, TF), recursive = T, showWarnings = F)
    pdf(sprintf("%s/%s_subset/%s_%s_vs_conservation_hpval%.2f_kapval%.2f.pdf", plotdir, TF, currmotif_score, currmotif_count, h_pval_cutoff, ka_pval_cutoff), width=10, height=4)
    print(p1 + p2 + plot_layout(ncol=2))
    dev.off()

    # ka and h vs motif score ------------------------------------------------------
    fig(16,4)
    p0 <- ggplot(merged.tsub, aes(x=!!sym(currmotif_count))) + geom_bar() + 
      ggtitle(sprintf("# motifs of all fitted %s vs %s diff peaks (total %d peaks)", group2, group1, nrow(merged.tsub)))
    p1 <- ggplot(data=merged.tsub.subset, aes(y=ka,x=!!sym(currmotif_score)))+ geom_point() + 
      theme(text = element_text(size=20)) + ylab("EC50") + xlab(currmotif_score) +
      geom_smooth(method='lm', formula= y~x) + stat_cor(aes(label =  ..rr.label..)) 
    p2 <- ggplot(data=merged.tsub.subset, aes(y=ka,x=motif_score_bin))+ geom_boxplot() + 
      theme(text = element_text(size=20), axis.text.x=element_text(angle=45, hjust=1)) + 
      ylab("EC50") + xlab(paste0("binned ",currmotif_score))
    p3 <- ggplot(data=merged.tsub.subset, aes(y=h,x=!!sym(currmotif_score)))+ geom_point() + 
      theme(text = element_text(size=20)) + ylab("hill coef") + xlab(currmotif_score) +
      geom_smooth(method='lm', formula= y~x) + stat_cor(aes(label =  ..rr.label..))
    p4 <- ggplot(data=merged.tsub.subset, aes(y=h,x=motif_score_bin))+ geom_boxplot() + 
      theme(text = element_text(size=20), axis.text.x=element_text(angle=45, hjust=1))+ 
      ylab("hill coef") + xlab(paste0("binned ",currmotif_score))
    
    # split by motif counts
    plots <- list()
    for (count in levels(merged.tsub.subset$motif_count_bin)){
      df <- merged.tsub.subset %>% filter(motif_count_bin==count)
      pp1 <- ggplot(data=df, aes(y=ka,x=!!sym(currmotif_score)))+ geom_point() + 
        theme(text = element_text(size=20)) + ylab("EC50") + xlab(currmotif_score) +
        geom_smooth(method='lm', formula= y~x) + stat_cor(aes(label =  ..rr.label..)) +
        ggtitle(paste0("motif_count=",count))
      pp2 <- ggplot(data=df, aes(y=ka,x=motif_score_bin))+ geom_boxplot() + 
        theme(text = element_text(size=20), axis.text.x=element_text(angle=45, hjust=1)) + 
        ylab("EC50") + xlab(paste0("binned ",currmotif_score)) +
        ggtitle(paste0("motif_count=",count))
      pp3 <- ggplot(data=df, aes(y=h,x=!!sym(currmotif_score)))+ geom_point() + 
        theme(text = element_text(size=20)) + ylab("hill coef") + xlab(currmotif_score) +
        geom_smooth(method='lm', formula= y~x) + stat_cor(aes(label =  ..rr.label..)) +
        ggtitle(paste0("motif_count=",count))
      pp4 <- ggplot(data=df, aes(y=h,x=motif_score_bin))+ geom_boxplot() + 
        theme(text = element_text(size=20), axis.text.x=element_text(angle=45, hjust=1))+ 
        ylab("hill coef") + xlab(paste0("binned ",currmotif_score)) +
        ggtitle(paste0("motif_count=",count))
      p <- pp1 + pp2 + pp3 + pp4 + plot_layout(ncol=4)
      plots[[count]] <- p
    }
    dir.create(sprintf("%s/%s_subset", plotdir, TF), recursive = T, showWarnings = F)
    pdf(sprintf("%s/%s_subset/%s_vs_ka_h_hpval%.2f_kapval%.2f.pdf", plotdir, TF, currmotif_score, h_pval_cutoff, ka_pval_cutoff), width=18, height=4)
    print(p1 + p2 + p3 + p4 + plot_layout(ncol=4))
    print(plots)
    dev.off()
    
    # ka and h vs motif counts ------------------------------------------------------
    fig(16,4)
    p0 <- ggplot(merged.tsub, aes(x=!!sym(currmotif_count))) + geom_bar() + 
      ggtitle(sprintf("# motifs of all fitted %s vs %s diff peaks (total %d peaks)", group2, group1, nrow(merged.tsub)))
    p1 <- ggplot(data=merged.tsub.subset, aes(y=ka,x=!!sym(currmotif_count)))+ geom_point() + 
      theme(text = element_text(size=20)) + ylab("EC50") + xlab(currmotif_count) +
      geom_smooth(method='lm', formula= y~x) + stat_cor(aes(label =  ..rr.label..)) 
    p2 <- ggplot(data=merged.tsub.subset, aes(y=ka,x=motif_count_bin))+ geom_boxplot() + 
      theme(text = element_text(size=20), axis.text.x=element_text(angle=45, hjust=1)) + 
      ylab("EC50") + xlab(paste0("binned ",currmotif_count))
    p3 <- ggplot(data=merged.tsub.subset, aes(y=h,x=!!sym(currmotif_count)))+ geom_point() + 
      theme(text = element_text(size=20)) + ylab("hill coef") + xlab(currmotif_count) +
      geom_smooth(method='lm', formula= y~x) + stat_cor(aes(label =  ..rr.label..))
    p4 <- ggplot(data=merged.tsub.subset, aes(y=h,x=motif_count_bin))+ geom_boxplot() + 
      theme(text = element_text(size=20), axis.text.x=element_text(angle=45, hjust=1))+ 
      ylab("hill coef") + xlab(paste0("binned ",currmotif_count))
    
    # split by motif score
    plots <- list()
    for (bin in levels(merged.tsub.subset$motif_score_bin)){
      df <- merged.tsub.subset %>% filter(motif_score_bin==bin)
      pp1 <- ggplot(data=df, aes(y=ka,x=!!sym(currmotif_count)))+ geom_point() + 
        theme(text = element_text(size=20)) + ylab("EC50") + xlab(currmotif_count) +
        geom_smooth(method='lm', formula= y~x) + stat_cor(aes(label =  ..rr.label..)) +
        ggtitle(paste0("motif_score=", bin))
      pp2 <- ggplot(data=df, aes(y=ka,x=motif_count_bin))+ geom_boxplot() + 
        theme(text = element_text(size=20), axis.text.x=element_text(angle=45, hjust=1)) + 
        ylab("EC50") + xlab(paste0("binned ",currmotif_count)) +
        ggtitle(paste0("motif_score=", bin))
      pp3 <- ggplot(data=df, aes(y=h,x=!!sym(currmotif_count)))+ geom_point() + 
        theme(text = element_text(size=20)) + ylab("hill coef") + xlab(currmotif_count) +
        geom_smooth(method='lm', formula= y~x) + stat_cor(aes(label =  ..rr.label..)) +
        ggtitle(paste0("motif_score=", bin))
      pp4 <- ggplot(data=df, aes(y=h,x=motif_count_bin))+ geom_boxplot() + 
        theme(text = element_text(size=20), axis.text.x=element_text(angle=45, hjust=1))+ 
        ylab("hill coef") + xlab(paste0("binned ",currmotif_count)) +
        ggtitle(paste0("motif_score=", bin))
      p <- pp1 + pp2 + pp3 + pp4 + plot_layout(ncol=4)
      plots[[bin]] <- p
    }
    dir.create(sprintf("%s/%s_subset", plotdir, TF), recursive = T, showWarnings = F)
    pdf(sprintf("%s/%s_subset/%s_vs_ka_h_hpval%.2f_kapval%.2f.pdf", plotdir, TF, currmotif_count, h_pval_cutoff, ka_pval_cutoff), width=18, height=5)
    print(p0 + plot_layout(ncol=4))
    print(p1 + p2 + p3 + p4 + plot_layout(ncol=4))
    print(plots)
    dev.off()
  }
  
  
}



