library(tidyverse)
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(VennDiagram)
library(data.table)

source("../utils/tf_funcs.R")

# peak type composition of all consensus peaks -------------------------------------------------------
outdir_main <- "../../output/02-atac/08/"
plotdir_main <- "../../plots/02-atac/08/"
dir.create(plotdir_main, recursive=T, showWarnings=F)
dir.create(outdir_main, recursive=T, showWarnings=F)

peaks <- read.table('../../output/02-atac/01/consensus_peaks_HEK293T.bed')
colnames(peaks) <- c("seqnames", "start", "end", "name", "score", "strand")
peaks <- GRanges(peaks)

addArchRGenome("hg38")
geneAnnot <- getArchRGenome(geneAnnotation = T, genomeAnnotation = F) %>% as.list
peaks_anno <- ArchR:::.fastAnnoPeaks(peaks, BSgenome = BSgenome.Hsapiens.UCSC.hg38, geneAnnotation = geneAnnot)
saveRDS(peaks_anno, paste0(outdir_main, "/peaks_anno_all_consensus.rds"))

p1 <- ggplot(peaks_anno %>% as.data.frame, aes(x=1, fill=peakType)) + geom_bar()
p2 <- ggplot(peaks_anno %>% as.data.frame, aes(x=1, fill=peakType)) + geom_bar(position="fill") + ylab("count%")
pdf(paste0(plotdir_main, "/peak_type_all_consensus.pdf"), width=3, height=6)
p1
p2
dev.off()

# peak type composition of differential peaks, TF group ---------------------------------------------
plotdir <- paste0(plotdir_main, "/tf_group")
outdir <- paste0(outdir_main, "/tf_group")
dir.create(plotdir, recursive=T, showWarnings=F)
dir.create(outdir, recursive=T, showWarnings=F)

atacdiff_dir <- "../../output/02-atac/04/deseq_data_TF/"
meta <- read_tsv("atac_sample_sheet.txt")

celltype <- "HEK293T"
cond2 <- "GFP"
log2fc.cutoff <- 0.58
padj.cutoff <- 0.05

submeta <- meta[meta$CellType==celltype,]
cond1_ls <- unique(submeta$TF)
cond1_ls <- cond1_ls[cond1_ls != cond2]

peaktype_all <- list()
for (cond1 in cond1_ls){
  message(cond1)
  atacdiff <- fread(paste0(atacdiff_dir, "/deseq_", cond1, "_vs_", cond2, "_filteredConsensus_full.tsv"))
  atacdiff_sig <- atacdiff %>% dplyr::filter(padj<=padj.cutoff & abs(log2FoldChange)>=log2fc.cutoff)
  peaks <- GRanges(atacdiff_sig)
  peaks_anno <- ArchR:::.fastAnnoPeaks(peaks, BSgenome = BSgenome.Hsapiens.UCSC.hg38, geneAnnotation = geneAnnot)
  saveRDS(peaks_anno, paste0(outdir, "/peaks_anno_diffpeak_", cond1, "_vs_", cond2, ".rds"))
  
  out <- peaks_anno %>% as.data.frame %>% dplyr::select(peakType) %>% table %>% as.data.frame %>% dplyr::rename(peakType=".", count="Freq") %>% dplyr::mutate(TF=cond1)
  peaktype_all[[cond1]] <- out
}
peaktype_all <- do.call(rbind, peaktype_all)
saveRDS(peaktype_all, paste0(outdir, "/peak_type_all.rds"))

p1 <- ggplot(peaktype_all, aes(x=TF, y=count, fill=peakType)) + geom_col(position="stack") + theme(axis.text.x=element_text(angle=90, hjust=1))
p2 <- ggplot(peaktype_all, aes(x=TF, y=count, fill=peakType)) + geom_col(position="fill") + ylab("count %") + theme(axis.text.x=element_text(angle=90, hjust=1))
pdf(paste0(plotdir, "/peak_type_diffpeak_all.pdf"), width=10, height=6)
p1
p2
dev.off()

# peak type composition of differential peaks, TF Dose ---------------------------------------------
plotdir <- paste0(plotdir_main, "/tf_dose")
outdir <- paste0(outdir_main, "/tf_dose")
dir.create(plotdir, recursive=T, showWarnings=F)
dir.create(outdir, recursive=T, showWarnings=F)

atacdiff_dir <- "../../output/02-atac/04/deseq_data_Plate_TF_Dose/"
meta <- read_tsv("atac_sample_sheet.txt")

celltype <- "HEK293T"
cond2 <- "GFP_1"
log2fc.cutoff <- 0.58
padj.cutoff <- 0.05

submeta <- meta[meta$CellType==celltype,]
cond1_ls <- unique(submeta$TF_Dose)
cond1_ls <- cond1_ls[cond1_ls != cond2]

peaktype_all <- list()
for (cond1 in cond1_ls){
  message(cond1)
  anno_out <- paste0(outdir, "/peaks_anno_diffpeak_", cond1, "_vs_", cond2, ".rds")
  if (!file.exists(anno_out)){
    atacdiff <- fread(Sys.glob(paste0(atacdiff_dir, "/deseq_*", cond1, "_vs_*", cond2, "_filteredConsensus_full.tsv")))
    atacdiff_sig <- atacdiff %>% dplyr::filter(padj<=padj.cutoff & abs(log2FoldChange)>=log2fc.cutoff)
    peaks <- GRanges(atacdiff_sig)
    if (length(peaks) > 1){
      peaks_anno <- ArchR:::.fastAnnoPeaks(peaks, BSgenome = BSgenome.Hsapiens.UCSC.hg38, geneAnnotation = geneAnnot)  
      saveRDS(peaks_anno, anno_out) 
    } else{
      next
    }
  } else{
    peaks_anno <- readRDS(anno_out)
  }
  
  
  if (length(peaks_anno) > 0){
    out <- peaks_anno %>% as.data.frame %>% dplyr::select(peakType) %>% table %>% as.data.frame() %>% dplyr::rename(peakType=".", count="Freq") %>% dplyr::mutate(TF_Dose=cond1, TF=strsplit(cond1, "_")[[1]][1],
                                                                                                                                                                  Plasmid_Dose=strsplit(cond1, "_")[[1]][2])
    peaktype_all[[cond1]] <- out  
  }
}

peaktype_all <- do.call(rbind, peaktype_all)
saveRDS(peaktype_all, paste0(outdir, "/peak_type_all.rds"))

p1 <- ggplot(peaktype_all, aes(x=Plasmid_Dose, y=count, fill=peakType)) + geom_col(position="stack") + theme(axis.text.x=element_text(angle=90, hjust=1)) + facet_wrap(~TF) 
p2 <- ggplot(peaktype_all, aes(x=Plasmid_Dose, y=count, fill=peakType)) + geom_col(position="fill") + ylab("count %") + theme(axis.text.x=element_text(angle=90, hjust=1)) + facet_wrap(~TF)
pdf(paste0(plotdir, "/peak_type_diffpeak_all_plasmid_dose.pdf"), width=10, height=10)
p1
p2
dev.off()

# peak type composition of differential peaks intersects, TF Dose ---------------------------------------------
plotdir <- paste0(plotdir_main, "/tf_dose_vennsets/")
outdir <- "../../output/02-atac/09/diff_peak_comp/" # this is output from the next script
dir.create(plotdir, recursive=T, showWarnings=F)
dir.create(outdir, recursive=T, showWarnings=F)

atacdiff_dir <- "../../output/02-atac/04/deseq_data_Plate_TF_Dose/"
cons <- readRDS(paste0(outdir_main, "/peaks_anno_all_consensus.rds"))
meta <- read_tsv("atac_sample_sheet.txt")

celltype <- "HEK293T"
log2fc.cutoff <- 0.58
padj.cutoff <- 0.05

submeta <- meta[meta$CellType==celltype,]
cond1_ls <- unique(submeta$TF)
cond1_ls <- cond1_ls[!cond1_ls %in% c("GFP", "negctl", "tetCTCF") ]

#for (i in seq_along(cond1_ls)){
for (i in 3:length(cond1_ls)){
## read differential data 
  TF <- cond1_ls[i]
  dosages <- c("1", "0.75", "0.5", "0.25","0.05")
  TF_dosages <- paste0(TF,"_",dosages)
  
  vennlist_file <- paste0(outdir, "/", TF, "_differential_peak_list.rds")
  if (!file.exists(vennlist_file)){
    a <- fread(Sys.glob(sprintf("%s/deseq_*%s_vs_*GFP_1_filteredConsensus_full.tsv", deseq_dir, TF_dosages[1]))) %>% dplyr::filter((abs(log2FoldChange)>log2fc.cutoff) & (padj<padj.cutoff))
    b <- fread(Sys.glob(sprintf("%s/deseq_*%s_vs_*GFP_1_filteredConsensus_full.tsv", deseq_dir, TF_dosages[2]))) %>% dplyr::filter((abs(log2FoldChange)>log2fc.cutoff) & (padj<padj.cutoff))
    c <- fread(Sys.glob(sprintf("%s/deseq_*%s_vs_*GFP_1_filteredConsensus_full.tsv", deseq_dir, TF_dosages[3]))) %>% dplyr::filter((abs(log2FoldChange)>log2fc.cutoff) & (padj<padj.cutoff))
    d <- fread(Sys.glob(sprintf("%s/deseq_*%s_vs_*GFP_1_filteredConsensus_full.tsv", deseq_dir, TF_dosages[4]))) %>% dplyr::filter((abs(log2FoldChange)>log2fc.cutoff) & (padj<padj.cutoff))
    e <- fread(Sys.glob(sprintf("%s/deseq_*%s_vs_*GFP_1_filteredConsensus_full.tsv", deseq_dir, TF_dosages[5]))) %>% dplyr::filter((abs(log2FoldChange)>log2fc.cutoff) & (padj<padj.cutoff))
    
    vennlist <- list(dose1=a$name, dose0.75=b$name, dose0.5=c$name, dose0.25=d$name, dose0.05=e$name)
    saveRDS(vennlist, paste0(outdir, "/", TF, "_differential_peak_list.rds"))
  } else{
    vennlist <- readRDS(vennlist_file)  
  }
  
  if (sum(lengths(vennlist) == 0)==5){
    next
  }
  # venn diagram of differential peaks 
  overlap <- get.venn.partitions(x=vennlist, force.unique = TRUE, keep.elements = TRUE,
                                 hierarchical = FALSE)
  
  # peak sets
  peakset1 <- overlap[overlap$`..set..`=="dose1∩dose0.75∩dose0.5∩dose0.25∩dose0.05", '..values..'] %>% unlist %>% unname
  subset1 <- cons[(cons$name %in% peakset1)]
  
  peakset2 <- overlap[overlap$`..set..`=='(dose1∩dose0.75∩dose0.5∩dose0.25)∖(dose0.05)', '..values..'] %>% unlist %>% unname
  subset2 <- cons[(cons$name %in% peakset2)]
  
  peakset3 <- overlap[overlap$`..set..`=='(dose1∩dose0.75∩dose0.5)∖(dose0.25∪dose0.05)', '..values..'] %>% unlist %>% unname
  subset3 <- cons[(cons$name %in% peakset3)]
  
  peakset4 <- overlap[overlap$`..set..`=='(dose1∩dose0.75)∖(dose0.5∪dose0.25∪dose0.05)', '..values..'] %>% unlist %>% unname
  subset4 <- cons[(cons$name %in% peakset4)]
  
  peakset5 <- overlap[overlap$`..set..`=='(dose1)∖(dose0.75∪dose0.5∪dose0.25∪dose0.05)', '..values..'] %>% unlist %>% unname
  subset5 <- cons[(cons$name %in% peakset5)]
  
  # construct the data frame for motif scores
  subset_ls <- list(subset1, subset2, subset3, subset4, subset5)
  df <- data.frame(peakType=lapply(subset_ls, function(n){n$peakType}) %>% unlist,
                   group=c(rep("d1&d0.75&d0.5&d0.25&d0.05",length(subset1)),
                           rep("(d1&d0.75&d0.5&d0.25)not(d0.05)", length(subset2)),
                           rep("(d1&d0.75&d0.5)not(d0.25|d0.05)", length(subset3)),
                           rep("(d1&d0.75)not(d0.5|d0.25|d0.05)", length(subset4)),
                           rep("(d1)not(d0.75|d0.5|d0.25|d0.05)", length(subset5))))
  
  p1 <- ggplot(df, aes(x=group, fill=peakType)) + geom_bar(position="stack") + theme_classic() + theme(axis.text.x=element_text(angle=90)) + ylab("count")
  p2 <- ggplot(df, aes(x=group, fill=peakType)) + geom_bar(position="fill")  + theme_classic() + theme(axis.text.x=element_text(angle=90)) + ylab("count %")

  
  pdf(paste0(plotdir, "/", TF, "_differential_venn_peaktype.pdf"), 
      width=6, height=6)
  lapply(paste0("p",1:2), function(n){print(get(n))})
  dev.off()
}
