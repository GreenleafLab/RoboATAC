library(tidyverse)
library(VennDiagram)
library(motifmatchr)
library(SummarizedExperiment)
library(ggpubr)
library(ggrepel)
library(ggrastr)
library(data.table)
library(ComplexHeatmap)


# 25 colors
colpal <- read_tsv("../utils/colpal_TF.tsv") %>% 
          column_to_rownames(var="TF")


# read TF independent inputs --------------------------------------------------------
cpm_df <- read.table("../../output/02-atac/01/cpm.tsv", row.names=1) # just use the cpm normalized peak count matrix

# read motif scores
motif.scores <- readRDS("../../output/02-atac/03/motif.scores.JASPAR2020.rds")

# add peak names 
cons <- readRDS("../../output/02-atac/08/peaks_anno_all_consensus.rds")
start(ranges(cons)) <- start(ranges(cons)) + 2
ov <- findOverlaps(cons, rowRanges(motif.scores), type="equal")

rowData(motif.scores) <- cons[ov@from[order(ov@to)]] %>% mcols
rownames(motif.scores) <- rowData(motif.scores)$name

# inputs -------------------------------------------------------------------------
# inputs for JASPAR motifs
motif_code_ls <- c("ALX4"="MA0681.2_PHOX2B",
                   "ELF1"="MA0473.3_ELF1",
                   "IRF4"="MA1419.1_IRF4",
                   "KLF1"="MA0039.4_KLF4",
                   "KLF4"="MA0039.4_KLF4",
                   "LEF1"="MA0768.1_LEF1",
                   "NR4A1"="MA1112.2_NR4A1",
                   "OCT4"="MA1115.1_POU5F1",
                   "PRDM1"="MA0508.3_PRDM1",
                   "SOX2"="MA0143.4_SOX2",
                   "SP4"="MA0685.1_SP4",
                   "SPI1"="MA0080.5_SPI1",
                   "TCF3"="MA0522.3_TCF3")

deseq_dir <- "../../output/02-atac/04/deseq_data_Plate_TF_Dose/"
plotdir <- "../../plots/02-atac/09/diff_peak_comp"
outdir <- "../../output/02-atac/09/diff_peak_comp"

dir.create(plotdir, showWarnings = F, recursive=T)
dir.create(outdir, showWarnings = F, recursive=T)


padj_cutoff <- 0.05
log2fc_cutoff <- 0.58

all_means <- list()
# loop through ----------------------------------------------------------
for (i in seq_along(motif_code_ls)){
  ## read differential data ----------------------------------------------------------------------
  TF <- names(motif_code_ls)[i]
  motif_code <- motif_code_ls[i]
  
  dosages <- c("1", "0.75", "0.5", "0.25","0.05")
  TF_dosages <- paste0(TF,"_",dosages)
  
  
  vennlist_file <- paste0(outdir, "/", TF, "_differential_peak_list.rds")
  if (!file.exists(vennlist_file)){
    a <- fread(Sys.glob(sprintf("%s/deseq_*%s_vs_*GFP_1_filteredConsensus_full.tsv", deseq_dir, TF_dosages[1]))) %>% dplyr::filter((abs(log2FoldChange)>log2fc_cutoff) & (padj<padj_cutoff))
    b <- fread(Sys.glob(sprintf("%s/deseq_*%s_vs_*GFP_1_filteredConsensus_full.tsv", deseq_dir, TF_dosages[2]))) %>% dplyr::filter((abs(log2FoldChange)>log2fc_cutoff) & (padj<padj_cutoff))
    c <- fread(Sys.glob(sprintf("%s/deseq_*%s_vs_*GFP_1_filteredConsensus_full.tsv", deseq_dir, TF_dosages[3]))) %>% dplyr::filter((abs(log2FoldChange)>log2fc_cutoff) & (padj<padj_cutoff))
    d <- fread(Sys.glob(sprintf("%s/deseq_*%s_vs_*GFP_1_filteredConsensus_full.tsv", deseq_dir, TF_dosages[4]))) %>% dplyr::filter((abs(log2FoldChange)>log2fc_cutoff) & (padj<padj_cutoff))
    e <- fread(Sys.glob(sprintf("%s/deseq_*%s_vs_*GFP_1_filteredConsensus_full.tsv", deseq_dir, TF_dosages[5]))) %>% dplyr::filter((abs(log2FoldChange)>log2fc_cutoff) & (padj<padj_cutoff))
    
    vennlist <- list(dose1=a$name, dose0.75=b$name, dose0.5=c$name, dose0.25=d$name, dose0.05=e$name)
    saveRDS(vennlist, paste0(outdir, "/", TF, "_differential_peak_list.rds"))
  } else{
    vennlist <- readRDS(vennlist_file)  
  }
  
  # venn diagram of differential peaks ----------------------------------------------------------------------
  overlap <- get.venn.partitions(x=vennlist, force.unique = TRUE, keep.elements = TRUE,
                                 hierarchical = FALSE)
  # overlap$`..set..`
  # 
  # overlap$`..count..`
  # 
  # rm(overrideTriple)
  # collist <- c("lightblue", "green", "red", "orange", "purple")
  # venn.diagram(x=vennlist, fill = collist,
  #              alpha = c(0.5, 0.5, 0.5, 0.5, 0.5), lwd =1, paste0(plotdir, "/", TF, "_differential_venn_diagram_equal.svg"), imagetype="svg",
  #              height=5, width=6, fontfamily = "sansserif", cat.fontfamily = "sansserif",resolution=300,
  #              disable.logging=T)
  # 

  # # this block doesnt work for more than 3 sets
  # overrideTriple <- T
  # venn.diagram(x=vennlist, fill = collist, 
  #              alpha = c(0.5, 0.5, 0.5, 0.5, 0.5), lwd =1, paste0(plotdir, "/", TF, "_differential_venn_diagram_proportional.svg"), imagetype="svg", 
  #              height=5, width=6, fontfamily = "sansserif", cat.fontfamily = "sansserif",resolution=300,
  #              disable.logging=T)
  
  
  # # upset plot ----------------------------------------------------------------------
  # m <- make_comb_mat(vennlist)
  # 
  # pdf(paste0(plotdir, "/", TF, "_upset.pdf"), width=8, height=4)
  # print(UpSet(m, comb_order = rev(order(comb_size(m))), row_title=paste0(TF, " differential peaks overlap")))
  # dev.off()

  # get motif scores and counts for peak sets ----------------------------------------
  # peaks that open at lowest dose
  peakset1 <- overlap[overlap$`..set..`=="dose1∩dose0.75∩dose0.5∩dose0.25∩dose0.05", '..values..'] %>% unlist %>% unname
  subset1 <- motif.scores[(rownames(motif.scores) %in% peakset1), motif_code]
  
  peakset2 <- overlap[overlap$`..set..`=='(dose1∩dose0.75∩dose0.5∩dose0.25)∖(dose0.05)', '..values..'] %>% unlist %>% unname
  subset2 <- motif.scores[(rownames(motif.scores) %in% peakset2), motif_code]
  
  peakset3 <- overlap[overlap$`..set..`=='(dose1∩dose0.75∩dose0.5)∖(dose0.25∪dose0.05)', '..values..'] %>% unlist %>% unname
  subset3 <- motif.scores[(rownames(motif.scores) %in% peakset3), motif_code]
  
  peakset4 <- overlap[overlap$`..set..`=='(dose1∩dose0.75)∖(dose0.5∪dose0.25∪dose0.05)', '..values..'] %>% unlist %>% unname
  subset4 <- motif.scores[(rownames(motif.scores) %in% peakset4), motif_code]
  
  peakset5 <- overlap[overlap$`..set..`=='(dose1)∖(dose0.75∪dose0.5∪dose0.25∪dose0.05)', '..values..'] %>% unlist %>% unname
  subset5 <- motif.scores[(rownames(motif.scores) %in% peakset5), motif_code]
  
  # # some optional quick checks
  # subset1
  # 
  # which(motifScores(subset1)>0) %>% length
  # 
  # genelist1 <- a[a$name %in% peakset1, "gene_name"] %>% table() %>% sort(decreasing = T)
  # genelist1 %>% head(20)
  
  # construct the data frame for motif scores
  subset_ls <- list(subset1, subset2, subset3, subset4, subset5)
  df <- data.frame(score=lapply(subset_ls, function(n){motifScores(n) %>% as.numeric}) %>% unlist,
                   count=lapply(subset_ls, function(n){motifCounts(n) %>% as.numeric}) %>% unlist,  
                   peakType=lapply(subset_ls, function(n){rowData(n)$peakType}) %>% unlist,
                   group=c(rep("d1&d0.75&d0.5&d0.25&d0.05",length(subset1)),
                           rep("(d1&d0.75&d0.5&d0.25)not(d0.05)", length(subset2)),
                           rep("(d1&d0.75&d0.5)not(d0.25|d0.05)", length(subset3)),
                           rep("(d1&d0.75)not(d0.5|d0.25|d0.05)", length(subset4)),
                           rep("(d1)not(d0.75|d0.5|d0.25|d0.05)", length(subset5))))
  dfwithmotif <- df[df$score != 0,]
# 
#   ggplot(dfwithmotif, aes(x=score, color=group)) + stat_ecdf(geom = "step") + theme_classic()+ ylab("cumulative density") + xlab("motif score")
#   ggplot(dfwithmotif, aes(x=count, color=group)) + stat_ecdf(geom = "step") + theme_classic()+ ylab("cumulative density") + xlab("motif count")
# 
#   p1 <- ggviolin(df, y="score", fill = "group", draw_quantiles = 0.5)
#   p2 <- ggdensity(df, x="score", fill="group", facet.by="group")
#   p3 <- gghistogram(df, x="score", y="density", fill="group", facet.by="group")
#   p4 <- ggviolin(dfwithmotif, y="score", fill = "group",draw_quantiles = 0.5) + ylab("score (non-zero only)")
#   p5 <- ggboxplot(dfwithmotif, y="score", fill = "group") + ylab("score (non-zero only)")
#   p6 <- ggdensity(dfwithmotif, x="score", fill="group", facet.by="group") + xlab("score (non-zero only)")
#   p7 <- gghistogram(dfwithmotif, x="score", y="density", fill="group", facet.by="group") + xlab("score (non-zero only)")
# 
#   p8 <- ggviolin(df, y="count", fill = "group",draw_quantiles = 0.5)
#   p9 <- ggboxplot(df, y="count", fill = "group")
#   p10 <- ggdensity(df, x="count", fill="group", facet.by="group")
#   p11 <- gghistogram(df, x="count", y="density", fill="group", facet.by="group")
#   p12 <- ggviolin(dfwithmotif, y="count", fill = "group",draw_quantiles = 0.5) + ylab("count (non-zero only)")
#   p13 <- ggboxplot(dfwithmotif, y="count", fill = "group") + ylab("count (non-zero only)")
#   p14 <- ggdensity(dfwithmotif, x="count", fill="group", facet.by="group") + xlab("count (non-zero only)")
#   p15 <- gghistogram(dfwithmotif, x="count", y="density", fill="group", facet.by="group") + xlab("count (non-zero only)")
# 
#   pdf(paste0(plotdir, "/", TF, "_differential_venn_motifscore_motifcount.pdf"),
#       width=8, height=4)
#   lapply(paste0("p",1:15), function(n){print(get(n))})
#   dev.off()
# 
#   # peak types
#   p16 <- ggplot(df, aes(x=group, fill=peakType)) + geom_bar(position="stack") + theme_classic() + theme(axis.text.x=element_text(angle=90)) + ylab("count")
#   p17 <- ggplot(df, aes(x=group, fill=peakType)) + geom_bar(position="fill")  + theme_classic() + theme(axis.text.x=element_text(angle=90)) + ylab("count %")
#   p18 <- ggplot(dfwithmotif, aes(x=group, fill=peakType)) + geom_bar(position="stack") + theme_classic() + theme(axis.text.x=element_text(angle=90)) + ylab("count (non-zero motif count only)")
#   p19 <- ggplot(dfwithmotif, aes(x=group, fill=peakType)) + geom_bar(position="fill")  + theme_classic() + theme(axis.text.x=element_text(angle=90)) + ylab("count % (non-zero motif count only)")
#   
#   # does promoter tend to have more/less better/worse motifs than distal regions for this TF? 
#   subset <- motif.scores[,motif_code]
#   all <- data.frame(count=motifCounts(subset) %>% as.numeric, score=motifScores(subset)%>% as.numeric, peakType=rowData(subset)$peakType)
#   allwithmotif <- all[all$score != 0,]
#   p20 <- ggplot(all, aes(x=peakType, y=count)) + rasterise(geom_boxplot(), dpi=150) + theme_classic()
#   p21 <- ggplot(all, aes(x=peakType, y=score)) + rasterise(geom_boxplot(), dpi=150) + theme_classic()
#   p22 <- ggplot(allwithmotif, aes(x=peakType, y=count)) + rasterise(geom_boxplot(), dpi=150) + theme_classic() + ylab("motif count (nonzero only)")
#   p23 <- ggplot(allwithmotif, aes(x=peakType, y=score)) + rasterise(geom_boxplot(), dpi=150) + theme_classic() + ylab("motif score (nonzero only)")
#   
#   pdf(paste0(plotdir, "/", TF, "_differential_venn_peaktype.pdf"), 
#       width=6, height=6)
#   lapply(paste0("p",16:23), function(n){print(get(n))})
#   dev.off()
  

  # prettier version summary plot
  mean_data <- df %>%
    group_by(group) %>%
    summarise(mean_score = mean(score),
              mean_count = mean(count),
              npeak_Promoter = sum(peakType=="Promoter"),
              npeak_Distal = sum(peakType=="Distal"),
              npeak_Exonic = sum(peakType=="Exonic"),
              npeak_Intronic = sum(peakType=="Intronic")) %>%
    mutate(pct_Promoter=npeak_Promoter/(npeak_Promoter + npeak_Distal + npeak_Exonic + npeak_Intronic) * 100,
           pct_Distal=npeak_Distal/(npeak_Promoter + npeak_Distal + npeak_Exonic + npeak_Intronic) * 100,
           pct_Exonic=npeak_Exonic/(npeak_Promoter + npeak_Distal + npeak_Exonic + npeak_Intronic) * 100,
           pct_Intronic=npeak_Intronic/(npeak_Promoter + npeak_Distal + npeak_Exonic + npeak_Intronic) * 100)
  
  mean_data_nonzero <- dfwithmotif %>%
    group_by(group) %>%
    summarise(mean_score = mean(score),
              mean_count = mean(count),
              npeak_Promoter = sum(peakType=="Promoter"),
              npeak_Distal = sum(peakType=="Distal"),
              npeak_Exonic = sum(peakType=="Exonic"),
              npeak_Intronic = sum(peakType=="Intronic")) %>%
    mutate(pct_Promoter=npeak_Promoter/(npeak_Promoter + npeak_Distal + npeak_Exonic + npeak_Intronic) * 100,
           pct_Distal=npeak_Distal/(npeak_Promoter + npeak_Distal + npeak_Exonic + npeak_Intronic) * 100,
           pct_Exonic=npeak_Exonic/(npeak_Promoter + npeak_Distal + npeak_Exonic + npeak_Intronic) * 100,
           pct_Intronic=npeak_Intronic/(npeak_Promoter + npeak_Distal + npeak_Exonic + npeak_Intronic) * 100)
  
  colnames(mean_data_nonzero) <- paste0(colnames(mean_data_nonzero), "_nonzero") 
  mean_data_nonzero <- mean_data_nonzero %>% dplyr::rename(group="group_nonzero")
    
  means <- merge(mean_data, mean_data_nonzero, by="group") %>% mutate(TF=TF)
  
  means <- means %>% mutate(mean_score_scale=scales::rescale(mean_score),
                            mean_score_nonzero_scale=scales::rescale(mean_score_nonzero),
                            mean_count_scale=scales::rescale(mean_count),
                            mean_count_nonzero_scale=scales::rescale(mean_count_nonzero),
                            pct_Promoter_scale=scales::rescale(pct_Promoter),
                            pct_Distal_scale=scales::rescale(pct_Distal),
                            pct_Exonic_scale=scales::rescale(pct_Exonic),
                            pct_Intronic_scale=scales::rescale(pct_Intronic),
                            pct_Promoter_nonzero_scale=scales::rescale(pct_Promoter_nonzero),
                            pct_Distal_nonzero_scale=scales::rescale(pct_Distal_nonzero),
                            pct_Exonic_nonzero_scale=scales::rescale(pct_Exonic_nonzero),
                            pct_Intronic_nonzero_scale=scales::rescale(pct_Intronic_nonzero))
  all_means[[i]] <- means
  
  # # score
  # p16 <- ggplot(df, aes(x = score, fill=group, color=group)) +
  #   geom_density(adjust = 2, alpha=0.3)  + ylab(paste0(motif_code, " motif score")) +
  #   geom_vline(data = mean_data, aes(xintercept = mean_score, color=group), linetype="dashed") +
  #   geom_text_repel(data = mean_data, aes(label = round(mean_score, 2), x = mean_score, y = 0)) +
  #   theme_classic() + theme(legend.position = "right", text=element_text(size=14))
  # 
  # p17 <- ggplot(dfwithmotif, aes(x = score, fill=group, color=group)) +
  #   geom_density(adjust = 2, alpha=0.3)  + ylab(paste0(motif_code, " motif score (nonzero only)")) +
  #   geom_vline(data = mean_data_nonzero, aes(xintercept = mean_score_nonzero, color=group), linetype="dashed") +
  #   geom_text_repel(data = mean_data_nonzero, aes(label = round(mean_score_nonzero, 2), x = mean_score_nonzero, y = 0)) +
  #   theme_classic() + theme(legend.position = "right", text=element_text(size=14))
  # 
  # # counts
  # p18 <- ggplot(df, aes(x = count, fill=group, color=group)) +
  #   geom_density(adjust = 3, alpha=0.3) + xlim(c(-1, 10)) + ylab(paste0(motif_code, " motif count")) +
  #   geom_vline(data = mean_data, aes(xintercept = mean_count, color=group), linetype="dashed") +
  #   geom_text_repel(data = mean_data, aes(label = round(mean_count, 2), x = mean_count, y = 0)) +
  #   theme_classic() + theme(legend.position = "right", text=element_text(size=14))
  # 
  # p19 <- ggplot(dfwithmotif, aes(x = count, fill=group, color=group)) +
  #   geom_density(adjust = 3, alpha=0.3) + xlim(c(-1, 10)) + ylab(paste0(motif_code, " motif count (nonzero only)")) +
  #   geom_vline(data = mean_data_nonzero, aes(xintercept = mean_count_nonzero, color=group), linetype="dashed") +
  #   geom_text_repel(data = mean_data_nonzero, aes(label = round(mean_count_nonzero, 2), x = mean_count_nonzero, y = 0)) +
  #   theme_classic() + theme(legend.position = "right", text=element_text(size=14))
  # 
  # 
  # pdf(paste0(plotdir, "/", TF, "_differential_venn_motifscore_motifcount_pretty.pdf"),
  #     width= 8, height=4)
  # lapply(paste0("p",16:19), function(n){print(get(n))})
  # dev.off()

  # # correlate peak CPM with motif scores and motif counts ----------------------------------
  # cpm_dfsub <- cpm_df[, grep(paste0(TF,"|GFP"),colnames(cpm_df))]
  # cpm_dfsub %>% head
  # 
  # df <- data.frame(motifScore=motifScores(subset1[peakset1]) %>% as.numeric,
  #                  motifCount=motifCounts(subset1[peakset1]) %>% as.numeric,
  #                  dose0cpm=cpm_dfsub[peakset1,grep("GFP",colnames(cpm_dfsub))] %>% rowMeans(na.rm=T),
  #                  dose005cpm=cpm_dfsub[peakset1,grep("005",colnames(cpm_dfsub))] %>% rowMeans(na.rm=T),
  #                  dose025cpm=cpm_dfsub[peakset1,grep("025",colnames(cpm_dfsub))] %>% rowMeans(na.rm=T),
  #                  dose050cpm=cpm_dfsub[peakset1,grep("050",colnames(cpm_dfsub))] %>% rowMeans(na.rm=T),
  #                  dose075cpm=cpm_dfsub[peakset1,grep("075",colnames(cpm_dfsub))] %>% rowMeans(na.rm=T),
  #                  dose100cpm=cpm_dfsub[peakset1,grep("100",colnames(cpm_dfsub))] %>% rowMeans(na.rm=T))
  # df %>% head
  # 
  # p1 <- ggscatter(data=df[df$motifScore>0,], y="dose0cpm", x="motifScore") + yscale("log10")
  # p2 <- ggscatter(data=df[df$motifScore>0,], y="dose005cpm", x="motifScore") + yscale("log10")
  # p3 <- ggscatter(data=df[df$motifScore>0,], y="dose025cpm", x="motifScore") + yscale("log10")
  # p4 <- ggscatter(data=df[df$motifScore>0,], y="dose050cpm", x="motifScore") + yscale("log10")
  # p5 <- ggscatter(data=df[df$motifScore>0,], y="dose075cpm", x="motifScore") + yscale("log10")
  # p6 <- ggscatter(data=df[df$motifScore>0,], y="dose100cpm", x="motifScore") + yscale("log10")
  # 
  # p7 <- ggboxplot(data=df[df$motifScore>0,], y="dose0cpm", x="motifCount") + yscale("log10")
  # p8 <- ggboxplot(data=df[df$motifScore>0,], y="dose005cpm", x="motifCount") + yscale("log10")
  # p9 <- ggboxplot(data=df[df$motifScore>0,], y="dose025cpm", x="motifCount") + yscale("log10")
  # p10 <- ggboxplot(data=df[df$motifScore>0,], y="dose050cpm", x="motifCount") + yscale("log10")
  # p11 <- ggboxplot(data=df[df$motifScore>0,], y="dose075cpm", x="motifCount") + yscale("log10")
  # p12 <- ggboxplot(data=df[df$motifScore>0,], y="dose100cpm", x="motifCount") + yscale("log10")
  # 
  # 
  # pdf(paste0(plotdir, "/", TF, "_correlate_peakcpm_vs_motifscore_motifcount.pdf"),
  #     width=6, height=6)
  # lapply(paste0("p",1:12), function(n){print(get(n))})
  # dev.off()


}


# plot all means -------------------------------------------------
all_means <- do.call(rbind, all_means)

saveRDS(all_means, paste0(outdir, "/all_means.rds"))
all_means$group <- factor(all_means$group, levels=c("d1&d0.75&d0.5&d0.25&d0.05", "(d1&d0.75&d0.5&d0.25)not(d0.05)", "(d1&d0.75&d0.5)not(d0.25|d0.05)", 
                                                    "(d1&d0.75)not(d0.5|d0.25|d0.05)", "(d1)not(d0.75|d0.5|d0.25|d0.05)"))

all_means_sub <- all_means[!all_means$TF %in% c("LEF1", "NR4A1", "ELF1", "PRDM1"),]
all_means_sub2 <- all_means[!all_means$TF %in% c("KLF1", "KLF4", "SP4"),]
all_means_sub3 <- all_means[!all_means$TF %in% c("KLF1", "KLF4", "SP4","LEF1", "NR4A1", "ELF1", "PRDM1"),]

#plot_df <- all_means
#plot_df <- all_means_sub
#plot_df <- all_means_sub2
plot_df <- all_means_sub3
BL_stallion <- colpal[unique(plot_df$TF), "color"]
  
p1 <- ggplot(plot_df, aes(x=group, y=mean_count, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1)) + scale_color_manual(values=BL_stallion)
p2 <- ggplot(plot_df, aes(x=group, y=mean_count_nonzero, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)
p3 <- ggplot(plot_df, aes(x=group, y=mean_count_scale, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)
p4 <- ggplot(plot_df, aes(x=group, y=mean_count_nonzero_scale, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)

p5 <- ggplot(plot_df, aes(x=group, y=mean_score, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)
p6 <- ggplot(plot_df, aes(x=group, y=mean_score_nonzero, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)
p7 <- ggplot(plot_df, aes(x=group, y=mean_score_scale, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)
p8 <- ggplot(plot_df, aes(x=group, y=mean_score_nonzero_scale, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)

p9 <- ggplot(plot_df, aes(x=group, y=pct_Promoter, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)
p10 <- ggplot(plot_df, aes(x=group, y=pct_Promoter_nonzero, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)
p11 <- ggplot(plot_df, aes(x=group, y=pct_Promoter_scale, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)
p12 <- ggplot(plot_df, aes(x=group, y=pct_Promoter_nonzero_scale, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)

p13 <- ggplot(plot_df, aes(x=group, y=pct_Distal, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)
p14 <- ggplot(plot_df, aes(x=group, y=pct_Distal_nonzero, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)
p15 <- ggplot(plot_df, aes(x=group, y=pct_Distal_scale, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)
p16 <- ggplot(plot_df, aes(x=group, y=pct_Distal_nonzero_scale, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)

p17 <- ggplot(plot_df, aes(x=group, y=pct_Exonic, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)
p18 <- ggplot(plot_df, aes(x=group, y=pct_Exonic_nonzero, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)
p19 <- ggplot(plot_df, aes(x=group, y=pct_Exonic_scale, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)
p20 <- ggplot(plot_df, aes(x=group, y=pct_Exonic_nonzero_scale, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)

p21 <- ggplot(plot_df, aes(x=group, y=pct_Intronic, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)
p22 <- ggplot(plot_df, aes(x=group, y=pct_Intronic_nonzero, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)
p23 <- ggplot(plot_df, aes(x=group, y=pct_Intronic_scale, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)
p24 <- ggplot(plot_df, aes(x=group, y=pct_Intronic_nonzero_scale, color=TF, group=TF)) + geom_line() + geom_point(size=2) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1))  + scale_color_manual(values=BL_stallion)


#pdf(paste0(plotdir, "/all_means.pdf"), height=6, width=6)
#pdf(paste0(plotdir, "/all_means_noLEF1_NR4A1_ELF1_PRDM1.pdf"), height=6, width=6)
#pdf(paste0(plotdir, "/all_means_noKLF1_KLF4_SP4.pdf"), height=6, width=6)
pdf(paste0(plotdir, "/all_means_noLEF1_NR4A1_ELF1_PRDM1_KLF1_KLF4_SP4.pdf"), height=6, width=6)
lapply(paste0("p",1:24), function(n){eval(parse(text=n))})
dev.off()


