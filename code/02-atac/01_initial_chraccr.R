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
  library(rtracklayer)
  library(data.table)
})

source("../utils/tf_funcs.R")

outdir <- "../../output/02-atac/01/"
plotdir <- "../../plots/02-atac/01/"
dir.create(plotdir, recursive=T, showWarnings=F)
dir.create(outdir, recursive=T, showWarnings=F)

# 25 colors
colpal <- read_tsv("../utils/colpal_TF.tsv") %>% 
  column_to_rownames(var="TF")

cmap_chromvar <- c("#3361A5", "#248AF3", "#14B3FF", "#88CEEF", "#C1D5DC", "#EAD397", "#FDB31A", "#E42A2A", "#A31D1D")

# Make dsa bam object -----------------------------------------------------
sampleAnnot <- read.csv("atac_sample_sheet.txt", sep="\t")
sampleAnnot <- sampleAnnot[sampleAnnot$bamFilenameFull != "#N/A",]
sampleAnnot <- sampleAnnot[sampleAnnot$CellType == "HEK293T",]


setConfigElement("annotationColumns", c("Plasmid_Dose", "TF"))
setConfigElement("filteringCovgCount", 1L)
setConfigElement("filteringCovgReqSamples", 0.25)
setConfigElement("filteringSexChroms", TRUE)
setConfigElement("normalizationMethod", "quantile")

consensuspeaks <- getPeakSet.snakeATAC.v2(sampleAnnot, filePrefixCol = "peakFilename", 
                          genome = "hg38", dataDir = "", sampleIdCol = "sampleName", 
                          replicateCol="TF_Dose", capMaxPeaks=150000)

# Filter dsa to remove the peaks corresponding to the transfected ----------
#exclude regions overlapping with gene body of the overexpressed TFs
# for ultima data, when sequenced deep enough, the reads are not dominated by these plasmid peaks, 
#                   no obvious peak to remove, we will still filter out these regions
exclTFs <- lapply(list.files("exclude_blacklist", full=T), read.csv)
exclTFs <- do.call(rbind, exclTFs)

blacklist <- GRanges(seqnames=exclTFs$chr, 
                     ranges=IRanges(start=exclTFs$start,end=exclTFs$end), 
                     strand = exclTFs$strand)

overlap <- findOverlaps(blacklist, consensuspeaks) # col 1 index in blacklist, col 2 index in consensuspeaks
overlap_id <- subjectHits(overlap) # List of all the indices in consensuspeaks
all_id <- 1:length(consensuspeaks) # A list of all indices for consensus peaks 1 to the length
consensuspeaks.f <- consensuspeaks[!(all_id %in% overlap_id)] # Subset concensuspeaks

dsa <- DsATAC.bam(sampleAnnot, "bamFilenameFull", "hg38", regionSets = list("filteredConsensus" = consensuspeaks.f), sampleIdCol="sampleName", 
                          pairedEnd=FALSE)
saveDsAcc(dsa, paste0(outdir, "/dsa_filt"))
gr <- getCoord(dsa,"filteredConsensus")

# export consensus peaks as bed
export.bed(gr, con=paste0(outdir,'/consensus_peaks_HEK293T.bed'))
bed <- read.table(paste0(outdir, '/consensus_peaks_HEK293T.bed'))

# add extra columns to make a 10 column bed file for chrombpnet model input
bed$V7 <- "."
bed$V8 <- "."
bed$V9 <- "."
bed$V10 <- 250
write.table(bed, paste0(outdir, '/consensus_peaks_HEK293T_10col.bed'), quote=F, sep="\t", row.names=F, col.names=F)

# Normalization --------------------------------------------------
dsa <- loadDsAcc(paste0(outdir, "/dsa_filt"))

dsa_norm <- transformCounts(dsa, method = "CPM")
saveDsAcc(dsa_norm, paste0(outdir, "/dsa_norm_cpm"))
write.table(data.frame(dsa_norm@counts$filteredConsensus, 
                       row.names = dsa_norm@coord$filteredConsensus$name), 
                     paste0(outdir, "/cpm.tsv"), row.names=TRUE)

dsa_norm <- transformCounts(dsa, method = "quantile")
saveDsAcc(dsa_norm, paste0(outdir, "/dsa_norm_quantile"))

dsa_norm <- transformCounts(dsa, method = "RPKM")
dsa_norm <- transformCounts(dsa_norm, method = "log2")
dsa_norm <- transformCounts(dsa_norm, method = "quantile")
saveDsAcc(dsa_norm, paste0(outdir, "/dsa_norm_rpkmlog2quantile"))

dsa_norm <- transformCounts(dsa, method = "vst")
saveDsAcc(dsa_norm, paste0(outdir, "/dsa_norm_vst"))


# PCA ---------------------------------------------------------------------
# now in pca.R
command <- "sbatch -p wjg,sfgf,biochem --mem-per-cpu=200g --time=24:00:00 --nodes=2 --job-name=pca_hek --out=slurm-logs/slurm-%j-hek_pca.out --wrap \"Rscript 02_pca.R\""
system(command)

# motif matching ---------------------------------------------------------------------
command <- "sbatch -p wjg,sfgf,biochem --mem-per-cpu=200g --time=24:00:00 --nodes=2 --job-name=hek_motif --out=slurm-logs/slurm-%j-hek_motif.out --wrap \"Rscript 03_motif_matching.R\""
system(command)

# Create ChrAccR Reports --------------------------------------------------
reportDir <- paste0(outdir, "/ChrAccR_reports")

dsa_norm <- loadDsAcc(paste0(outdir, '/dsa_norm_rpkmlog2quantile/'))
setConfigElement("chromVarRegionTypes", c("filteredConsensus"))
setConfigElement("chromVarMotifs", getMotifs(database=JASPAR2020::JASPAR2020))
createReport_exploratory(dsa_norm, reportDir)

dsa <- loadDsAcc(paste0(outdir, "/dsa_filt"))
setConfigElement("regionTypes",c("filteredConsensus"))
setConfigElement("differentialColumns", c("TF"))
createReport_differential(dsa, reportDir)


# replot chromvar heatmap ------------------------------------------------
chromvar <- readRDS(paste0(reportDir, "/exploratory_data/chromVarDev_filteredconsensus.rds"))
library(ComplexHeatmap)

# filter for top 100 variable motifs
vars <- chromVAR::computeVariability(chromvar) %>% dplyr::arrange(desc(variability))
top100 <- vars %>% head(100) %>% rownames

# # alternatively, filter for top 100 non-JUN/FOS motifs
#top100 <- vars[!grepl("JUN|FOS", rownames(vars)),] %>% head(100) %>% rownames
# top100 <- vars[!grepl("JUN|FOS", rownames(vars)),] %>% head(200) %>% rownames

mostVarIdx <- which(rownames(chromvar@assays@data$deviations) %in% top100)
match(top100, rownames(chromvar@assays@data$deviations))

ha <- columnAnnotation(TF=chromvar$TF, Plasmid_Dose=chromvar$Plasmid_Dose, PlateNum=as.character(chromvar$PlateNum),
                       col=list(TF=setNames(colpal$color, rownames(colpal)),
                                PlateNum=setNames(c("#F8766D", "#CD9600", "#00B8E7"), c("1", "2", "3")),
                                Plasmid_Dose=setNames(viridisLite::viridis(6), c("0", "0.05", "0.25", "0.5", "0.75", "1"))))


mat <- chromvar@assays@data$z[mostVarIdx,]

color_mapping_1 <- cmap_chromvar
break_points <- quantile(mat %>% as.matrix, probs = seq(0.001, 0.999, length.out = length(cmap_chromvar)))
color_mapping_2 <- circlize::colorRamp2(break_points, cmap_chromvar)
color_mapping_3 <- circlize::colorRamp2(seq(-150, 150, length.out = length(cmap_chromvar)), cmap_chromvar)
color_mapping_4 <- circlize::colorRamp2(seq(-100, 100, length.out = length(cmap_chromvar)), cmap_chromvar)
#color_mapping_4 <- circlize::colorRamp2(seq(-20, 20, length.out = length(cmap_chromvar)), cmap_chromvar)

ht1 <- Heatmap(mat,
      name="ChromVAR \ndeviations \n(z-score)",
      top_annotation = ha,
      cluster_rows=T,
      cluster_columns=T,
      column_names_gp = gpar(fontsize = 2),
      column_labels = chromvar$sampleName,
      row_title = paste0(dim(mat)[1], " top variable motifs"),
      row_names_gp = gpar(fontsize=6),
      row_labels = rownames(mat),
      column_title = paste0(dim(mat)[2], " samples"),
      col = color_mapping_1
)

ht2 <- Heatmap(mat,
      name="ChromVAR \ndeviations \n(z-score)",
      top_annotation = ha,
      cluster_rows=T,
      cluster_columns=T,
      column_names_gp = gpar(fontsize = 2),
      column_labels = chromvar$sampleName,
      row_title = paste0(dim(mat)[1], " top variable motifs"),
      row_names_gp = gpar(fontsize=6),
      row_labels = rownames(mat),
      column_title = paste0(dim(mat)[2], " samples"),
      col = color_mapping_2
)

ht3 <- Heatmap(mat,
      name="ChromVAR \ndeviations \n(z-score)",
      top_annotation = ha,
      cluster_rows=T,
      cluster_columns=T,
      column_names_gp = gpar(fontsize = 2),
      column_labels = chromvar$sampleName,
      row_title = paste0(dim(mat)[1], " top variable motifs"),
      row_names_gp = gpar(fontsize=6),
      row_labels = rownames(mat),
      column_title = paste0(dim(mat)[2], " samples"),
      col = color_mapping_3
)

ht4 <- Heatmap(mat,
      name="ChromVAR \ndeviations \n(z-score)",
      top_annotation = ha,
      cluster_rows=T,
      cluster_columns=T,
      column_names_gp = gpar(fontsize = 2),
      column_labels = chromvar$sampleName,
      row_title = paste0(dim(mat)[1], " top variable motifs"),
      row_names_gp = gpar(fontsize=6),
      #row_names_gp = gpar(fontsize=3),
      row_labels = rownames(mat),
      column_title = paste0(dim(mat)[2], " samples"),
      col = color_mapping_4
)

pdf(paste0(plotdir, "/chromvar_heatmap.pdf"), height=10, width=8)
#pdf(paste0(plotdir, "/chromvar_heatmap_noJUNFOS.pdf"), height=10, width=8)
# pdf(paste0(plotdir, "/chromvar_heatmap_top200_noJUNFOS.pdf"), height=10, width=8)
plot(ht1)
plot(ht2)
plot(ht3)
plot(ht4)
dev.off()


# replot chromvar heatmap per plate ------------------------------------------------
for (plate in unique(chromvar$PlateNum)){
      row_textsize <- 3
      top_n <- 200
      excl_JUNFOS <- T
      #outname <- paste0(plotdir, "/chromvar_heatmap_plate", plate, ".pdf")
      #outname <- paste0(plotdir, "/chromvar_heatmap_noJUNFOS_plate", plate, ".pdf")
      outname <- paste0(plotdir, "/chromvar_heatmap_top200_noJUNFOS_plate", plate, ".pdf")

      # filter for top 100 variable motifs
      subchromvar <- chromvar[,chromvar$PlateNum==plate]

      # recompute z scores
      subchromvar@assays@data$z <- subchromvar@assays@data$deviations %>% t %>% scale %>% t
      vars <- chromVAR::computeVariability(subchromvar) %>% dplyr::arrange(desc(variability))
      if (excl_JUNFOS){
            # filter for top non-JUN/FOS motifs
            top100 <- vars[!grepl("JUN|FOS", rownames(vars)),] %>% head(top_n) %>% rownames
      } else{
            top100 <- vars %>% head(top_n) %>% rownames
      }
         
      mostVarIdx <- which(rownames(subchromvar@assays@data$deviations) %in% top100)
      #match(top100, rownames(subchromvar@assays@data$deviations))

      names(BL_stallion) <- subchromvar$TF %>% unique

      ha <- columnAnnotation(TF=subchromvar$TF, Plasmid_Dose=subchromvar$Plasmid_Dose, PlateNum=as.character(subchromvar$PlateNum),
                        col=list(TF=BL_stallion[1:length(subchromvar$TF %>% unique)],
                                    PlateNum=setNames(viridisLite::viridis(3), c("1", "2", "3"))))


      mat <- subchromvar@assays@data$z[mostVarIdx,]

      color_mapping_1 <- cmap_chromvar
      break_points <- quantile(mat %>% as.matrix, probs = seq(0.001, 0.999, length.out = length(cmap_chromvar)))
      color_mapping_2 <- circlize::colorRamp2(break_points, cmap_chromvar)
      color_mapping_3 <- circlize::colorRamp2(seq(-150, 150, length.out = length(cmap_chromvar)), cmap_chromvar)
      color_mapping_4 <- circlize::colorRamp2(seq(-5, 5, length.out = length(cmap_chromvar)), cmap_chromvar)

      ht1 <- Heatmap(mat,
            name="ChromVAR \ndeviations \n(z-score)",
            top_annotation = ha,
            cluster_rows=T, cluster_columns=T,
            column_names_gp = gpar(fontsize = 2), column_labels = subchromvar$sampleName,
            row_title = paste0(dim(mat)[1], " top variable motifs"),
            row_names_gp = gpar(fontsize=row_textsize),
            row_labels = rownames(mat),
            column_title = paste0(dim(mat)[2], " samples"),
            col = color_mapping_1
      )

      ht2 <- Heatmap(mat,
            name="ChromVAR \ndeviations \n(z-score)",
            top_annotation = ha,
            cluster_rows=T, cluster_columns=T,
            column_names_gp = gpar(fontsize = 2), column_labels = subchromvar$sampleName,
            row_title = paste0(dim(mat)[1], " top variable motifs"),
            row_names_gp = gpar(fontsize=row_textsize),
            row_labels = rownames(mat),
            column_title = paste0(dim(mat)[2], " samples"),
            col = color_mapping_2
      )

      ht3 <- Heatmap(mat,
            name="ChromVAR \ndeviations \n(z-score)",
            top_annotation = ha,
            cluster_rows=T, cluster_columns=T,
            column_names_gp = gpar(fontsize = 2), column_labels = subchromvar$sampleName,
            row_title = paste0(dim(mat)[1], " top variable motifs"),
            row_names_gp = gpar(fontsize=row_textsize),
            row_labels = rownames(mat),
            column_title = paste0(dim(mat)[2], " samples"),
            col = color_mapping_3
      )

      ht4 <- Heatmap(mat,
            name="ChromVAR \ndeviations \n(z-score)",
            top_annotation = ha,
            cluster_rows=T, cluster_columns=T,
            column_names_gp = gpar(fontsize = 2), column_labels = subchromvar$sampleName,
            row_title = paste0(dim(mat)[1], " top variable motifs"),
            row_names_gp = gpar(fontsize=row_textsize),
            row_labels = rownames(mat),
            column_title = paste0(dim(mat)[2], " samples"),
            col = color_mapping_4
      )

      pdf(outname, height=10, width=8)
      plot(ht1)
      plot(ht2)
      plot(ht3)
      plot(ht4)
      dev.off()

}