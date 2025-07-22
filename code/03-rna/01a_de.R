# Differential analysis of RNAseq counts estimated by kallisto
###################### import packages ######################
suppressWarnings({
  library(DESeq2)
  library(biomaRt)
  library(sleuth)
  library(tidyr)
  library(tximport)
  library(RColorBrewer)
  library(viridis)
  library(ggplot2)
  library(ggpubr)
  library(ComplexHeatmap)
  library(EnhancedVolcano)
  library(ggrastr)
  library(tidyverse)
  library(data.table)
})

# color palette ---------------------------------------------------------
# 25 colors
colpal <- read_tsv("../utils/colpal_TF.tsv") %>% 
  column_to_rownames(var="TF")

source("../utils/tf_funcs.R")

###################### read data ######################
outdir <- "../../output/03-rna/01/"
plotdir <- "../../plots/03-rna/01/"

dir.create(outdir, recursive=T, showWarnings=F)
dir.create(plotdir, recursive=T, showWarnings=F)

s2c <- read_tsv("rna_sample_sheet.txt")

# change names to be consistent with ATAC
s2c[s2c$TF=="EF1a-SPI1", "TF"] <- "SPI1"
s2c[s2c$TF=="EF1a-KLF1", "TF"] <- "KLF1"
s2c[s2c$TF=="tet-CTCF", "TF"] <- "tetCTCF"
s2c$TF <- as.factor(s2c$TF)
BL_stallion <- colpal[levels(s2c$TF), "color"]

# generate transcript id to gene id conversion map
# convert transcript id to gene name

# import all kallisto count data
t2g <- read_tsv("../../data/rna_resources/t2g_synKLF1XBP1FOXO1.tsv") # a transcript id to gene id gene name dictionary, including a synthetic id for the codon optimized KLF1/XBP1/FOXO1 sequence
t2g$target_id <- lapply(t2g$target_id, function(n){strsplit(n, split="\\.")[[1]][1]}) %>% unlist
raw_data <- tximport(s2c$rnaPath, type="kallisto",tx2gene=t2g, ignoreTxVersion=TRUE)
colnames(raw_data$counts) <- s2c$sampleName
colnames(raw_data$abundance) <- s2c$sampleName

write.table(raw_data$abundance, paste0(outdir, "/raw_TPM_abundance.tsv"), sep="\t")
write.table(raw_data$counts, paste0(outdir, "/raw_counts.tsv"), sep="\t")

# plot PCA ---
#data <- round(as.data.frame(raw_data$counts))
data <- round(as.data.frame(fread(paste0(outdir, "/raw_counts.tsv"))) %>% column_to_rownames("V1"))
dds <- DESeqDataSetFromMatrix(countData = data, colData = s2c, design = ~ TF_Dose)
#dds <- DESeqDataSetFromMatrix(countData = data, colData = s2c, design = ~ TF)
dds <- vst(dds)


rv <- rowVars(assay(dds))
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
pca <- prcomp(t(assay(dds)[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], name = colnames(dds))


p1 <- muRtools::getDimRedPlot(d, annot=s2c, colorCol="TF", shapeCol="Plasmid_Dose",ptSize=5, colScheme=BL_stallion) + 
                theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)) + xlab(paste0("PC1 (", round(percentVar[1]*100),"%)")) + ylab(paste0("PC2 (", round(percentVar[2]*100),"%)"))

p2 <- muRtools::getDimRedPlot(d, annot=s2c["TF"], colorCol="TF", ptSize=5, colScheme=BL_stallion) + 
                theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)) + xlab(paste0("PC1 (", round(percentVar[1]*100),"%)")) + ylab(paste0("PC2 (", round(percentVar[2]*100),"%)"))

p3 <- muRtools::getDimRedPlot(d, annot=s2c["Plasmid_Dose"], colorCol="Plasmid_Dose", ptSize=5) + 
                theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)) + xlab(paste0("PC1 (", round(percentVar[1]*100),"%)")) + ylab(paste0("PC2 (", round(percentVar[2]*100),"%)"))

p4 <- muRtools::getDimRedPlot(d, annot=s2c["PlateNum"], colorCol="PlateNum", ptSize=5) + 
                theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)) + xlab(paste0("PC1 (", round(percentVar[1]*100),"%)")) + ylab(paste0("PC2 (", round(percentVar[2]*100),"%)"))

p5 <- muRtools::getDimRedPlot(d, annot=s2c["CellType"], colorCol="CellType", ptSize=5) + 
                theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)) + xlab(paste0("PC1 (", round(percentVar[1]*100),"%)")) + ylab(paste0("PC2 (", round(percentVar[2]*100),"%)"))

pdf(paste0(plotdir, "/PCA_all.pdf"), width=10, height=8)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
dev.off()

# plot PCA for each cell type
for (celltype in c("HEK293T")){
  subs2c <- s2c[s2c$CellType==celltype,]
  subs2c$TF <- as.factor(subs2c$TF)
  BL_stallion <- colpal[levels(subs2c$TF), "color"]
  subdata <- data[,subs2c$sampleName]

  dds <- DESeqDataSetFromMatrix(countData = subdata, colData = subs2c, design = ~ TF_Dose)
  dds <- vst(dds)
  suffix <- ""

  rv <- rowVars(assay(dds))
  select <- order(rv, decreasing = TRUE)[seq_len(min(500,
          length(rv)))]
  pca <- prcomp(t(assay(dds)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], name = colnames(dds))


  p1 <- muRtools::getDimRedPlot(d, annot=subs2c, colorCol="TF", shapeCol="Plasmid_Dose",ptSize=5, colScheme=BL_stallion) + 
                  theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)) + xlab(paste0("PC1 (", round(percentVar[1]*100),"%)")) + ylab(paste0("PC2 (", round(percentVar[2]*100),"%)"))

  p2 <- muRtools::getDimRedPlot(d, annot=subs2c["TF"], colorCol="TF", ptSize=5, colScheme=BL_stallion) + 
                  theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)) + xlab(paste0("PC1 (", round(percentVar[1]*100),"%)")) + ylab(paste0("PC2 (", round(percentVar[2]*100),"%)"))

  p3 <- muRtools::getDimRedPlot(d, annot=subs2c["Plasmid_Dose"], colorCol="Plasmid_Dose", ptSize=5) + 
                  theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)) + xlab(paste0("PC1 (", round(percentVar[1]*100),"%)")) + ylab(paste0("PC2 (", round(percentVar[2]*100),"%)"))

  p4 <- muRtools::getDimRedPlot(d, annot=subs2c["PlateNum"], colorCol="PlateNum", ptSize=5) + 
                  theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)) + xlab(paste0("PC1 (", round(percentVar[1]*100),"%)")) + ylab(paste0("PC2 (", round(percentVar[2]*100),"%)"))

  p5 <- getDimRedPlot_BL(d, annot=subs2c, colorCol="TF", alphaCol="Plasmid_Dose", colScheme=BL_stallion, ptSize=5) + 
    theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)) + xlab(paste0("PC1 (", round(percentVar[1]*100),"%)")) + ylab(paste0("PC2 (", round(percentVar[2]*100),"%)"))
  
  p6 <- getDimRedPlot_BL(d, annot=subs2c, colorCol="TF", sizeCol="Plasmid_Dose", colScheme=BL_stallion) + 
    theme_classic() + theme(text=element_text(size=20),panel.border = element_rect(color ="black",fill = NA, size = 1)) + xlab(paste0("PC1 (", round(percentVar[1]*100),"%)")) + ylab(paste0("PC2 (", round(percentVar[2]*100),"%)"))
  
  
  pdf(paste0(plotdir, "/PCA_", celltype, suffix, ".pdf"), width=10, height=8)
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  print(p6)
  dev.off()
}


###################### differential analysis (single deseq object) ######################
# add gene description to gene annotation
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                          dataset = "hsapiens_gene_ensembl",
                          host = 'https://ensembl.org')
id2name <- biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name","description"),
                          mart = mart, 
                          filters = "ensembl_gene_id",
                          values = rownames(data),
                          useCache = FALSE)
rownames(id2name) <- id2name$ensembl_gene_id
id2name["ENSG90000105610",] <- c("ENSG90000105610", "KLF1-syn", "KLF1 expressed from codon-optimized synthetic plasmid sequence")
id2name["ENSG90000100219",] <- c("ENSG90000100219", "XBP1-syn", "XBP1 expressed from codon-optimized synthetic plasmid sequence")
id2name["ENSG90000150907",] <- c("ENSG90000150907", "FOXO1-syn", "FOXO1 expressed from codon-optimized synthetic plasmid sequence")
write.table(id2name, paste0(outdir, "/id2name.tsv"), sep="\t")
id2name <- read.table(paste0(outdir, "/id2name.tsv"), sep="\t")

# differential 
for (celltype in c("HEK293T")){
  subs2c <- s2c[s2c$CellType==celltype,]
  subdata <- data[,subs2c$sampleName]

  cond2 <- "GFP" # denominator
  log2fc.cutoff <- 0.58
  padj.cutoff <- 0.05
  
  cond1_ls <- unique(subs2c$TF)
  cond1_ls <- cond1_ls[!cond1_ls==cond2] 
  
  subs2c$TF <- factor(subs2c$TF, levels=c(cond2, cond1_ls))
  
  dds <- DESeqDataSetFromMatrix(countData = subdata, colData = subs2c, design = ~ TF)
  suffix <- ""
    
  dds <- DESeq(dds)
  
  dir.create(paste0(outdir, "/tf_group/", celltype), recursive=T, showWarnings=F)
  #saveRDS(dds, paste0(outdir, "/tf_group/", celltype, "/dds_", celltype, suffix, ".rds"))
  #dds <- readRDS(paste0(outdir, "/tf_group/", celltype, "/dds_", celltype, suffix, ".rds"))
  for (cond1 in cond1_ls){
    cond1 <- gsub("-", ".", cond1) # deseq automatically converts hyphens to dots
    p <- plotDDS(dds, cond1, cond2, id2name, "TF", log2fc.cutoff=log2fc.cutoff, padj.cutoff=padj.cutoff, 
                  labels=c(cond1, "SPI1", "KLF1-syn", "FOXO1-syn", "XBP1-syn"), outdir=paste0(outdir, "/tf_group/", celltype, suffix))
    dir.create(paste0(plotdir, "/tf_group/", celltype), recursive=T, showWarnings=F)
    pdf(paste0(plotdir, "/tf_group/", celltype, "/MA_Volcano_plot_",cond1, "_vs_", cond2,"_log2fc",log2fc.cutoff, "_p", padj.cutoff, suffix, ".pdf"))
    print(p[[1]])
    print(p[[2]])
    dev.off()
  }
}



###################### differential analysis (each TF dose vs control) ######################
cond2 <- "GFP_1" # denominator
log2fc.cutoff <- 0.58
padj.cutoff <- 0.05

for (celltype in c("HEK293T")){
  for (plate in c(1,2,3)){
    message(paste0("plate ", plate))
    subs2c <- s2c[(s2c$CellType==celltype) & (s2c$PlateNum==plate),]
    subdata <- data[,subs2c$sampleName]

    cond1_ls <- unique(subs2c$TF_Dose)
    cond1_ls <- cond1_ls[cond1_ls!=cond2]
    subs2c$TF_Dose <- factor(subs2c$TF_Dose, levels=c(cond2, cond1_ls))
     
    out <- paste0(outdir, "/tf_dose/", celltype)
    plot <- paste0(plotdir, "/tf_dose/", celltype)
    dir.create(out, recursive=T, showWarnings=F)
    dir.create(plot, recursive=T, showWarnings=F)

    dds_obj <- paste0(out, "/dds_", celltype, "_plate", plate, ".rds")
    if (!file.exists(dds_obj)){
      dds <- DESeqDataSetFromMatrix(countData = subdata, colData = subs2c, design = ~TF_Dose)
      dds <- DESeq(dds)
      saveRDS(dds, dds_obj)
    } else{
      #dds <- readRDS(dds_obj)
      message("dds object exists. Skip DESeq.")
    }
    
    for (cond1 in cond1_ls){
      cond1 <- gsub("-", ".", cond1) # deseq automatically converts hyphens to dots

      labels_array_string <- sprintf("c(\"%s\",\"KLF1-syn\",\"FOXO1-syn\",\"XBP1-syn\")", strsplit(cond1, split="_")[[1]][1])
      command <- sprintf("sbatch -p wjg,biochem,sfgf --mem-per-cpu=100g --time=2:00:00 --job-name=%s --out=logs/de/slurm-%%j-%s --wrap \"Rscript scripts/01b_job_de.R %s %s %s %s %s %s %s %s %s %s %s %s\"",
                          paste0(celltype, "_", cond1, "_vs_", cond2), paste0(celltype, "_", cond1, "_vs_", cond2),
                          dds_obj, paste0(outdir, "/id2name.tsv"), cond1, cond2, "TF_Dose", 0.58, 0.05, out, plot, 
                          labels_array_string, "T", 50)
      
      print(command)
      system(command)
    }
  }
}